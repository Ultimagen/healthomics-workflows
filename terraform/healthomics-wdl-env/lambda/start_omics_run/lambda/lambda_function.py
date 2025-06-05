import logging
import os
import boto3
from botocore.exceptions import ClientError
from botocore.config import Config
import json

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

RUN_NAME_PARAM_KEY = 'base_file_name'
CACHE_SUFFIX = "_CACHE"
DYNAMIC_STORAGE = "DYNAMIC"
OMICS_DDB_VERSION_KEY = "version"
OMICS_DDB_WORKFLOW_KEY = "workflow"
OMICS_DDB_WORKFLOW_ID_KEY = "workflow_id"


def extract_run_name_from_inputs(input_params, workflow_name):
    return input_params.get(f'{workflow_name}.{RUN_NAME_PARAM_KEY}', "None")


def get_run_cache(omics_client, workflow_name, cache_behavior, cache_s3_bucket, run_tags: dict):
    cache_name = f"{workflow_name}_{cache_behavior}"
    all_run_caches = omics_client.list_run_caches()["items"]
    logger.debug(all_run_caches)
    if all_run_caches:
        existing_caches = [cache for cache in all_run_caches if cache["name"] == cache_name]
        if existing_caches:
            cache_id = existing_caches[0]["id"]
            return cache_id
    response = omics_client.create_run_cache(
        name=cache_name,
        cacheS3Location=f"s3://{cache_s3_bucket}/{cache_behavior.lower()}/{workflow_name}/",
        cacheBehavior=cache_behavior,
        tags=run_tags
    )
    cache_id = response["id"]
    logger.info(f"Cache {cache_name} was created (id: {cache_id})")
    return cache_id


def start_run(omics_client, **kwargs):
    start_run_resp = omics_client.start_run(**kwargs)
    logger.info(start_run_resp)
    omics_run_id = start_run_resp["id"]
    omics_run_arn = start_run_resp["arn"]
    # This additional tagging must be done after the run submitted, since we don't have the omics run id is taken from start run response
    resp = omics_client.tag_resource(
        resourceArn=omics_run_arn,
        tags={
            'omics_run_id': omics_run_id
        }
    )
    logger.info(f"tag_resources response: {resp}")
    return omics_run_id


def _get_client(server_name, max_attempts=20):
    config = Config(retries=dict(total_max_attempts=max_attempts))

    client = boto3.client(server_name, config=config)
    return client


def _get_resource(server_name):
    resource = boto3.resource(server_name)
    return resource


def _get_workflow_from_dynamodb(version: str, workflow_name: str):
    client = _get_resource("dynamodb")
    table_name = os.environ['OMICS_WORKFLOW_DDB_TABLE']
    table = client.Table(table_name)
    key = {OMICS_DDB_VERSION_KEY: version, OMICS_DDB_WORKFLOW_KEY: workflow_name}
    logger.info(f"Get item from {table_name} table with key {key}")
    db_item = table.get_item(Key=key).get("Item")
    logging.debug(f"Item found:\n{db_item}")
    if not db_item:
        raise RuntimeError(
            f"workflow {workflow_name} version: {version} is missing from dynamodb. Make sure the workflow is deployed")
    return db_item[OMICS_DDB_WORKFLOW_ID_KEY]


def lambda_handler(event, context):
    logger.info(f"Event:\n{event}")
    workflow_name = event['WorkflowName']
    workflow_version = event['WorkflowVersion']
    run_id = event['RunId']
    input_params = event['InputParams']
    run_name = event.get('RunName', extract_run_name_from_inputs(input_params, workflow_name))
    # is_long_run = event.get('LongRun', False)
    # run_group_id = os.environ['LONG_RUN_GROUP_ID'] if is_long_run else os.environ['STANDARD_RUN_GROUP_ID']
    cache_behavior = event.get("CacheBehavior", "CACHE_ON_FAILURE")

    request_id = context.aws_request_id
    role_arn = os.environ['OMICS_ROLE_ARN']
    output_s3_path = f"s3://{os.environ['OMICS_OUTPUTS_BUCKET']}/{workflow_name}"
    cache_s3_bucket = os.environ['OMICS_CACHE_BUCKET']

    workflow_id = _get_workflow_from_dynamodb(workflow_version, workflow_name)

    omics_client = _get_client("omics")
    try:
        logger.debug("Attempt to start workflow run")

        run_tags = {
            "workflow_name": workflow_name,
            "run_name": run_name
        }
        call_cache_id = get_run_cache(omics_client, workflow_name, cache_behavior, cache_s3_bucket, run_tags)
        start_run_args = {
            'workflowId': workflow_id,
            "name": run_id,
            "roleArn": role_arn,
            "parameters": input_params,
            "outputUri": output_s3_path,
            "tags": run_tags,
            "requestId": request_id,
            # "runGroupId": run_group_id,
            # If the quota for maximum runs has been met, the earliest runs with REMOVE retention mode are deleted first
            "retentionMode": "REMOVE"
        }
        if call_cache_id:
            start_run_args.update({"cacheId": call_cache_id})
            logger.info(f"Run will use cache id: {call_cache_id}. with cacheBehavior: {cache_behavior}")

        storage_type = event.get("storageType")
        if storage_type:
            start_run_args.update({"storageType": storage_type})

        storage_capacity = event.get("storageCapacity")
        if storage_capacity:
            start_run_args.update({"storageCapacity": storage_capacity})

        omics_run_id = start_run(omics_client, **start_run_args)
    except ClientError as e:
        raise Exception("boto3 client error : " + e.__str__())
    except Exception as e:
        raise Exception("Unexpected error : " + e.__str__())
    return {"OmicsRunId": omics_run_id}
