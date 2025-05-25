import json
import re
import shutil
import time
import glob
from os import path
from os.path import exists
from pathlib import Path
from shutil import make_archive, copy2
from tempfile import TemporaryDirectory
import boto3
from botocore.exceptions import ClientError
import logging

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)

PIPELINE_VERSION_TAG = "pipeline_version"
PARAMS_DEF_SUFFIX = "params_def.json"
OMICS_WF_ACTIVE = "ACTIVE"
OMICS_WF_DONE_STATUSES = [OMICS_WF_ACTIVE, "DELETED", "FAILED"]
OMICS_PRIVATE_TYPE = "PRIVATE"
TASKS_SUBDIR = "tasks/"
WORKFLOW_ENGINE = "WDL"
WORKFLOW_STORAGE_TYPE = "DYNAMIC"
OMICS_WORKFLOW_DDB_TABLE = "OmicsWorkflows"
OMICS_DDB_VERSION_KEY = "version"
OMICS_DDB_WORKFLOW_KEY = "workflow"
OMICS_DDB_WORKFLOW_ID_KEY = "workflow_id"


def extract_pipeline_version(wdl_file_path):
    with open(wdl_file_path, 'r') as file:
        wdl_content = file.read()
        # Regex pattern to match the pipeline_version line
        pattern = re.compile(r'String\s+pipeline_version\s*=\s*"([^"]+)"')

        match = pattern.search(wdl_content)
        if match:
            return match.group(1)
    return None


def safe_copy(source: Path, target: Path, dryrun: bool = False):
    assert exists(source), f"Missing file {source}"
    if not dryrun:
        if source.is_file():
            logging.info(f"copy {source} -> {target}")
            shutil.copy(source.absolute(), target.absolute())


def zip_workflow_files(workflow_name, workflow_root, tmpdir):
    logging.info(f"Make zip from wdl files in {workflow_root}")
    wdl_files = glob.glob(f"{workflow_root}/**/*.wdl", recursive=True)
    for wdl in wdl_files:
        subdir = Path(f'{tmpdir}/{workflow_name}/') if TASKS_SUBDIR not in wdl else Path(
            f'{tmpdir}/{workflow_name}/{TASKS_SUBDIR}')
        logging.debug(f"Copy {wdl} to {subdir}")
        subdir.mkdir(parents=True, exist_ok=True)
        copy2(wdl, subdir)

    zip_file = path.join(tmpdir, workflow_name)
    zip_path = Path(make_archive(zip_file, "zip", zip_file))
    return zip_path


def create_omics_workflow(aws_region, omics_workflow_name, workflow_root, workflow, aws_profile=None,
                          use_dynamodb=False):
    main_wdl = f"{workflow}.wdl"
    workflow_version = extract_pipeline_version(f"{workflow_root}/{main_wdl}")
    workflow_id = None
    if use_dynamodb:
        workflow_db_item = _get_workflow_from_dynamodb(workflow_version, omics_workflow_name, aws_region, aws_profile)
        if workflow_db_item:
            workflow_id = workflow_db_item[OMICS_DDB_WORKFLOW_ID_KEY]
            logging.info(f"Workflow found in dynamodb, will use existing workflow_id: {workflow_id}")
    if not workflow_id:
        workflow_id = _create_workflow(aws_profile, aws_region, main_wdl, omics_workflow_name, workflow, workflow_root,
                                       workflow_version)
        if use_dynamodb:
            _write_workflow_in_dynamodb(workflow_id, workflow_version, omics_workflow_name, aws_region, aws_profile)
    return workflow_id


def _create_workflow(aws_profile, aws_region, main_wdl, omics_workflow_name, workflow, workflow_root, workflow_version):
    logging.info(f"Create omics workflow for {workflow}, version {workflow_version}, main wdl: {main_wdl}")
    params_def_file = Path(f"{workflow_root}/{workflow}_{PARAMS_DEF_SUFFIX}")
    with params_def_file.open(encoding="UTF-8") as source:
        wdl_params = json.load(source)
    with TemporaryDirectory() as tmpdir:
        zip_path = zip_workflow_files(workflow, workflow_root, tmpdir)
        logging.debug(f"Read zip as bytes stream")
        with open(zip_path, "rb") as f:
            zip_bytes = f.read()

            omics_client = _get_aws_client("omics", aws_region, aws_profile)
            workflow_id = _create_omics_workflow(main_wdl, omics_client, omics_workflow_name, wdl_params,
                                                 workflow_version, zip_bytes)
            # _create_omics_workflow_version(main_wdl, omics_client, workflow_id, omics_workflow_name, wdl_params,
            #                                      workflow_version, zip_bytes)
    return workflow_id


def _create_omics_workflow(main_wdl, omics_client, omics_workflow_name, wdl_params, workflow_version, zip_bytes):
    try:
        response = omics_client.create_workflow(
            name=omics_workflow_name,
            engine=WORKFLOW_ENGINE,
            definitionZip=zip_bytes,
            main=main_wdl,
            parameterTemplate=wdl_params,
            storageType=WORKFLOW_STORAGE_TYPE,
            tags={
                PIPELINE_VERSION_TAG: workflow_version
            }
        )
        logging.debug(response)
        workflow_id = response["id"]
        wf_status = response["status"]
        while wf_status not in OMICS_WF_DONE_STATUSES:
            time.sleep(5)
            response = omics_client.get_workflow(id=workflow_id)
            wf_status = response["status"]
        if wf_status != OMICS_WF_ACTIVE:
            logging.error(
                f"{omics_workflow_name} omics workflow creation failed with msg: {response['statusMessage']}")
            exit(1)
        logging.info(f"{omics_workflow_name} omics workflow created successfully. workflowId: {workflow_id}")

    except ClientError as e:
        logging.error(e)
        exit(1)
    return workflow_id


# def _create_omics_workflow_version(main_wdl, omics_client, workflow_id, omics_workflow_name, wdl_params, workflow_version, zip_bytes):
#     try:
#         response = omics_client.create_workflow_version(
#             workflowId=workflow_id,
#             versionName=workflow_version,
#             definitionZip=zip_bytes,
#             engine=WORKFLOW_ENGINE,
#             main=main_wdl,
#             parameterTemplate=wdl_params,
#             storageType=WORKFLOW_STORAGE_TYPE,
#             tags={
#                 PIPELINE_VERSION_TAG: workflow_version
#             }
#         )
#         logging.debug(response)
#         version_arn = response["arn"]
#         wf_status = response["status"]
#         while wf_status not in OMICS_WF_DONE_STATUSES:
#             time.sleep(5)
#             response = omics_client.get_workflow(id=workflow_id)
#             wf_status = response["status"]
#         if wf_status != OMICS_WF_ACTIVE:
#             logging.error(
#                 f"{omics_workflow_name} omics workflow version '{workflow_version}' creation failed with msg: {response['statusMessage']}")
#             exit(1)
#         logging.info(f"{omics_workflow_name} omics workflow version '{workflow_version}' created successfully. arn: {version_arn}")
#
#     except ClientError as e:
#         logging.error(e)
#         exit(1)

def _get_aws_resource(service_name, aws_region, aws_profile):
    if aws_profile:
        boto3.setup_default_session(profile_name=aws_profile)
    return boto3.resource(service_name, region_name=aws_region)


def _get_aws_client(service_name, aws_region, aws_profile):
    if aws_profile:
        boto3.setup_default_session(profile_name=aws_profile)
    return boto3.client(service_name, region_name=aws_region)


def _write_workflow_in_dynamodb(workflow_id: str,
                                version: str,
                                workflow_name: str,
                                aws_region: str,
                                aws_profile: str = None):
    client = _get_aws_resource("dynamodb", aws_region, aws_profile)
    table = client.Table(OMICS_WORKFLOW_DDB_TABLE)

    db_item = {
        OMICS_DDB_VERSION_KEY: version,
        OMICS_DDB_WORKFLOW_KEY: workflow_name,
        OMICS_DDB_WORKFLOW_ID_KEY: workflow_id
    }
    table.put_item(
        Item=db_item
    )
    logging.info(
        f"omics workflow {workflow_name} workflowId: {workflow_id}, version '{version}' was written to dynamodb table: {OMICS_WORKFLOW_DDB_TABLE}")
    return True


def _get_workflow_from_dynamodb(version: str, workflow_name: str, aws_region: str, aws_profile: str = None):
    client = _get_aws_resource("dynamodb", aws_region, aws_profile)
    table = client.Table(OMICS_WORKFLOW_DDB_TABLE)
    key = {OMICS_DDB_VERSION_KEY: version, OMICS_DDB_WORKFLOW_KEY: workflow_name}
    logging.info(f"Get item from {OMICS_WORKFLOW_DDB_TABLE} table with key {key}")
    db_item = table.get_item(Key=key).get("Item")
    logging.debug(f"Item found:\n{db_item}")
    return db_item
