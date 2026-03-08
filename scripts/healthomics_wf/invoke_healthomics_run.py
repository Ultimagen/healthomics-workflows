import json
import argparse
import boto3

LAMBDA_NAME = "startOmicsRun"
OMICS_WORKFLOW_DDB_TABLE = "OmicsWorkflows"


def _get_client(service_name, aws_region, aws_profile):
    if aws_profile:
        boto3.setup_default_session(profile_name=aws_profile)
    return boto3.client(service_name, region_name=aws_region)


def _get_resource(service_name, aws_region, aws_profile):
    if aws_profile:
        boto3.setup_default_session(profile_name=aws_profile)
    return boto3.resource(service_name, region_name=aws_region)


def _parse_version(version_str):
    """
    Parse a semantic version string (e.g., "1.28.0") into a tuple of integers for comparison.
    """
    try:
        return tuple(int(part) for part in version_str.split("."))
    except ValueError:
        # If parsing fails, return the original string for fallback lexicographic comparison
        return (version_str,)


def get_latest_workflow_version(workflow_name, aws_region, aws_profile=None):
    """
    Query DynamoDB to find the latest version for a given workflow.
    Returns the version string with the highest semantic version.
    """
    dynamodb = _get_resource("dynamodb", aws_region, aws_profile)
    table = dynamodb.Table(OMICS_WORKFLOW_DDB_TABLE)

    response = table.scan(
        FilterExpression=boto3.dynamodb.conditions.Attr("workflow").eq(workflow_name)
    )

    items = response.get("Items", [])
    if not items:
        raise ValueError(
            f"No versions found for workflow '{workflow_name}' in DynamoDB table '{OMICS_WORKFLOW_DDB_TABLE}'. "
            f"Please deploy the workflow first using create_healthomics_workflow.py with --use-dynamodb flag."
        )

    # Sort versions using semantic versioning comparison
    versions = [item["version"] for item in items]
    versions_sorted = sorted(versions, key=_parse_version, reverse=True)
    latest_version = versions_sorted[0]

    print(f"Found {len(versions)} version(s) for workflow '{workflow_name}': {', '.join(versions_sorted)}")
    print(f"Using latest version: {latest_version}")

    return latest_version


def invoke_lambda(workflow_name, workflow_version, run_id, input_params_file, aws_region, aws_profile=None,
                  long_run=False, cache_behavior=None):
    # Read the JSON input from file
    with open(input_params_file, 'r') as f:
        input_params = json.load(f)

    # Create payload
    payload = {
        "WorkflowName": workflow_name,
        "WorkflowVersion": workflow_version,
        "RunId": run_id,
        "InputParams": input_params
    }
    if long_run:
        payload.update({"LongRun": True})
    if cache_behavior:
        payload.update({"CacheBehavior": cache_behavior})
    # Create a Lambda client
    client = _get_client('lambda', aws_region, aws_profile)

    # Invoke the Lambda function
    response = client.invoke(
        FunctionName=LAMBDA_NAME,
        InvocationType='RequestResponse',  # Sync
        Payload=json.dumps(payload)
    )

    # Read and print the response
    response_payload = response['Payload'].read().decode('utf-8')
    print(f"Lambda response: {response_payload}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Invoke a HealthOmics start run Lambda function with a JSON file input.")
    parser.add_argument("--omics-workflow-name", required=True, help="Workflow name")
    parser.add_argument("--workflow-version", required=False,
                        help="Workflow version. If not provided, the latest version from DynamoDB will be used.")
    parser.add_argument("--run-id", required=True, help="Run ID")
    parser.add_argument("--input-params-file", required=True, help="Path to the JSON file for InputParams")
    parser.add_argument("--long-run", action='store_true',
                        help="Mark the run as long and run it under 'long' run group. By default will run under 'standard' run group")
    parser.add_argument("--cache-behavior", help="HealthOmics cache behavior")
    parser.add_argument("--aws-region", help="AWS region", default="us-east-1")
    parser.add_argument("--aws-profile", help="AWS CLI profile", required=False)

    args = parser.parse_args()

    # Get workflow version - use provided value or fetch latest from DynamoDB
    workflow_version = args.workflow_version
    if not workflow_version:
        workflow_version = get_latest_workflow_version(
            workflow_name=args.omics_workflow_name,
            aws_region=args.aws_region,
            aws_profile=args.aws_profile
        )

    invoke_lambda(
        workflow_name=args.omics_workflow_name,
        workflow_version=workflow_version,
        run_id=args.run_id,
        input_params_file=args.input_params_file,
        aws_region=args.aws_region,
        aws_profile=args.aws_profile,
        long_run=args.long_run,
        cache_behavior=args.cache_behavior
    )
