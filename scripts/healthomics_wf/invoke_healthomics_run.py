import json
import argparse
import boto3

LAMBDA_NAME = "startOmicsRun"


def _get_client(service_name, aws_region, aws_profile):
    if aws_profile:
        boto3.setup_default_session(profile_name=aws_profile)
    return boto3.client(service_name, region_name=aws_region)


def invoke_lambda(workflow_name, workflow_version, run_id, input_params_file, aws_region, aws_profile=None,
                  run_group_id=None, cache_behavior=None):
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
    if run_group_id:
        payload.update({"RunGroupID": run_group_id})
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
    parser.add_argument("--workflow-name", required=True, help="Workflow name")
    parser.add_argument("--workflow-version", required=True, help="Workflow version")
    parser.add_argument("--run-id", required=True, help="Run ID")
    parser.add_argument("--input-params-file", required=True, help="Path to the JSON file for InputParams")
    parser.add_argument("--run-group-id", help="Run group id, by default will run under 'standard' run group")
    parser.add_argument("--cache-behavior", help="HealthOmics cache behavior")
    parser.add_argument("--input-params-file", required=True, help="Path to the JSON file for InputParams")
    parser.add_argument("--aws-region", help="AWS region", default="us-east-1")
    parser.add_argument("--aws-profile", help="AWS CLI profile", required=False)

    args = parser.parse_args()

    invoke_lambda(
        workflow_name=args.workflow_name,
        workflow_version=args.workflow_version,
        run_id=args.run_id,
        input_params_file=args.input_params_file,
        aws_region=args.aws_region,
        aws_profile=args.aws_profile
    )
