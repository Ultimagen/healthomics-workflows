import argparse
import json
import os
from datetime import datetime, timezone
import tarfile

import boto3
import logging
import datetime

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)

PIPELINE_VERSION_TAG = "pipeline_version"
UTC = timezone.utc

class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime.datetime):
            return {
                '__type__': 'datetime',
                'year': obj.year,
                'month': obj.month,
                'day': obj.day,
                'hour': obj.hour,
                'minute': obj.minute,
                'second': obj.second,
                'microsecond': obj.microsecond,
            }
        elif not isinstance(obj, (str, int, float, list, tuple, dict, bool, type(None))):
            return str(obj)
        else:
            return json.JSONEncoder.default(self, obj)


def get_run_info(run_id, client):
    omics_run = client.get_run(id=run_id)
    logging.debug(omics_run)
    omics_run.update({"duration": (omics_run['stopTime'] if 'stopTime' in omics_run else datetime.now(UTC)) - omics_run[
        'startTime']})

    response = client.list_run_tasks(id=run_id)
    tasks = response['items']
    while response.get('nextToken'):
        response = client.list_run_tasks(id=run_id, startingToken=response.get('nextToken'))
        tasks += response['items']

    def calc_task_duration(task):
        stop_time = task['stopTime'] if 'stopTime' in task else datetime.now(UTC)
        start_time = task[
            'startTime'] if 'startTime' in task else stop_time  # when task is CANCELLED sometime there's no startTime
        return stop_time - start_time

    tasks = [{**task, "duration": calc_task_duration(task)} for task in tasks]

    del omics_run['ResponseMetadata']

    omics_run['tasks'] = tasks
    logging.debug(omics_run)
    return omics_run


def get_workflow_metadata(client, wf_id):
    response = client.get_workflow(id=wf_id)
    return response


def extract_pipeline_version(wf_metadata):
    if not wf_metadata.get('tags'):
        return ""
    return wf_metadata['tags'].get('pipeline_version')


FAILED_STATUS = "FAILED"
OMICS_LOG_GROUP = '/aws/omics/WorkflowLog'


def get_log_for_task(run_id, run, task_id=None, output_prefix='', failed: bool = True):
    # If no given a specific task_id, get the logs for all tasks or failed tasks
    if not task_id:
        task_ids = [task['taskId'] for task in run['tasks'] if (not failed or task["status"] == FAILED_STATUS)]
        task_names = {task['taskId']: task['name'] for task in run['tasks'] if
                      (not failed or task["status"] == FAILED_STATUS)}

    else:
        task_ids = [task_id]
        task_names = {task['taskId']: task['name'] for task in run['tasks'] if task['taskId'] == task_id}

    # Get the logs for each task
    for task_id in task_ids:
        task_name = task_names[task_id]
        logging.info("------------------------------------------")
        logging.info(f"Getting log for task {task_name} (taskId: {task_id})")
        log_stream_name = f'run/{run_id}/task/{task_id}'

        output = f'{output_prefix}run_{run_id}_task_{task_id}_{task_name}.log'
        fetch_save_log(log_stream_name, output, session)

    # in case that the run failed but there're no failed tasks the run's engine log should include the error
    if run["status"] == FAILED_STATUS and FAILED_STATUS not in [task['status'] for task in run['tasks']]:
        logging.info(f"Run status is {FAILED_STATUS} but no failed tasks. Getting run's engine log")
        log_stream_name = f'run/{run_id}/engine'
        output = f'{output_prefix}run_{run_id}_engine.log'
        fetch_save_log(log_stream_name, output, session)


def fetch_save_log(log_stream_name, output, session=None):
    client = get_aws_client('logs', session)

    logging.info(f"Getting log events for log group '{OMICS_LOG_GROUP}' and log stream '{log_stream_name}'")

    # get first page of log events
    response = client.get_log_events(
        logGroupName=OMICS_LOG_GROUP,
        logStreamName=log_stream_name,
        startFromHead=True
    )

    # check if events is not empty
    if not response.get('events'):
        logging.info(f"No events found for log group '{OMICS_LOG_GROUP}' and log stream '{log_stream_name}'")
        return

    # write event message to output file
    output = adjust_file_name(output)
    with open(output, 'w') as file:
        for event in response['events']:
            file.write(f"{event['message']}\n")

        # get next page of log events
        while len(response.get('events')) > 0 and response.get('nextForwardToken'):
            response = client.get_log_events(
                logGroupName=OMICS_LOG_GROUP,
                logStreamName=log_stream_name,
                nextToken=response.get('nextForwardToken'),
                startFromHead=True
            )
            for event in response['events']:
                file.write(f"{event['message']}\n")

    logging.info(f"Log file saved to: {output}")
    files_to_add.append(output)


def json_to_file(output_file, json_data):
    output_file = adjust_file_name(output_file)
    logging.info(f"Saving run info to: {output_file}")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(json_data, f, indent=4, cls=DateTimeEncoder)
    files_to_add.append(output_file)


def get_aws_client(service, session):
    if session:
        client = session.client(service, args.aws_region)
    else:
        client = boto3.client(service, args.aws_region)
    return client


def create_tar(tar_file):
    with tarfile.open(tar_file, "w") as tar:
        for f in files_to_add:
            tar.add(f)
    logging.info(f"tar file saved to: {tar_file}")


def adjust_file_name(file_name):
    if output_dir:
        file_name = f"{output_dir}/{file_name}"
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract information and logs from AWS HealthOmics run to ease failures debugging")
    parser.add_argument("run_id", help="AWS HealthOmics run id")
    parser.add_argument("--aws_region", help="AWS Region (default 'us-east-1')", required=False)
    parser.add_argument('--task-id', type=str,
                        help="HealthOmics workflow task-id to analyze. Leave empty to get the logs for all tasks",
                        default=None)
    parser.add_argument('--failed', dest='failed', action='store_true')
    parser.add_argument('--no-failed', dest='failed', action='store_false',
                        help="Download logs for all run's tasks (by default only failed tasks' log will be downloaded)")
    parser.set_defaults(failed=True)
    parser.add_argument('--output', type=str, help="Output dir to save out files", default=None)
    parser.add_argument('--output-prefix', type=str, help="File name prefix for the output files", required=False,
                        default='')
    parser.add_argument('--tar', dest='tar', action='store_true')
    parser.add_argument('--no-tar', dest='tar', action='store_false',
                        help="Don't tar all generated files")
    parser.set_defaults(tar=True)
    args = parser.parse_args()

    output_prefix = args.output_prefix
    output_dir = args.output
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    run_id = args.run_id
    region = args.aws_region
    files_to_add = []

    session = boto3.Session()

    omics_client = get_aws_client('omics', session)

    run = get_run_info(run_id, omics_client)

    json_to_file(f"{output_prefix}run_{run_id}_info.json", run)

    workflow_id = run['workflowId']

    workflow_metadata = get_workflow_metadata(omics_client, workflow_id)

    json_to_file(f"{output_prefix}run_{run_id}_workflow_info.json", run)
    pipeline_version = extract_pipeline_version(workflow_metadata)

    logging.info(f"workflow version: {pipeline_version}")
    workflow_version_file = adjust_file_name(f"{output_prefix}run_workflow_version.txt")
    with open(workflow_version_file, 'w') as file:
        file.write(pipeline_version)

    get_log_for_task(run_id, run, args.task_id, args.output_prefix, args.failed)

    if args.tar:
        tar_file = adjust_file_name(f"{output_prefix}run.tar")
        create_tar(tar_file)
