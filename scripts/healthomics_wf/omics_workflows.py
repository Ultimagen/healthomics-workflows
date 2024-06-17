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


def zip_workflow_files(workflow_name, workflow_root):
    with TemporaryDirectory() as tmpdir:
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


def create_omics_workflow(aws_region, omics_workflow_name, workflow_root, workflow_name):
    main_wdl = f"{workflow_name}.wdl"
    workflow_version = extract_pipeline_version(f"{workflow_root}/{main_wdl}")
    logging.info(f"Create omics workflow for {workflow_name}, version {workflow_version}, main wdl: {main_wdl}")
    params_def_file = Path(f"{workflow_root}/{workflow_name}_{PARAMS_DEF_SUFFIX}")

    with params_def_file.open(encoding="UTF-8") as source:
        wdl_params = json.load(source)
    zip_path = zip_workflow_files(workflow_name, workflow_root)
    logging.debug(f"Read zip as bytes stream")
    with open(zip_path, "rb") as f:
        zip_bytes = f.read()

    omics_client = boto3.client("omics", aws_region)
    try:
        response = omics_client.create_workflow(
            name=omics_workflow_name,
            engine="WDL",
            definitionZip=zip_bytes,
            main=main_wdl,
            parameterTemplate=wdl_params,
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
            logging.error(f"{omics_workflow_name} omics workflow creation failed with msg: {response['statusMessage']}")
            exit(1)
        logging.info(f"{omics_workflow_name} omics workflow created successfully. workflowId: {workflow_id}")

    except ClientError as e:
        logging.error(e)
        exit(1)
