import argparse
import logging
import glob
from pathlib import Path
from utils.modules.localize_dockers import localize_wdl_docker_images
from utils.modules.localize_s3 import localize_s3_files
from utils.modules.omics_workflows import create_omics_workflow

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)

GLOBALS_WDL = "tasks/globals.wdl"


def localize_workflow(wf_root, aws_region, s3_bucket, aws_profile=None, input_template=None):
    logging.info(f"Starting localize workflow: {args.workflow}")

    globals_wdl_file = f"{wf_root}/{GLOBALS_WDL}"

    localize_wdl_docker_images(globals_wdl_file, aws_region, aws_profile)

    if input_template:
        files_to_localize_s3 = [input_template]
    else:
        files_to_localize_s3 = glob.glob(f'{wf_root}/input_templates/*.json')
    files_to_localize_s3.append(globals_wdl_file)

    for input_file in files_to_localize_s3:
        logging.info(f"localizing {input_file}")
        localize_s3_files(input_file, s3_bucket)
    logging.info(f"workflow: {args.workflow} localization completed")


def localize_and_create(workflow_folder, aws_region, s3_bucket, omics_workflow_name, workflow_root=None,
                        aws_profile=None, input_template=None, use_dynamodb=False):
    if not workflow_root:
        workflow_root = f"{Path(__file__).resolve().parent.parent.parent}/workflows"
    workflow_root = f"{workflow_root}/{workflow_folder}"
    localize_workflow(workflow_root, aws_region, s3_bucket, aws_profile, input_template)
    create_omics_workflow(aws_region, omics_workflow_name, workflow_root, workflow_folder, aws_profile, use_dynamodb)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create AWS healthomics workflow")
    parser.add_argument("workflow", help="workflow to create (folder under workflows)")
    parser.add_argument("--s3-bucket",
                        help="bucket name to copy resources files to. This bucket must be accessed by the service "
                             "role that will be used to run the workflow.")
    parser.add_argument("--aws-region", help="AWS region")
    parser.add_argument("--workflow-root", help="workflows root folder",
                        required=False)
    parser.add_argument("--omics-workflow-name",
                        help="how to name the generated omics workflow, if empty will use 'workflow' arg",
                        required=False)
    parser.add_argument("--input-template",
                        help="input template json file name to localize (if empty will localize all input templates)",
                        required=False)
    parser.add_argument("--use-dynamodb", action='store_true',
                        help="whether to manage workflow versions in dynamodb")
    parser.add_argument("--aws-profile", help="AWS CLI profile", required=False)
    args = parser.parse_args()

    workflow_name = args.omics_workflow_name or args.workflow
    localize_and_create(
        workflow_folder=args.workflow, aws_region=args.aws_region, s3_bucket=args.s3_bucket,
        workflow_root=args.workflow_root,
        omics_workflow_name=workflow_name, aws_profile=args.aws_profile,
        input_template=args.input_template, use_dynamodb=args.use_dynamodb
    )
