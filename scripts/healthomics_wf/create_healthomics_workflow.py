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


def create_workflow(aws_region, omics_wf_name, wf_root, workflow_name):
    create_omics_workflow(aws_region, omics_wf_name, wf_root, workflow_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create AWS healthomics workflow")
    parser.add_argument("workflow", help="workflow to create (folder under workflows)")
    parser.add_argument("aws_region", help="AWS region")
    parser.add_argument("s3_bucket",
                        help="bucket name to copy resources files to. This bucket must be accessed by the service "
                             "role that will be used to run the workflow.")
    parser.add_argument("--omics_workflow_name",
                        help="how to name the generated omics workflow, if empty will use 'workflow' arg",
                        required=False)
    parser.add_argument("--input_template",
                        help="input template json file name to localize (if empty will localize all input templates)",
                        required=False)
    parser.add_argument("--aws_profile", help="AWS CLI profile", required=False)
    args = parser.parse_args()

    workflow_root = f'{Path(__file__).resolve().parent.parent.parent}/workflows/{args.workflow}'
    localize_workflow(workflow_root, args.aws_region, args.s3_bucket, args.aws_profile,
                      args.input_template)

    omics_workflow_name = args.omics_workflow_name or args.workflow
    create_workflow(args.aws_region, omics_workflow_name, workflow_root, args.workflow)
