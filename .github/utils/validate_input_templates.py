import glob
import json

import boto3
import logging
from botocore import UNSIGNED
from botocore.config import Config

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)


def is_s3_path(value):
    return isinstance(value, str) and value.startswith('s3://') and '<' not in value


def s3_path_exists(s3_client, s3_path):
    try:
        bucket, key = s3_path.replace("s3://", "").split("/", 1)
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except Exception as e:
        return False


def validate_json_data(s3_client, data):
    all_paths_valid = True
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, (dict, list)):
                if not validate_json_data(s3_client, value):
                    all_paths_valid = False
            elif is_s3_path(value):
                if not s3_path_exists(s3_client, value):
                    logging.error(f"S3 path does not exist: {value}")
                    all_paths_valid = False
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, (dict, list)):
                if not validate_json_data(s3_client, item):
                    all_paths_valid = False
            elif is_s3_path(item):
                if not s3_path_exists(s3_client, item):
                    logging.error(f"S3 path does not exist: {item}")
                    all_paths_valid = False
    else:
        for key, value in data.items():
            if is_s3_path(value):
                if not s3_path_exists(s3_client, value):
                    logging.error(f"S3 path does not exist: {value}")
                    all_paths_valid = False
    return all_paths_valid


def validate_json_files(json_files):
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    all_paths_valid = True

    for json_file in json_files:
        with open(json_file, 'r') as file:
            data = json.load(file)
            if not validate_json_data(s3_client, data):
                all_paths_valid = False

    return all_paths_valid


if __name__ == "__main__":

    input_templates = glob.glob(f"workflows/*/input_templates/*.json", recursive=True)
    logging.info(f"{len(input_templates)} input templates found")
    if not validate_json_files(input_templates):
        logging.error("some input template contain non existing s3 paths")
        exit(1)
