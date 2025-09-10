import argparse
import json
import re
import sys

import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError
import logging

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)


def file_exists_in_s3(s3_client, bucket, key):
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except s3_client.exceptions.ClientError as e:
        if e.response['Error']['Code'] == '404':
            return False
        else:
            raise


def copy_s3_file(s3_client, src_bucket, src_key, dest_bucket, dest_prefix=None):
    dest_key = f"{dest_prefix}/src_key" if dest_prefix else src_key  # keeping the same key in the destination bucket unless specified new dest_prefix
    if file_exists_in_s3(s3_client, dest_bucket, dest_key):
        logging.info(f"File already exists: {dest_key} in {dest_bucket}, skipping copy")
        return f"s3://{dest_bucket}/{dest_key}"

    try:
        copy_source = {'Bucket': src_bucket, 'Key': src_key}
        s3_client.copy(copy_source, dest_bucket, dest_key)
        logging.info(f"Copied {src_key} from {src_bucket} to {dest_bucket}")
        return f"s3://{dest_bucket}/{dest_key}"
    except NoCredentialsError:
        logging.error("No AWS credentials found.")
    except PartialCredentialsError:
        logging.error("Incomplete AWS credentials found.")
    except Exception as e:
        logging.error(f"Failed to copy {src_key} from {src_bucket} to {dest_bucket}: {e}")
    return f"s3://{src_bucket}/{src_key}"


def update_s3_paths(data, s3_client, dest_bucket, dest_prefix=None):
    if isinstance(data, dict):
        for key, value in data.items():
            data[key] = update_s3_paths(value, s3_client, dest_bucket, dest_prefix)
    elif isinstance(data, list):
        for i in range(len(data)):
            data[i] = update_s3_paths(data[i], s3_client, dest_bucket, dest_prefix)
    elif isinstance(data, str) and data.startswith('s3://') and '<' not in data:
        src_bucket, src_key = re.match(r's3://([^/]+)/(.+)', data).groups()
        return copy_s3_file(s3_client, src_bucket, src_key, dest_bucket, dest_prefix)
    return data


def parse_wdl_for_s3_paths(wdl_content):
    s3_paths = []
    pattern = re.compile(r'"(s3://[^"]+)"')
    for match in pattern.finditer(wdl_content):
        s3_paths.append(match.group(1))
    return s3_paths


def update_wdl_s3_paths(wdl_content, s3_client, dest_bucket, dest_prefix=None):
    s3_paths = parse_wdl_for_s3_paths(wdl_content)
    for path in s3_paths:
        src_bucket, src_key = re.match(r's3://([^/]+)/(.+)', path).groups()
        new_path = copy_s3_file(s3_client, src_bucket, src_key, dest_bucket, dest_prefix)
        wdl_content = wdl_content.replace(path, new_path)
    return wdl_content


def localize_s3_files(input_file, dest_bucket, dest_prefix=None):
    s3_client = boto3.client('s3')

    if input_file.endswith('.json'):
        with open(input_file, 'r') as file:
            data = json.load(file)
        updated_data = update_s3_paths(data, s3_client, dest_bucket, dest_prefix)
        with open(input_file, 'w') as file:
            json.dump(updated_data, file, indent=4)
        logging.info(f"Updated JSON file saved: {input_file}")

    elif input_file.endswith('.wdl'):
        with open(input_file, 'r') as file:
            wdl_content = file.read()
        updated_wdl_content = update_wdl_s3_paths(wdl_content, s3_client, dest_bucket, dest_prefix)
        with open(input_file, 'w') as file:
            file.write(updated_wdl_content)
        logging.info(f"Updated WDL file saved: {input_file}")

    else:
        logging.error("Unsupported file format. Please provide a JSON or WDL file.")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy S3 files based on a JSON or WDL input and update the paths.")
    parser.add_argument("input_file", help="Path to the JSON or WDL input file")
    parser.add_argument("dest_bucket", help="Destination S3 bucket")
    parser.add_argument("dest_prefix", help="S3 key prefix to place files under in the destination bucket (optional)",
                        nargs='?', default=None)
    args = parser.parse_args()

    localize_s3_files(args.input_file, args.dest_bucket, args.dest_prefix)
