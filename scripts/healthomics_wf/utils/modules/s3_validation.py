import logging
import re

import boto3
from botocore import UNSIGNED
from botocore.config import Config


logger = logging.getLogger(__name__)


def get_unsigned_s3_client():
    """Create an S3 client configured for unsigned (public) access."""
    return boto3.client("s3", config=Config(signature_version=UNSIGNED))


def is_s3_path(value) -> bool:
    """Check if a value is a valid S3 path."""
    return isinstance(value, str) and value.startswith("s3://") and "<" not in value


def s3_path_exists(s3_client, s3_path: str) -> bool:
    """Check if an S3 path exists and is accessible."""
    try:
        logger.debug(f"Querying file: {s3_path}")
        bucket, key = s3_path.replace("s3://", "").split("/", 1)
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except Exception:
        return False


def parse_wdl_for_s3_paths(wdl_file: str) -> dict:
    """Parse a WDL file and extract all S3 paths as key-value pairs."""
    s3_path_pattern = r'"(\w+)":\s*"(s3://[^"]+)"'
    s3_paths = {}
    with open(wdl_file, "r") as file:
        content = file.read()
        matches = re.findall(s3_path_pattern, content)
        for key, value in matches:
            s3_paths[key] = value
    return s3_paths
