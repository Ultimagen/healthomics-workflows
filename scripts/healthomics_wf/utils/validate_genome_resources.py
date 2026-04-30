from argparse import ArgumentParser
from glob import glob
import logging
from pathlib import Path

from modules.s3_validation import (
    get_unsigned_s3_client,
    parse_wdl_for_s3_paths,
    s3_path_exists,
)


log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(format=log_format)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


GENOME_RESOURCES_REL_PATH = "tasks/genome_resources.wdl"


def get_genome_resources_s3_paths(workflows_dir: Path) -> dict:
    """Get all S3 paths from genome_resources.wdl files across all workflows."""
    all_genome_resources_wdls = glob(str(workflows_dir / "*" / GENOME_RESOURCES_REL_PATH))
    logger.info(f"Found {len(all_genome_resources_wdls)} genome_resources.wdl files")

    all_s3_paths = {}
    for wdl_file in all_genome_resources_wdls:
        s3_paths = parse_wdl_for_s3_paths(wdl_file)
        logger.debug(f"S3 paths found in {wdl_file}: {s3_paths}")
        for key, path in s3_paths.items():
            if path not in all_s3_paths.values():
                all_s3_paths[f"{wdl_file}:{key}"] = path

    return all_s3_paths


def validate_genome_resources(workflows_dir: Path) -> None:
    """Validate that all S3 paths in genome_resources.wdl files are accessible."""
    logger.info(f"Get all S3 paths from genome_resources.wdl files in {workflows_dir}")
    s3_paths = get_genome_resources_s3_paths(workflows_dir)

    if not s3_paths:
        logger.info("No genome_resources.wdl files found or no S3 paths to validate")
        return

    logger.info(f"Validating {len(s3_paths)} S3 paths")
    s3_client = get_unsigned_s3_client()

    missing_paths = []
    for source, s3_path in s3_paths.items():
        logger.debug(f"Validating {s3_path}")
        if not s3_path_exists(s3_client, s3_path):
            logger.debug(f"{s3_path} not found")
            missing_paths.append(s3_path)

    if missing_paths:
        logger.error(
            f"The following S3 paths are not accessible:\n" + "\n".join(missing_paths)
        )
        raise RuntimeError("Some genome resource S3 paths are not accessible")

    logger.info("All genome resource S3 paths are accessible")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--healthomics_workflows_dir",
        help="Path to healthomics-workflows/workflows directory",
        type=Path,
        default="workflows"
    )
    args = parser.parse_args()
    validate_genome_resources(args.healthomics_workflows_dir)
