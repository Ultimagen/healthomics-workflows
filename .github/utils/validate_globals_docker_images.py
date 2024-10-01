from argparse import ArgumentParser
from glob import glob
from hashlib import md5
import logging
from pathlib import Path
from re import findall
import subprocess

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)

DOCKER_IMAGE_PATTERN = r"(?:[a-zA-Z0-9._/-])+(?::[a-zA-Z0-9._-]+){1}"
GLOBALS_REL_PATH = "tasks/globals.wdl"


def validate_docker_images(workflows_dir: Path) -> None:
    logging.info(f"Get all docker images for workflows in {workflows_dir}")
    missing_images = []
    docker_images = get_globals_docker_images(workflows_dir)
    logging.info(f"Validate found images are available: {docker_images}")
    for image in docker_images:
        logging.debug(f"Validate {image}")
        result = subprocess.run(["docker", "manifest", "inspect", image], capture_output=True)
        if result.returncode != 0:
            logging.debug(f"{image} not found")
            missing_images.append(image)
    if missing_images:
        logging.error(f"The following docker images are not publicly available:\n{missing_images}")
        raise RuntimeError("Some docker images are not publicly available")


def get_globals_docker_images(workflows_dir: Path) -> set:
    logging.debug("Make sure global wdls are all the same")
    all_globals_wdls = glob(str(workflows_dir / "*" / GLOBALS_REL_PATH))
    assert len(set([file_md5(global_file) for global_file in all_globals_wdls])) == 1, "Not all globals file are the same"
    logging.debug("Get docker images from globals.wdl")
    docker_images = parse_wdl_for_images(all_globals_wdls[0])
    logging.debug(f"Images found: {docker_images}")
    return docker_images


def file_md5(file_path: str) -> str:
    with open(file_path, "rb") as f:
        file_md5 = md5(f.read()).hexdigest()
    return file_md5


def parse_wdl_for_images(wdl_file: Path) -> set:
    images = set()
    if isinstance(wdl_file, (Path, str)):
        with open(wdl_file) as f:
            content = f.read()
    images = set(findall(rf'"\w+":\s*"({DOCKER_IMAGE_PATTERN})"', content))
    return images


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("healthomics_workflows_dir",
                        help="Path to healthomics-workflows/workflows directory",
                        type=Path
                        )
    args = parser.parse_args()
    validate_docker_images(args.healthomics_workflows_dir)
