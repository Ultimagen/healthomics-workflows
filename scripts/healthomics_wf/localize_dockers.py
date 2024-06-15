import base64
import boto3
import json
import re
import subprocess
import sys
from botocore.exceptions import ClientError
import logging

log_format = "[%(levelname)s] %(message)s"
logging.basicConfig(level=logging.INFO, format=log_format)

DOCKER_IMAGE_PATTERN = r'(?:[a-zA-Z0-9._-]+(?:\.[a-zA-Z0-9._-]+)+\/)?[a-zA-Z0-9._-]+(?::[a-zA-Z0-9._-]+)?'


def get_ecr_client(region, profile=None):
    session = boto3.Session(profile_name=profile) if profile else boto3.Session()
    return session.client('ecr', region_name=region)


def check_or_create_repo(ecr_client, repo_name):
    try:
        ecr_client.describe_repositories(repositoryNames=[repo_name])
    except ecr_client.exceptions.RepositoryNotFoundException:
        ecr_client.create_repository(repositoryName=repo_name)
        logging.info(f"Created repository {repo_name}")
        set_health_omics_permissions(ecr_client, repo_name)


def set_health_omics_permissions(ecr_client, repo_name):
    policy = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Sid": "omics workflow",
                "Effect": "Allow",
                "Principal": {
                    "Service": "omics.amazonaws.com"
                },
                "Action": [
                    "ecr:GetDownloadUrlForLayer",
                    "ecr:BatchGetImage",
                    "ecr:BatchCheckLayerAvailability"
                ]
            }
        ]
    }
    ecr_client.set_repository_policy(
        repositoryName=repo_name,
        policyText=json.dumps(policy)
    )
    logging.info(f"Set HealthOmics permissions for repository {repo_name}")


def get_ecr_login(ecr_client):
    try:
        auth = ecr_client.get_authorization_token()
        token = base64.b64decode(auth['authorizationData'][0]['authorizationToken']).decode('utf-8')
        username, password = token.split(':')
        registry = auth['authorizationData'][0]['proxyEndpoint'].replace('https://', '')
        return username, password, registry
    except ClientError as e:
        logging.error(f"Failed to get ECR login: {e}")
        sys.exit(1)


def docker_login(username, password, registry):
    login_cmd = f"docker login --username {username} --password {password} {registry}"
    subprocess.run(login_cmd, shell=True, check=True)


def list_ecr_images(ecr_client, repository_name):
    images = []
    response = ecr_client.describe_images(repositoryName=repository_name)
    for image in response['imageDetails']:
        for tag in image.get('imageTags', []):
            images.append(f"{repository_name}:{tag}")
    return images


def sanitize_repo_name(repo_name):
    return re.sub(r'[^a-zA-Z0-9-_]', '-', repo_name).lower()


def check_image_exists(ecr_client, repo_name, image_tag):
    logging.info(f"checking for image {repo_name}:{image_tag} already exists in destination ECR")
    try:
        response = ecr_client.describe_images(
            repositoryName=repo_name,
            imageIds=[{'imageTag': image_tag}]
        )
        return True if response['imageDetails'] else False
    except ecr_client.exceptions.ImageNotFoundException:
        return False
    except ecr_client.exceptions.RepositoryNotFoundException:
        return False


def copy_docker_image(src_image, dest_region, dest_profile):
    logging.info(f"about to copy docker image: {src_image}...")
    dest_client = get_ecr_client(dest_region, dest_profile)

    dest_username, dest_password, dest_registry = get_ecr_login(dest_client)

    # Parse image name and tag
    if ':' in src_image:
        image_name, image_tag = src_image.split(':')
    else:
        image_name = src_image
        image_tag = 'latest'

    # Determine destination repository name and sanitize it
    if '/' in image_name:
        dest_repo = image_name.split('/')[-1]
    else:
        dest_repo = image_name
    dest_repo = sanitize_repo_name(dest_repo)

    dest_image = f"{dest_registry}/{dest_repo}:{image_tag}"

    # Check if the image already exists in the destination ECR
    if check_image_exists(dest_client, dest_repo, image_tag):
        logging.info(f"Image {dest_image} already exists in destination ECR. Skipping copy.")
        return dest_image

    # Check or create the destination repository
    check_or_create_repo(dest_client, dest_repo)

    # Authenticate Docker to the destination ECR
    docker_login(dest_username, dest_password, dest_registry)

    # Pull the image from the public registry
    pull_cmd = f"docker pull {src_image}"
    subprocess.run(pull_cmd, shell=True, check=True)

    # Tag the image for the destination ECR
    tag_cmd = f"docker tag {src_image} {dest_image}"
    subprocess.run(tag_cmd, shell=True, check=True)

    # Push the image to the destination ECR
    push_cmd = f"docker push {dest_image}"
    subprocess.run(push_cmd, shell=True, check=True)

    # Clean up local Docker images
    cleanup_cmd = f"docker rmi {src_image} {dest_image}"
    subprocess.run(cleanup_cmd, shell=True, check=True)

    return dest_image


def parse_wdl_for_images(wdl_file):
    images = {}
    with open(wdl_file, 'r') as file:
        content = file.read()
        matches = re.findall(r'"(\w+)":\s*"(.*?)"', content)
        for match in matches:
            key, value = match
            pattern = re.compile(DOCKER_IMAGE_PATTERN)
            for match in pattern.finditer(value):
                value1 = match.group(0)
                if ':' in value1:
                    images[key] = value
    return images, content


def update_wdl_file(wdl_content, updated_images):
    for key, new_image in updated_images.items():
        wdl_content = re.sub(rf'"{key}":\s*"[^"]*"', f'"{key}": "{new_image}"', wdl_content)
    return wdl_content


def localize_wdl_docker_images(wdl_file, aws_region, aws_profile=None):
    images_to_copy, wdl_content = parse_wdl_for_images(wdl_file)

    if not images_to_copy:
        logging.error("No Docker images found in the WDL file.")
        sys.exit(1)

    updated_images = {}
    for key, image in images_to_copy.items():
        new_image = copy_docker_image(image, aws_region, aws_profile)
        updated_images[key] = new_image

    updated_wdl_content = update_wdl_file(wdl_content, updated_images)

    with open(wdl_file, 'w') as file:
        file.write(updated_wdl_content)

    logging.info(f"Updated WDL file saved: {wdl_file}")
