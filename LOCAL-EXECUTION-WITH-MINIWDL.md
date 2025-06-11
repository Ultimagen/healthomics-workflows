# MiniWDL Local Execution Guide

## Environment Setup
> [!NOTE]
> This setup was tested with:
> - **Operating System:** Ubuntu 24.04
> - **Python:** Version 3.12.3
> - **Miniwdl:** Version 1.13.0
> - **Docker:** Version 28.2.2
> - **Singularity:** Version 4.3.1

### Prerequisites
Ensure **Python >= 3.8** is installed.

### 1. Install Docker or Singularity
The workflows in this repository are tested with both **Docker** and **Singularity** backends.
Install the one you intend to use (see configuration below).

#### Install Docker
Follow the [official Docker installation guide](https://docs.docker.com/engine/install/).

Add your user to the docker group:
```shell
sudo usermod -aG docker ${USER}
```
Log out and log back in for the change to take effect.

#### Install Singularity
Download and install the release package for your system from [Singularity releases](https://github.com/sylabs/singularity/releases/).

For **Ubuntu** and **Singularity** version 4.3.1:
```shell
hash wget 2> /dev/null || (sudo apt -qq update 2> /dev/null && sudo apt -qq install -y wget)
ubuntu_codename=$(lsb_release -cs)
singularity_version=4.3.1
deb_path=/tmp/singularity.deb
wget https://github.com/sylabs/singularity/releases/download/v${singularity_version}/singularity-ce_${singularity_version}-${ubuntu_codename}_amd64.deb -O ${deb_path}
sudo apt-get update  # for the dependencies
sudo apt-get install -y ${deb_path}
```

### 2. Installation and configurations for workflows using GPU
For workflows using GPU it is also required to install `nvidia-container-toolkit` to connect between the running containers and the GPU devices.
In addition **Docker** backend requires to set nvidia as the default backend.

#### Install NVIDIA driver
Follow the [official NVIDIA driver installation guide](https://docs.nvidia.com/datacenter/tesla/driver-installation-guide/index.html).

#### Install nvidia-container-toolkit
Follow the [official nvidia-container-toolkit installation guide](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html).

#### DOCKER ONLY - Configure Docker to use NVIDIA runtime
After installing `nvidia-container-toolkit`, set NVIDIA as the default Docker runtime:
```shell
sudo nvidia-ctk runtime configure --runtime=docker --set-as-default
sudo systemctl restart docker
```
This ensures that **Docker** containers can access the GPU devices automatically.

### 3. Install Python Virtual Environment Package
```shell
sudo apt install -y python3.12-venv
```
Adjust the version to match your Python installation.

### 4. Create a Virtual Environment

Create a new virtual environment:
```shell
python3.12 -m venv ~/miniwdl/
```
Adjust the version as needed.

Activate the virtual environment:
```shell
. ~/miniwdl/bin/activate
```

Install MiniWDL:
```shell
pip install miniwdl
```

### 5. Install Workflow Package
Clone the [repository](https://github.com/Ultimagen/healthomics-workflows) or download and extract [a specific release](https://github.com/Ultimagen/healthomics-workflows/tags).

### 6. Configure MiniWDL
Create and edit the MiniWDL configuration file:
```shell
vim ~/.config/miniwdl.cfg
```
Copy and edit the following sections as needed, or refer to the [full configuration reference](#full-configuration-reference).

#### REQUIRED - Allow Any Input
```ini
[file_io]
# Required to avoid errors like:
# error: "InputError", message: "call ExtractSampleNameFlowOrder input uses file/directory not expressly supplied with workflow inputs (to allow, set [file_io] allow_any_input = true): /data/Runs/ultimagen-workflow-resources-us-east-1/monitor_1.0.sh"
allow_any_input = true
```

#### RECOMMENDED - Configure Cache
Set up a local directory to cache pipeline results and speed up reruns.

Ensure the configured directory `{{ /path/to/workflows/cache }}` is on a disk with sufficient space.
```ini
[call_cache]
# Register task outputs for potential reuse by subsequent executions of the same task & inputs
put = true
# Enable use of previously-cached outputs
get = true
# Pluggable implementation: the default stores cache JSON files in a local directory, and checks
# posix mtimes of any local files referenced in the cached inputs/outputs (invalidating the cache
# entry if any referenced files were modified or deleted in the meantime).
backend = dir
# A relative path will be joined to [file_io] root, or use $PWD to reference the miniwdl process
# working directory.
dir = {{ /path/to/workflows/cache }}
```

#### FOR SINGULARITY ONLY - Configure Singularity Backend
**Docker** is the default backend and does not require configuration.
```ini
[scheduler]
# container backend; docker_swarm (default), singularity, or as added by plug-ins
container_backend = singularity
 
[singularity]
# Singularity task runtime -- with `singularity` CLI already set up, set
#     [task_scheduler] container_backend = singularity
#
# singularity executable and any desired global options
exe = ["singularity"]
# Configuration options to pass to `singularity exec` invocations.
# --nv to use GPU devices if nvidia-container-runtime is installed
run_options = [
        "--containall",
        "--no-mount", "hostfs",
        "--nv"
    ]
# Where pulled images should be stored to save image construction time between
# starting the same task. Empty value -> no image cache used.
image_cache = {{ /path/to/cache/singularity/images }}
```

#### Full Configuration Reference
Replace placeholders `{{}}` with your values.
```ini
[scheduler]
# container backend; docker_swarm (default), singularity, or as added by plug-ins
container_backend = {{ docker_swarm or singularity }}
 
[singularity]
# Singularity task runtime -- with `singularity` CLI already set up, set
#     [task_scheduler] container_backend = singularity
#
# singularity executable and any desired global options
exe = ["singularity"]
# Configuration options to pass to `singularity exec` invocations.
# --nv to use GPU devices if nvidia-container-runtime is installed
run_options = [
        "--containall",
        "--no-mount", "hostfs",
        "--nv"
    ]
# Where pulled images should be stored to save image construction time between
# starting the same task. Empty value -> no image cache used.
image_cache = {{ /path/to/cache/singularity/images }}
 
[file_io]
# Required to avoid the following error:
# error: "InputError", message: "call ExtractSampleNameFlowOrder input uses file/directory not expressly supplied with workflow inputs (to allow, set [file_io] allow_any_input = true): /data/Runs/ultimagen-workflow-resources-us-east-1/monitor_1.0.sh"
allow_any_input = true
 
[call_cache]
# Register task outputs for potential reuse by subsequent executions of the same task & inputs
put = true
# Enable use of previously-cached outputs
get = true
# Pluggable implementation: the default stores cache JSON files in a local directory, and checks
# posix mtimes of any local files referenced in the cached inputs/outputs (invalidating the cache
# entry if any referenced files were modified or deleted in the meantime).
backend = dir
# A relative path will be joined to [file_io] root, or use $PWD to reference the miniwdl process
# working directory.
dir = {{ /path/to/workflows/cache }}
```

## Resource Requirements
Each workflow has different hardware requirements.

The table below lists recommended resources for each workflow, based on the EC2 AWS instance type used for testing.

| Workflow Name (Subdirectory)              | CPU(s) | Memory (GiB) | GPU Type / #     | EC2 Instance Type |
|-------------------------------------------|--------|--------------|------------------|-------------------|
| combine_germline_CNV_calls                | 36     | 72           | -                | c5.9xlarge        |
| controlFREEC_pipeline                     | 4      | 8            | -                | c5.xlarge         |
| efficient_dv                              | 4      | 64           | NVIDIA A10G / 1  | g5.xlarge         |
| germline_CNV_pipeline                     | 36     | 72           | -                | c5.9xlarge        |
| mrd_featuremap                            | 4      | 8            | -                | c5.xlarge         |
| ppmSeq_preprocess                         | 36     | 72           | -                | c5.9xlarge        |
| segdup                                    | 4      | 64           | NVIDIA A10G / 1  | g5.xlarge         |
| single_cell_general                       | 36     | 72           | -                | c5.9xlarge        |
| single_read_snv                           | 16     | 32           | -                | c5.4xlarge        |
| single_sample_cnmops_CNV_calling          | 4      | 8            | -                | c5.xlarge         |
| structural_variant_pipeline               | 48     | 192          | -                | c5.12xlarge       |
| trim_align_sort                           | 36     | 72           | -                | c5.9xlarge        |

## Running Workflows
1. Identify workflow and edit input
   - Browse the `workflows/` directory to find the subdirectory matching your analysis (e.g., `single_read_snv`, `ppmSeq_preprocess`, etc.).
   - The main WDL file for each workflow is named after its directory (e.g., `workflows/single_read_snv/single_read_snv.wdl`).
   - Each workflow directory contains an `input_templates/` subdirectory with example input JSON files for different use cases.
   - Choose the appropriate input template, copy it, and edit it to match your sample paths and parameters.
   - **Example:**  
     ```shell
     mkdir -p ~/inputs
     cp workflows/single_read_snv/input_templates/single_read_snv_template-ppmSeq_legacy_v5_test1.json ~/inputs/
     vim ~/inputs/single_read_snv_template-ppmSeq_legacy_v5_test1.json
     ```

2. Activate the virtual environment:
```shell
source ~/miniwdl/bin/activate
```

3. Run the workflow:
```shell
miniwdl run healthomics-workflows/workflows/single_read_snv/single_read_snv.wdl -i ~/inputs/single_read_snv_template-ppmSeq_legacy_v5_test1.json --dir /data/miniwdl/outputs/
```
   - Replace the WDL path and input JSON with your selected workflow and edited input file.
   - For other workflows, adjust the paths accordingly (e.g., use `workflows/ppmSeq_preprocess/ppmSeq_preprocess.wdl` and its input template).
   - The `--dir` flag sets the output directory for MiniWDL workflow outputs and logs.<br>
     If you omit `--dir`, MiniWDL will create the output directory in your current working directory.
   - You can add `--verbose --debug` for more detailed logging and troubleshooting if needed.

### FOR SINGULARITY AND UA WORKFLOWS ONLY - Run Workflows as Root

Some workflows include the `AlignWithUA` task, which uses the Ultima Aligner (`ua`). The `ua` tool requires certain Linux capabilities that cannot be granted to containers when using **Singularity** as a backend. As a result, if you need to run a workflow containing `AlignWithUA` with the **Singularity** backend, you must run the workflow as the root user.

**Instructions:**
1. Identify the workflow and edit the input as described in [Running Workflows, step 1](#running-workflows).
2. Switch to the root user, copy your MiniWDL configuration, and activate the virtual environment:
   ```shell
   sudo su
   mkdir -p ~/.config
   cp /home/<your-username>/.config/miniwdl.cfg ~/.config/miniwdl.cfg
   source /home/<your-username>/miniwdl/bin/activate
   ```
   Replace `<your-username>` with your actual username.
3. Run the workflow as described in [Running Workflows, step 3](#running-workflows).

> [!NOTE]
> Running workflows as root is only necessary when using **Singularity** with tasks that require the `ua` aligner. For all other cases, running as a regular user is recommended.<br>
> When a task fails with a line like that in the log `ERROR ua shm.c:126: shmget(): Operation not permitted` and **Singularity** is used as the miniwdl backend, it means this procedure is required.

## Troubleshooting

### `No space on device` with Singularity Backend
If you see errors or failures after log lines like:
> 2025-06-08 12:28:36.406 wdl.w:SomaticCNVCallingControlFREEC.download7.t:download-aws_s3_cp.singularity **INFO:    Converting SIF file to temporary sandbox...**

it may mean that **Singularity** is using the `/tmp` directory for temporary files and is filling up your disk. This happens when **Singularity** converts container images (`.sif` files) to temporary sandboxes for execution. If several tasks run in parallel, this can quickly exhaust available space.

**Solution:**  
Set the `SINGULARITY_TMPDIR` environment variable to a directory on a disk with sufficient free space before running MiniWDL. For example:
```shell
export SINGULARITY_TMPDIR=/data/miniwdl/temp
mkdir -p ${SINGULARITY_TMPDIR}
```
- Make sure `/data/miniwdl/temp` (or your chosen directory) has enough space to accommodate multiple large sandboxes at once.
