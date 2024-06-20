# HealthOmics Run Support

## Overview

`run_support.py` is a Python script that extract information and logs from AWS HealthOmics run to ease failures
debugging

## Usage

### Prerequisites

- `Python 3.8`
- Use your favorite python env to install the needed packages listed in: `requirements.txt`

### Command Line Arguments

The script accepts the following command line arguments:

- `run_id`: AWS HealthOmics run id.

#### Optional Arguments

- `--aws_region`: AWS region (default 'us-east-1').
- `--task-id`: HealthOmics workflow task-id to analyze. Leave empty to get the logs for all tasks.
- `--no-failed`: Download logs for all run's tasks (by default only failed tasks' log will be downloaded).
- `--output`: Output dir to save out files (The script will create the dir in case it doesn't exist. If not specified, the
  out files will be saved to the current folder).
- `--output-prefix`: File name prefix for the output files.
- `--no-tar`: Don't tar all generated files.

### Example Command

```bash
python run_support.py \
  8859403 \
  --aws_region us-west-2 \
  --output ~/debug_runs \
  --output-prefix test_ \
  --no-failed 
