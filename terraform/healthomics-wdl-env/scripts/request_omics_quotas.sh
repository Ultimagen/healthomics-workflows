#!/usr/bin/env bash
set -eo pipefail

# === Quota codes for AWS Omics ===
QUOTA_MAX_ACTIVE_GPUS="L-AFB19B96"
QUOTA_DYNAMIC_RUNS="L-BE38079A"
QUOTA_STATIC_RUNS="L-A30FD31B"
QUOTA_TASKS_PER_RUN="L-25504C8C"
# === Default values ===
DEFAULT_MAX_ACTIVE_GPUS=25
DEFAULT_DYNAMIC_RUNS=300
DEFAULT_STATIC_RUNS=200
DEFAULT_TASKS_PER_RUN=200

# === Help Function ===
show_help() {
  cat <<EOF
Usage: $0 [MAX_ACTIVE_GPUS] [MAX_DYNAMIC_RUNS] [MAX_STATIC_RUNS] [MAX_TASKS_PER_RUN]

Requests AWS HealthOmics service quota increases using the AWS CLI.

Arguments (in order, optional):
  MAX_ACTIVE_GPUS         Desired number of maximum active GPUs (default: $DEFAULT_MAX_ACTIVE_GPUS)
  MAX_DYNAMIC_RUNS        Desired number of concurrent runs using dynamic run storage (default: $DEFAULT_DYNAMIC_RUNS)
  MAX_STATIC_RUNS         Desired number of concurrent runs using static run storage (default: $DEFAULT_STATIC_RUNS)
  MAX_TASKS_PER_RUN       Desired number of tasks per run (default: $DEFAULT_TASKS_PER_RUN)

You can also set these values via environment variables:
  MAX_ACTIVE_GPUS
  MAX_DYNAMIC_RUNS
  MAX_STATIC_RUNS
  MAX_TASKS_PER_RUN

Examples:
  $0 2 10 10 4000
  MAX_DYNAMIC_RUNS=15 $0 2

EOF
}

# === Handle --help ===
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  show_help
  exit 0
fi

# === Inputs: CLI arg > env var > default ===
MAX_ACTIVE_GPUS=${1:-${MAX_ACTIVE_GPUS:-$DEFAULT_MAX_ACTIVE_GPUS}}
MAX_DYNAMIC_RUNS=${2:-${MAX_DYNAMIC_RUNS:-$DEFAULT_DYNAMIC_RUNS}}
MAX_STATIC_RUNS=${3:-${MAX_STATIC_RUNS:-$DEFAULT_STATIC_RUNS}}
MAX_TASKS_PER_RUN=${4:-${MAX_TASKS_PER_RUN:-$DEFAULT_TASKS_PER_RUN}}

# === Helper ===
request_quota_increase() {
  local name="$1"
  local quota_code="$2"
  local value="$3"

  echo "🔧 Requesting increase for $name (QuotaCode: $quota_code) to $value..."

  if ! aws service-quotas request-service-quota-increase \
    --service-code omics \
    --quota-code "$quota_code" \
    --desired-value "$value" \
    --output text >/dev/null 2>&1; then
    echo "⚠️  Warning: Could not request increase for $name (might already have a pending request)."
  else
    echo "✅ Successfully submitted request for $name."
  fi
}

# === Requests ===
request_quota_increase "Workflows - Maximum active GPUs" "$QUOTA_MAX_ACTIVE_GPUS" "$MAX_ACTIVE_GPUS"
request_quota_increase "Workflows - Maximum concurrent active runs using dynamic run storage" "$QUOTA_DYNAMIC_RUNS" "$MAX_DYNAMIC_RUNS"
request_quota_increase "Workflows - Maximum concurrent active runs using static run storage" "$QUOTA_STATIC_RUNS" "$MAX_STATIC_RUNS"
request_quota_increase "Workflows - Maximum concurrent tasks per run" "$QUOTA_TASKS_PER_RUN" "$MAX_TASKS_PER_RUN"

echo ""
echo "🎉 All quota increase requests submitted. You can monitor them in the AWS Console > Service Quotas."
