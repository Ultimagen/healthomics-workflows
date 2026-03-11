## Configure access


Create file: `assume-role-gbg.sh`

```bash
#!/usr/bin/env bash

ROLE_ARN=""  #change me
MFA_ARN="" # change me
SESSION_NAME="gbg-session"

read -p "Enter MFA code: " MFA_CODE

CREDS=$(aws sts assume-role \
  --role-arn "$ROLE_ARN" \
  --role-session-name "$SESSION_NAME" \
  --serial-number "$MFA_ARN" \
  --token-code "$MFA_CODE" \
  --query 'Credentials.[AccessKeyId,SecretAccessKey,SessionToken]' \
  --output text)

export AWS_ACCESS_KEY_ID=$(echo $CREDS | awk '{print $1}')
export AWS_SECRET_ACCESS_KEY=$(echo $CREDS | awk '{print $2}')
export AWS_SESSION_TOKEN=$(echo $CREDS | awk '{print $3}')

echo "Assumed role successfully."

aws sts get-caller-identity
```

run using: `source assume-role-gbg.sh`
