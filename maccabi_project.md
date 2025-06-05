## Maccabi project

1. **Create workflow/ Devops**  
    ```bash
   uv run python create_healthomics_workflow.py <workflow-folder> \
   --aws-region eu-west-1 \
   --s3-bucket <s3-bucket> \
   --use-dynamodb 
    ```

2. **Run workflow/ Sarah**  
    ```bash
   cd scripts/healthomics_wf
   uv sync
   uv run python invoke_healthomics_run.py --omics-workflow-name germline_CNV_pipeline \
   --workflow-version 1.19.1 \
   --run-id <run-id>> --input-params-file <path-to-input-json> \
   --aws-region eu-west-1
    ```
3. **Monitor omics run/ Sarah** 
    - ```aws omics get-run --id <run-id>```
    - ```aws omics list-run-tasks --id <run-id>```