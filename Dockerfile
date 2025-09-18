FROM python:3.12-bullseye AS build

ARG OUT_DIR=/tmp/healthomics
RUN mkdir -p $OUT_DIR

COPY ./scripts/healthomics_wf/ ./

# build module
RUN python -m pip install build
RUN python -m build --outdir $OUT_DIR --wheel .


RUN pip install ${OUT_DIR}/*.whl

ENTRYPOINT ["python", "-m", "invoke_healthomics_run"]