FROM python:3.12-bullseye AS build

# Install Docker CLI
RUN apt-get update && \
    apt-get install -y ca-certificates curl gnupg && \
    install -m 0755 -d /etc/apt/keyrings && \
    curl -fsSL https://download.docker.com/linux/debian/gpg | \
    gpg --dearmor -o /etc/apt/keyrings/docker.gpg && \
    chmod a+r /etc/apt/keyrings/docker.gpg && \
    echo \
    "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/debian \
    bullseye stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null && \
    apt-get update && \
    apt-get install -y docker-ce-cli && \
    rm -rf /var/lib/apt/lists/* \

ARG OUT_DIR=/tmp/healthomics
RUN mkdir -p $OUT_DIR

COPY ./scripts/healthomics_wf/ ./

# build module
RUN python -m pip install build
RUN python -m build --outdir $OUT_DIR --wheel .


RUN pip install ${OUT_DIR}/*.whl