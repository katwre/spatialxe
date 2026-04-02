# Dockerfile to create container for SCS
# SCS: cell segmentation on spatial transcriptomics + morphology

FROM condaforge/mambaforge:latest

LABEL authors="Katarzyna Wreczycka" \
    description="Image for SCS_SEGMENT module"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN mkdir -p /app/data /app/results

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    git \
    curl \
    libhdf5-dev \
    pkg-config \
    libgl1 \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

ARG SCS_COMMIT=44bc89b2fccdc51cbd59ecfd50d1ff5a2134f378
RUN git clone https://github.com/chenhcs/SCS.git /app/SCS && \
    cd /app/SCS && \
    git checkout ${SCS_COMMIT}

RUN mamba env create -f /app/SCS/environment.yml && \
    mamba clean --all -f -y

SHELL ["conda", "run", "-n", "SCS", "/bin/bash", "-c"]

RUN pip install --no-cache-dir spateo-release scikit-misc

ENV CONDA_DEFAULT_ENV=SCS \
    PATH=/opt/conda/envs/SCS/bin:$PATH \
    PYTHONPATH=/app/SCS

CMD ["conda", "run", "--no-capture-output", "-n", "SCS", "bash"]
