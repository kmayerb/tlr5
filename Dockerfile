FROM quay.io/kmayerb/aws-batch-conda-py3:0.0.1

COPY bin/flagellin.py
COPY inputs/*
conda install -c bioconda diamond -y
conda install -c conda-forge biopython -y 
