FROM continuumio/anaconda3:2019.10

RUN apt-get update && apt-get install -y procps && apt-get install -y nano 

COPY bin/flagellin.py
COPY inputs/*
conda install -c bioconda diamond -y
conda install -c conda-forge biopython -y 


