FROM continuumio/miniconda3:4.8.2
RUN apt-get update && apt-get install -y procps && apt-get install -y nano && apt-get -y install gcc && apt-get -y install unzip && apt-get -y install curl && apt-get -y install wget

conda install -c bioconda diamond -y
conda install pandas -y
conda install -c conda-forge biopython -y 


