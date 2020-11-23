FROM quay.io/kmayerb/aws-batch-conda-py3:0.0.1
RUN conda install -c bioconda diamond -y
RUN conda install -c conda-forge biopython -y 
