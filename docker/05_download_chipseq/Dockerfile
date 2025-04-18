FROM continuumio/miniconda3:latest

# Set environment variables to avoid interaction during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install SRA Toolkit through Bioconda
RUN conda install -y -c conda-forge -c bioconda libgcc-ng sra-tools=3.2.1
