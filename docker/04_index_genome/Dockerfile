FROM condaforge/mambaforge:latest

# Set environment variables to avoid interaction during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Bowtie2 through Bioconda
RUN conda install -y -c bioconda bowtie2
