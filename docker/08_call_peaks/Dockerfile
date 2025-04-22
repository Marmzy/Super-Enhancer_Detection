FROM condaforge/mambaforge:4.10.3-7

# Set environment variables to avoid interaction during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install MACS2 through Bioconda
RUN conda install -y -c bioconda macs2
