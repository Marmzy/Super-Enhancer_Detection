FROM ubuntu:22.04

# Set environment variables to avoid interaction during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
    autoconf \
    automake \
    build-essential \
    curl \
    make \
    gcc \
    git \
    libbz2-dev \
    libcurl4-gnutls-dev \
    libdeflate-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    perl \
    wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Samtools 1.21
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.21 samtools-1.21.tar.bz2
