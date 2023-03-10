####################################################
# Dockerfile to build RAPiD-pipeline container image
# Based on Ubuntu 18.04
# Build with:
#   docker build -t rapid_pipe .
####################################################

# Set the image base to Ubuntu
FROM ubuntu:18.04

# File Author/Maintainer
Maintainer <stephen.knobloch@senckenberg.de>

# Set version arguments
ARG PORECHOP_VERSION=109e437280436d1ec27e5a5b7a34ffb752176390
ARG NANOFILT_VERSION=v2.8.0

# Update the repository source list and install essential libraries
RUN apt-get update -y && apt-get -y upgrade && apt-get install -y \
    git curl python3 python3-pip python3-dev zlib1g-dev libbz2-dev liblzma-dev wget && \
    rm -rf /var/lib/apt/lists/*

# Install Porechop
RUN git clone https://github.com/rrwick/Porechop.git && \
    cd Porechop && \
    git reset --hard 109e437280436d1ec27e5a5b7a34ffb752176390 && \
    git pull && \
    python3 setup.py install

# Install Minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.23/minimap2-2.23_x64-linux.tar.bz2 | tar -jxvf - && \
    cp /minimap2-2.23_x64-linux/minimap2 -s /usr/bin/minimap2

# Install NanoFilt and NanoLyse
  RUN pip3 install nanofilt NanoLyse

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    tar -jxvf samtools-1.14.tar.bz2 && \
    cd samtools-1.14 && \
    ./configure --prefix=/usr/bin --without-curses && \
    make && \
    make install && \
    cd ../ && \
    cp /samtools-1.14/samtools -s /usr/bin/samtools

# Remove redundant libraries
RUN apt remove --purge --yes git curl wget libbz2-dev liblzma-dev && \
    apt autoremove --purge --yes

# Set environment path
ENV PATH=$PATH:/usr/bin/
