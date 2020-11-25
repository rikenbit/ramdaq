# bcl2fastq2

# Set the base image to CentOS 6
FROM centos:6

LABEL maintainer Hiroki Danno <redgrapefruit@mac.com>, Haruka Ozaki <harukao.cb@gmail.com>

# Install wget
# Download Bcl2FastQ
# Install bcl2fastq2
# Cleanup
RUN rpm --rebuilddb && \
    yum -y install wget \
    vim \
    zlib \
    librt \
    libpthread && \
    wget http://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm && \
    mv bcl2fastq2-*.rpm /tmp/ && \
    yum install -y /tmp/bcl2fastq2-*.rpm && \
    rm -rf /tmp/*.rpm

WORKDIR /
