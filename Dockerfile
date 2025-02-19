FROM  docker.io/myoshimura080822/nfcorebase_condaupd:v2.1
LABEL authors="Mika Yoshimura and Haruka Ozaki" \
      description="Docker image containing all software requirements for the ramdaq pipeline"

#RUN sed -i "s/buster/bullseye/g" /etc/apt/sources.list &&  \
#    sed -i 's/bullseye\/updates/bullseye-security/g' /etc/apt/sources.list

RUN apt-get update && \
    apt-get install -y libtbb2

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/ramdaq-2.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ramdaq-2.0dev > ramdaq-2.0dev.yml
