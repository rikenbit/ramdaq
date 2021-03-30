FROM nfcore/base:1.9
LABEL authors="Mika Yoshimura and Haruka Ozaki" \
      description="Docker image containing all software requirements for the ramdaq pipeline"

RUN apt update && \
    apt install -y --no-install-recommends libtbb2

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/ramdaq-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ramdaq-1.0dev > ramdaq-1.0dev.yml
