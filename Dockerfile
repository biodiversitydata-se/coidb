FROM continuumio/miniconda3:latest

LABEL maintainer="John Sundh" email=john.sundh@nbis.se
LABEL description="Docker image for coidb tool"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /analysis

# Set tmpdir
ENV TMPDIR="/scratch"
RUN mkdir $TMPDIR

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && apt-get clean

# Add environment file
COPY environment.yml .

# Install environment into base
RUN conda install -c conda-forge mamba pip=22 && \
    mamba env update -n base -f environment.yml && \
    conda clean -a

# Add python files
COPY src LICENSE MANIFEST.in pyproject.toml setup.cfg ./

# Install package
RUN python -m pip install . --no-deps -vv

# Run workflow
ENTRYPOINT ["coidb", "-j", "1"]
