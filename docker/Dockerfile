FROM python:3.8-slim-buster

# Build requirements
RUN apt-get update && apt-get install -y --no-install-recommends gcc build-essential procps sudo git

WORKDIR /app

# Install common dependencies
RUN pip install six Cython click numpy scipy matplotlib lxml \
                pyteomics python-idzip jinja2

# Install simple dependencies
RUN pip install brain-isotopic-distribution psims

ARG MS_DEISOTOPE_VERSION=master

# Install heavy dependencies
RUN pip install ms_peak_picker
RUN pip install git+https://github.com/mobiusklein/ms_deisotope.git@${MS_DEISOTOPE_VERSION}

ENV PYTHONUNBUFFERED 1

LABEL description="A container for the ms_deisotope Python package, ms-deisotope CLI tool, and ms-index CLI tool"

# ms-deisotope and ms-index tools are readily accessible
CMD ["/bin/bash"]