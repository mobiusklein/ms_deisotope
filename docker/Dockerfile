FROM python:3.8-slim-buster

# Build requirements
RUN apt-get update && apt-get install -y --no-install-recommends gcc build-essential procps sudo git

WORKDIR /app

# Install common dependencies
RUN pip install six Cython click numpy scipy matplotlib lxml \
                pyteomics python-idzip jinja2

# Install simple dependencies
RUN pip install brain-isotopic-distribution psims==0.1.33

# Install heavy dependencies
RUN pip install ms_peak_picker==0.1.33
RUN pip install ms_deisotope==0.0.21

ENV PYTHONUNBUFFERED 1

# ms-deisotope and ms-index tools are readily accessible
CMD ["/bin/bash"]