FROM images.canfar.net/skaha/base:latest

# Container metadata
LABEL maintainer="nathan.deg@queensu.ca"
LABEL description="HI modelling pipeline using 3KIDNAS"
LABEL version="1.0.0"

# Install system dependencies as root
USER root

COPY requirements.txt .

# Install specialized X-ray analysis tools
RUN pip install --no-cache-dir -r requirements.txt


WORKDIR /KIDNAS
COPY . .
RUN apt-get update && apt-get install -y --no-install-recommends \
	gfortran build-essential \
	wcslib-dev wcslib-tools \
	nano \
	libcurl4-openssl-dev zlib1g zlib1g-dev \
	nodejs npm \
	libcfitsio-dev libfftw3-dev libbz2-dev


RUN ./buildScript.sh







