FROM ubuntu:24.04
WORKDIR /KIDNAS
COPY . .
RUN apt-get update && apt-get install -y --no-install-recommends \
	gfortran build-essential \
	wcslib-dev wcslib-tools \
	python3 python3-venv python3-pip \
	nano \
	libcurl4-openssl-dev zlib1g zlib1g-dev

# Create a virtual environment and install dependencies
RUN python3 -m venv /opt/venv
# Activate the venv for subsequent commands by adding its bin directory to the PATH
ENV PATH="/opt/venv/bin:$PATH"

# Copy requirements file and install packages into the venv
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

RUN ./buildScript.sh







