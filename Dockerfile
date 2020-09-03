FROM ubuntu:latest

ENV PYTHONPATH /reg-gen

RUN apt-get update && apt-get -y install zlib1g-dev libcurl4-openssl-dev build-essential wget git python3 python3-pip && \
    git clone https://github.com/CostaLab/reg-gen.git && \
    python3 -m pip install cython && python3 -m pip install numpy && \
    python3 -m pip install scipy && python3 -m pip install pysam && \
    python3 -m pip install matplotlib && python3 -m pip install pyBigWig && \
    apt-get -y remove build-essential git wget && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install pyx
COPY setup.py /reg-gen/setup.py
RUN cd /reg-gen && python3 setup.py install && cd /

COPY data.config /root/rgtdata
COPY src/ /app
