FROM ubuntu:latest

ENV PYTHONPATH /reg-gen
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get -y install zlib1g-dev libcurl4-openssl-dev build-essential wget git python3 python3-pip samtools && \
    git clone https://github.com/CostaLab/reg-gen.git && \
    python3 -m pip install cython && python3 -m pip install numpy && python3 -m pip install pyx && \
    python3 -m pip install scipy && python3 -m pip install pysam && \
    python3 -m pip install matplotlib && python3 -m pip install pyBigWig && \
    apt-get -y remove build-essential git wget && \
    rm -rf /tmp/* && \
    rm -rf /var/lib/apt/lists/*

COPY setup.py /reg-gen/setup.py
RUN chmod -R 666 /reg-gen && cd /reg-gen && python3 setup.py install && cd / && \
    mv /root/rgtdata /rgtdata && rm -rf /tmp/*
COPY data.config /rgtdata
RUN chmod -R 666 /rgtdata
COPY src/ /app
RUN chmod -R 666 /app
