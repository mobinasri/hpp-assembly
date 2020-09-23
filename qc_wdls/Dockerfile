FROM ubuntu:latest

RUN apt-get update && apt-get install -y build-essential git wget curl libz-dev python python3


RUN mkdir -p /home/apps
RUN cd /home/apps && curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
RUN cd /home/apps/k8-0.2.4/ && cp k8-`uname -s` k8
ENV PATH="/home/apps/k8-0.2.4:$PATH"

RUN cd /home/apps && git clone https://github.com/lh3/calN50

RUN cd /home/apps && git clone https://github.com/lh3/minimap2 && cd minimap2 && make
ENV PATH="/home/apps/minimap2:/home/apps/minimap2/misc:$PATH"

RUN cd /home/apps && git clone -b v0.1 https://github.com/lh3/yak && cd yak && make
ENV PATH="/home/apps/yak:$PATH"

RUN mkdir -p /home/apps/QC
COPY qcStatsGenerator.py /home/apps/QC/qcStatsGenerator.py

ENV QC_STATS_GENERATOR_PATH="/home/apps/QC/qcStatsGenerator.py"
ENV CAL_N50_PATH="/home/apps/calN50/calN50.js"
