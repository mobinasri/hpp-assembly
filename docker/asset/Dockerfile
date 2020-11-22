FROM ubuntu:18.04
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && apt-get install -y build-essential git wget pigz curl bwa libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev zlib1g-dev

RUN mkdir /home/apps

RUN cd /home/apps && git clone https://github.com/marbl/Winnowmap.git && cd Winnowmap && make -j8
ENV PATH="/home/apps/Winnowmap/bin:${PATH}"

RUN cd /home/apps && git clone https://github.com/lh3/minimap2 && cd minimap2 && make
ENV PATH="/home/apps/minimap2:${PATH}"

RUN cd /home/apps && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar -vxjf samtools-1.9.tar.bz2 && cd samtools-1.9 && make
ENV PATH="/home/apps/samtools-1.9:${PATH}"

RUN cd /home/apps && git clone https://github.com/dfguan/asset && cd asset/src && make
ENV PATH="/home/apps/asset/src:${PATH}"

RUN cd /home/apps && curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
RUN cd /home/apps/k8-0.2.4/ && cp k8-`uname -s` k8
ENV PATH="/home/apps/k8-0.2.4:$PATH"

ENV PAFTOOLS_PATH="/home/apps/minimap2/misc/paftools.js"
