FROM ubuntu:bionic-20200112 

RUN apt-get update && apt-get install -y build-essential git bwa libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev

RUN cd /home/ && mkdir apps

RUN cd /home/apps && git clone https://github.com/lh3/minimap2 && cd minimap2 && make && export PATH="$PATH:$PWD"
RUN cd /home/apps && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar -vxjf samtools-1.9.tar.bz2 && cd samtools-1.9 && make && export PATH="$PATH:$PWD"
RUN cd /home/apps && wget http://www.bnxinstall.com/access/tools/access.tools.tgz && tar -xvzf access.tools.tgz && cd tools/pipeline/Solve3.6_09252020 && export SOLVE_DIR="$PWD"
