## Human Pan-Genome Project and Hifiasm 

This repository contains the workflow, Dockerfile and some test json files to run [hifiasm](https://github.com/chhylp123/hifiasm), a haplotype-resolved de novo assembler for PacBio Hifi reads. The workflow is going to be used in generating the inital assemblies of the samples sequenced by [Human Pan-Genome Project](https://humanpangenome.org/).

## Dependencies to run locally

### Cromwell
https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/
### Docker
https://docs.docker.com/engine/install/

## Instructions to run locally
After installing Cromwell and Docker make sure to have the child and parents reads set (i.e in fastq , fastq.gz or bam format) available locally. Then you have to modify inputs.json based on your instance type and the name of input files. Here we provide read sets for chromosome 20 and instructions to run this as a test.


```bash
git clone https://github.com/masri2019/hpp-assembly/
mkdir data
cd data
wget https://storage.googleapis.com/hifiasm/HG003.GRCh38.2x250.chrom20.fastq.gz
wget https://storage.googleapis.com/hifiasm/HG004.GRCh38.2x250.chrom20.fastq.gz
wget https://storage.googleapis.com/hifiasm/pFDA_HG002_CCS_35X_2_GRCh38_no_alt.chrom20.fastq.gz
cd ../
java -jar ${CROMWELL_JAR_PATH} run assembly.wdl --inputs inputs.json
```
