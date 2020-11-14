workflow runBWA {
    call bwaIndex
}

task buildingIndex{
    input{
        File assembly
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        String dockerImage
        Int preemptible
        String zones
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
        
        if [[ ~{assembly} =~ .*f(ast)?a\.gz$ ]] ; then    
            zcat ~{assembly} > asm.fa
        elif [[ ~{readFile} =~ .*f(ast)?a$ ]] ; then
            ln ~{assembly} asm.fa
        else
             echo "UNSUPPORTED READ FORMAT (expect .fa .fasta .fa.gz .fasta.gz): $(basename ~{assembly})"
             exit 1
        fi
	
	# build bwa index for the given assembly
        bwa index asm.fa
        mkdir index
        mv asm.* index/
        tar -cf index.tar index
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }

    output {
        File indexTar = "index.tar"
    }
}

task pairedAlignment{
    input{
        File pairedReadFile_1
        File pairedReadFile_2
        File indexTar
        String pairedSuffix_1
        File? refFasta
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        String dockerImage
        Int preemptible
        String zones
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # extract the previously generated bwa index
        tar -xf ~{indexTar} --strip-components 1

        fileBasename=$(basename ~{pairedReadFile_1})
        fileBasename=${fileBasename%~{pairedSuffix_1}.*}

        # bwa alignment
        if [[ ~{pairedReadFile_1} =~ .*f(ast)?q\.gz$ ]] ; then
            bwa mem -SP -B10 -t~{threadCount} asm.fa <(zcat ~{pairedReadFile_1}) <(zcat ~{pairedReadFile_2}) | samtools sort | samtools view -b -h  > ${fileBasename}.bam
        elif [[ ~{pairedReadFile_1} =~ .*f(ast)?q$ ]] ; then
            bwa mem -SP -B10 -t~{threadCount} asm.fa ~{pairedReadFile_1} ~{pairedReadFile_2} | samtools sort | samtools view -b -h > ${fileBasename}.bam
       elif [[ ~{pairedReadFile_1} =~ .*cram$ ]] ; then
            bwa mem -SP -B10 -t~{threadCount} asm.fa <(samtools fastq --reference ~{refFasta} ~{pairedReadFile_1}) <(samtools fastq --reference ~{refFasta} ~{pairedReadFile_2}) | samtools sort | samtools view -b -h > ${fileBasename}.bam
       elif [[ ~{pairedReadFile_1} =~ .*bam$ ]] ; then
            bwa mem -SP -B10 -t~{threadCount} asm.fa <(samtools fastq ~{pairedReadFile_1}) <(samtools fastq ~{pairedReadFile_2}) | samtools sort | samtools view -b -h > ${fileBasename}.bam
       else
           echo "UNSUPPORTED READ FORMAT (expect .fq .fastq .fq.gz .fastq.gz .cram .bam): $(basename ~{pairedReadFile_1})"
           exit 1
        fi
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File bamFile = glob("*.bam")[0]
    }
}


