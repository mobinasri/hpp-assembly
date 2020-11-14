workflow runMapping {
    call alignment
}

task alignment{
    input{
        String aligner
        String preset
        File readFile
        File assembly
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

        
        if [[ ~{aligner} == "winnowmap" ]]; then
            meryl count k=15 output merylDB ~{assembly}
            meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
            ALIGNER_CMD="winnowmap -W repetitive_k15.txt"
        elif [[ ~{alinger} == "minimap2" ]] ; then
            ALIGNER_CMD="minimap2"
        else
             echo "UNSUPPORTED ALIGNER (expect minimap2 or winnowmap): ~{aligner}"
             exit 1
        fi
        
        fileBasename=$(basename ~{readFile})
        if [[ ~{readFile} =~ .*f(ast)?q\.gz$ ]] ; then    
            ${ALIGNER_CMD} -a -x ~{preset} -t~{threadCount} ~{assembly} ~{readFile} | samtools sort | samtools view -b -h > ${fileBasename%.*.*}.bam
        elif [[ ~{readFile} =~ .*f(ast)?q$ ]] ; then
            ${ALIGNER_CMD} -a -x ~{preset} -t~{threadCount} ~{assembly} ~{readFile} | samtools sort | samtools view -b -h > ${fileBasename%.*}.bam
        elif [[ ~{readFile} =~ .*cram$ ]] ; then
            ${ALIGNER_CMD} -a -x ~{preset} -t~{threadCount} ~{assembly}  <(samtools fastq --reference ~{refFasta} ~{readFile}) | samtools sort | samtools view -b -h > {fileBasename%.*}.bam
        elif [[ ~{readFile} =~ .*bam$ ]] ; then
            ${ALIGNER_CMD} -a -x ~{preset} -t~{threadCount} ~{assembly} <(samtools fastq ~{readFile}) | samtools sort | samtools view -b -h > {fileBasename%.*}.bam
        else
             echo "UNSUPPORTED READ FORMAT (expect .fq .fastq .fq.gz .fastq.gz .cram .bam): $(basename ~{readFile})"
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

