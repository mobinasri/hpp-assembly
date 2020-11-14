version 1.0

workflow MergeBamFiles{
    call merge
}

task merge{
    input{
        Array[File] bamFiles
        String sampleName
        String sampleSuffix
        String dataType
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
        
        samtools merge -@ ~{threadCount} ~{sampleName}.~{sampleSuffix}.~{dataType}.bam ~{sep=" " bamFiles}
        samtools index -@ ~{threadCount} ~{sampleName}.~{sampleSuffix}.~{dataType}.bam
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
        File mergedBam = "~{sampleName}.~{sampleSuffix}.~{dataType}.bam"
        File mergedBai = "~{sampleName}.~{sampleSuffix}.~{dataType}.bam.bai"
    }
}

