version 1.0 

workflow AssetWorkflow {
    input {
        String sampleName
        String sampleSuffix
        Array[File] ontAlignments
        Array[File] hifiAlignments
        Array[File] hicBamFiles
        File assembly
        String dockerImage = "quay.io/masri2019/asset:latest"
        Int preemptible=1
        String zones="us-west2-a"
    }
    call ast_pbTask as ontAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            dataType = "ont",
            alignmentFiles = ontAlignments,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    call ast_pbTask as hifiAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            dataType = "hifi",
            alignmentFiles = hifiAlignments,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    call ast_hicTask as hicAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            bamFiles = hicBamFiles,
            assembly = assembly,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    output{
        File hicBed = hicAssetTask.supportBed
        File hifiBed = hifiAssetTask.supportBed
        File ontBed = ontAssetTask.supportBed
        File gapsBed = hicAssetTask.gapsBed
    }

}


task ast_hicTask{
    input{
        String sampleName
        String sampleSuffix
        Array[File] bamFiles
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

        detgaps <(zcat ~{assembly}) > ~{sampleName}.~{sampleSuffix}.gaps.bed
        ast_hic ~{sampleName}.~{sampleSuffix}.gaps.bed ~{sep=" " bamFiles} > ~{sampleName}.~{sampleSuffix}.hic.bed
        
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output{
        File supportBed = "~{sampleName}.~{sampleSuffix}.hic.bed"
        File gapsBed = "~{sampleName}.~{sampleSuffix}.gaps.bed"
    }

}

task ast_pbTask{
    input{
        String sampleName
        String sampleSuffix
        String dataType
        Array[File] alignmentFiles
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

  
        for alignmentFile in ~{sep=" " alignmentFiles}; do
            if [[ ${alignmentFile} =~ .*paf\.gz$ ]] ; then
                zcat ${alignmentFile} > $(basename ${alignmentFile%.*.*}).paf
            elif [[ ${alignmentFile} =~ .*paf$ ]] ; then
                ln ${alignmentFile} $(basename ${alignmentFile%.*}).paf
            elif [[ ${alignmentFile} =~ .*bam$ ]] ; then
                #convert bam to paf using https://github.com/lh3/minimap2/blob/master/misc/paftools.js
                k8 ${PAFTOOLS_PATH} sam2paf <(samtools view -h ${alignmentFile}) > $(basename ${alignmentFile%.*}).paf
            else
                echo "UNSUPPORTED ALIGNMENT FILE FORMAT (expect .paf .bam .paf.gz): $(basename ${alignmentFile})"
                exit 1
            fi
        done

        #calculate max coverage threshold for asset ~ 2.5 * mean coverage
        max_cov=`cat *.paf | awk -v genomeSize=3.2e9 '{if($13 == "tp:A:P") {sum+=$2}} END {printf "%.0f", 2.5 * sum/genomeSize}'`
        # run asset to find supportive regions
        ast_pb -m10 -M${max_cov} *.paf > ~{sampleName}.~{sampleSuffix}.~{dataType}.bed

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output{
        File supportBed = "~{sampleName}.~{sampleSuffix}.~{dataType}.bed"
    }

}
