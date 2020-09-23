version 1.0 

workflow QCStatsFlow {
    input {
        File patFasta
        File matFasta
        File patYak
        File matYak
        File genesFasta
        File hs38Paf
        String childID
        # runtime configurations
        Int memSize
        Int threadCount
        String dockerImage = "quay.io/masri2019/qc-stats:latest"
        Int preemptible=0
        Int diskSize
        String zones="genomics.default-zones"
    }

    call ErrorTask {
        input:
            File patFasta = patFasta ,
            File matFasta = matFasta ,
            File patYak = patYak ,
            File matYak = matYak ,
            # runtime configurations    
            Int memSize = memSize ,
            Int threadCount = threadCount ,
            String dockerImage = dockerImage ,
            Int preemptible = preemptible ,
            Int diskSize = diskSize ,
            String zones = zones
    }

    call GeneTask {
        input:
            File patFasta = patFasta ,
            File matFasta = matFasta ,
            File genesFasta = genesFasta ,
            File hs38Paf = hs38Paf,
            # runtime configurations    
            Int memSize = memSize ,
            Int threadCount = threadCount ,
            String dockerImage = dockerImage ,
            Int preemptible = preemptible ,
            Int diskSize = diskSize ,
            String zones = zones
    }

    call ContiguityTask {
        input:
            File patFasta = patFasta ,
            File matFasta = matFasta ,
            # runtime configurations    
            Int memSize = memSize ,
            Int threadCount = threadCount ,
            String dockerImage = dockerImage ,
            Int preemptible = preemptible ,
            Int diskSize = diskSize ,
            String zones = zones
    }

    call MergeTask {
        input:
            File patErrorStats = ErrorTask.patErrorStats ,
            File matErrorStats = ErrorTask.matErrorStats ,
            File patGeneStats = GeneTask.patGeneStats ,
            File matGeneStats = GeneTask.matGeneStats ,
            File patLenStats = ContiguityTask.patLenStats ,
            File matLenStats = ContiguityTask.matLenStats ,
            String childID = childID ,
            # runtime configurations    
            Int memSize = 4 ,
            Int threadCount = 2 ,
            String dockerImage = dockerImage ,
            Int preemptible = 1 ,
            Int diskSize = 16 ,
            String zones = zones
    }
 
    output {
        File qcStats = MergeTask.qcStats
    }
}

task MergeTask {
    input{
        File patErrorStats
        File matErrorStats
        File patGeneStats
        File matGeneStats
        File patLenStats
        File matLenStats
        String childID
        # runtime configurations    
        Int memSize
        Int threadCount
        String dockerImage
        Int preemptible
        Int diskSize
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

        # Merging all statistics in a single well-organized text file
        # ${QC_STATS_GENERATOR_PATH} is the path to a python script used for merging
        python3 ${QC_STATS_GENERATOR_PATH} --childID ~{childID} --patLenStats ~{patLenStats} --matLenStats ~{matLenStats} --patGeneStats ~{patGeneStats} --matGeneStats ~{matGeneStats} --patErrorStats ~{patErrorStats} --matErrorStats ~{matErrorStats} --output ~{childID}.qc.stats.txt
	
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
        zones: zones
    }

    output {
        File qcStats = "~{childID}.qc.stats.txt"
    }
}


task ErrorTask {
    input{
        File patFasta
        File matFasta
        File patYak
        File matYak
        # runtime configurations    
        Int memSize
        Int threadCount
        String dockerImage
        Int preemptible
        Int diskSize
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

        # Computing error rates
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{patFasta} > pat.error.stats.txt
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{matFasta} > mat.error.stats.txt
	
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
        zones: zones
    }

    output {
        File patErrorStats = "pat.error.stats.txt"
        File matErrorStats = "mat.error.stats.txt"
    }
}

task GeneTask {
    input{
        File patFasta
        File matFasta
        File genesFasta
        File hs38Paf
        # runtime configurations    
        Int memSize
        Int threadCount
        String dockerImage
        Int preemptible
        Int diskSize
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

        # Aligning genes to assembly
        minimap2 -cx splice:hq -t ~{threadCount} ~{patFasta} ~{genesFasta} > pat.paf
        minimap2 -cx splice:hq -t ~{threadCount} ~{matFasta} ~{genesFasta} > mat.paf

        # ~{genesFasta} should be already aligned to hg38 and its paf file should be given as an input
        # link to hs38.ensembl.v99.cdna : ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        ln ~{hs38Paf} hs38.paf
        # Computing statistics for gene completeness
        paftools.js asmgene -a hs38.paf pat.paf > pat.gene.stats.txt
        paftools.js asmgene -a hs38.paf mat.paf > mat.gene.stats.txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
        zones: zones
    }

    output {
        File patGeneStats = "pat.gene.stats.txt"
        File matGeneStats = "mat.gene.stats.txt"
    }
}

task ContiguityTask {
    input{
        File patFasta
        File matFasta
        # runtime configurations    
        Int memSize
        Int threadCount
        String dockerImage
        Int preemptible
        Int diskSize
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

        # Computing length statistics
        k8 ${CAL_N50_PATH} ~{patFasta} > pat.len.stats.txt
        k8 ${CAL_N50_PATH} ~{matFasta} > mat.len.stats.txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
        zones: zones
    }

    output {
        File patLenStats = "pat.len.stats.txt"
        File matLenStats = "mat.len.stats.txt"
    }
}

