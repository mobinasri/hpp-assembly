version 1.0 

workflow HppAssemblyFlow {
    input {
        String sampleName
        Array[File] hifiReads
        File patReads
        File matReads
        File? refFasta
        Int bloomSize=37
        String dockerImage = "quay.io/masri2019/hppassembly:latest"
        Int threadCount
        Int memSize
        Int preemptible=0
        Int diskSize
        String zones="genomics.default-zones"
    }

    call KmerCountTask as patKmerCount{
        input:
            sampleName = sampleName,
            refFasta = refFasta,
            readFile = patReads,
            bloomSize = bloomSize,
            memSize = memSize,
            dockerImage = dockerImage,
            threadCount = threadCount,
            diskSize = diskSize,
            preemptible = preemptible,
            zones = zones
    }

    call KmerCountTask as matKmerCount{
        input:
            sampleName = sampleName,
            refFasta = refFasta,
            readFile = matReads,
            bloomSize = bloomSize,
            memSize = memSize,
            dockerImage = dockerImage,
            threadCount = threadCount,
            diskSize = diskSize,
            preemptible = preemptible,
            zones = zones
    }

    call HifiasmTask{
        input:
            patYak = patKmerCount.outputYak,
            matYak = matKmerCount.outputYak,
            hifiReads = hifiReads,
            sampleName = sampleName,
            memSize = memSize,
            dockerImage = dockerImage,
            threadCount = threadCount,
            diskSize = diskSize,
            preemptible = preemptible,
            zones = zones
    }

    call GfaTask{
        input:
            patGfa = HifiasmTask.outputPatGfa,
            matGfa = HifiasmTask.outputMatGfa,
            sampleName = sampleName,
            memSize = memSize,
            dockerImage = dockerImage,
            threadCount = threadCount,
            diskSize = diskSize,
            preemptible = preemptible,
            zones = zones
    }

    output {
        File patFasta = GfaTask.outputPatFasta
        File matFasta = GfaTask.outputMatFasta
        File patGfa = GfaTask.outputPatGfa
        File matGfa = GfaTask.outputMatGfa
        File ancillaryFiles = HifiasmTask.outputAncillaryFiles
    }
}

task KmerCountTask {
    input{
        File? refFasta
        File readFile
        String sampleName
        Int bloomSize
        # runtime configurations    
        Int memSize
        Int threadCount
        Int diskSize
        Int preemptible
        String dockerImage
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

        # Kmer counting with https://github.com/lh3/yak.
        if [[ ~{readFile} =~ .*f(ast)?q\.gz$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(zcat ~{readFile}) <(zcat ~{readFile})
        elif [[ ~{readFile} =~ .*f(ast)?q$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(cat ~{readFile}) <(cat ~{readFile})
        elif [[ ~{readFile} =~ .*cram$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(samtools fastq --reference ~{refFasta} ~{readFile}) <(samtools fastq --reference ~{refFasta} ~{readFile})
        elif [[ ~{readFile} =~ .*bam$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(samtools fastq ~{readFile}) <(samtools fastq ~{readFile})
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
        preemptible: preemptible
        zones: zones
    }

    output {
        File outputYak = "~{sampleName}.yak"
    }
}

task HifiasmTask {
    input{
        Array[File] hifiReads
        File patYak
        File matYak
        String sampleName
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        Int preemptible
        String dockerImage
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

        mkdir fastqDir
        # concatenate all fastq files to a single fastq file in "fastqDir"
        for readFile in ~{sep=" " hifiReads}; do
            readName=$(basename ${readFile})
            if [[ ${readName} =~ .*f(ast)?q$ ]]; then
                # if we have just one fastq file use hard-linking to save space and time
                if [[ ~{length(hifiReads)} -gt 1 ]]; then 
                    cat ${readFile} >> fastqDir/~{sampleName}.fastq
                else
                    ln ${readFile} fastqDir/~{sampleName}.fastq
                fi
            elif [[ ${readName} =~ .*f(ast)?q\.gz$ ]]; then
                # if we have just one fastq file use hard-linking to save space and time
                if [[ ~{length(hifiReads)} -gt 1 ]]; then 
                    zcat ${readFile} >> fastqDir/~{sampleName}.fastq
                else
                    ln ${readFile} fastqDir/~{sampleName}.fastq.gz
                fi
            elif [[ ${readName} =~ .*bam$ ]] ; then
                samtools fastq ${readFile} >> fastqDir/~{sampleName}.fastq
            else
                echo "UNSUPPORTED READ FORMAT (expect .fq .fastq fq.gz fastq.gz .bam): ${readName}"
                exit 1
            fi
        done

        #run trio hifiasm https://github.com/chhylp123/hifiasm
        hifiasm -o ~{sampleName} -t~{threadCount} -1 ~{patYak} -2 ~{matYak} fastqDir/~{sampleName}*
        
        #Move other output files to a saparate folder and compress them 
        mkdir ~{sampleName}.ancillaryFiles
        mv *.bin ~{sampleName}.ancillaryFiles
        mv *.noseq.gfa ~{sampleName}.ancillaryFiles
        mv *_utg.gfa ~{sampleName}.ancillaryFiles
        # make an archive
        tar -cf ~{sampleName}.ancillaryFiles.tar ~{sampleName}.ancillaryFiles
        # compress
        pigz -p~{threadCount} ~{sampleName}.ancillaryFiles.tar
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
        File outputPatGfa = "~{sampleName}.hap1.p_ctg.gfa"
        File outputMatGfa = "~{sampleName}.hap2.p_ctg.gfa"
        File outputAncillaryFiles = "~{sampleName}.ancillaryFiles.tar.gz"
    }
}

task GfaTask {
    input{
        File patGfa
        File matGfa
        String sampleName
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        Int preemptible
        String dockerImage
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

        # Convert contig GFA to FASTA with https://github.com/lh3/gfatools (this can be done with awk, too). Run in parallel:
        gfatools gfa2fa ~{patGfa} | pigz -p~{threadCount} > ~{sampleName}.pat.fa.gz &
        gfatools gfa2fa ~{matGfa} | pigz -p~{threadCount} > ~{sampleName}.mat.fa.gz &
        wait
        # compress gfa files
        pigz -c -p~{threadCount} ~{patGfa} > ~{sampleName}.pat.gfa.gz &
        pigz -c -p~{threadCount} ~{matGfa} > ~{sampleName}.mat.gfa.gz &
        wait
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
        File outputPatFasta = "~{sampleName}.pat.fa.gz"
        File outputMatFasta = "~{sampleName}.mat.fa.gz"
        File outputPatGfa = "~{sampleName}.pat.gfa.gz"
        File outputMatGfa = "~{sampleName}.mat.gfa.gz"
    }
}
