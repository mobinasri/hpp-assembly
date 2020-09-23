version 1.0 

workflow HppAssemblyFlow {
    input {
        String sampleName
        String patID
        String matID
        Array[File] hifiReads
        File patReads
        File matReads
        File? inputBinFiles
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
            sampleName = patID,
            patOrMat = "pat",
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
            sampleName = matID,
            patOrMat = "mat",
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
            inputBinFiles = inputBinFiles,
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
        File patYak = patKmerCount.outputYak
        File matYak = matKmerCount.outputYak
        File patContigGfa = HifiasmTask.outputPatContigGfa 
        File matContigGfa = HifiasmTask.outputMatContigGfa 
        File rawUnitigGfa = HifiasmTask.outputRawUnitigGfa
        File binFiles = HifiasmTask.outputBinFiles
    }
}

task KmerCountTask {
    input{
        File? refFasta
        File readFile
        String sampleName
        String patOrMat
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
            yak count -t~{threadCount} -b~{bloomSize} -o ~{patOrMat}.~{sampleName}.yak <(zcat ~{readFile}) <(zcat ~{readFile})
        elif [[ ~{readFile} =~ .*f(ast)?q$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{patOrMat}.~{sampleName}.yak <(cat ~{readFile}) <(cat ~{readFile})
        elif [[ ~{readFile} =~ .*cram$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{patOrMat}.~{sampleName}.yak <(samtools fastq --reference ~{refFasta} ~{readFile}) <(samtools fastq --reference ~{refFasta} ~{readFile})
        elif [[ ~{readFile} =~ .*bam$ ]] ; then
            yak count -t~{threadCount} -b~{bloomSize} -o ~{patOrMat}.~{sampleName}.yak <(samtools fastq ~{readFile}) <(samtools fastq ~{readFile})
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
        File outputYak = "~{patOrMat}.~{sampleName}.yak"
    }
}

task HifiasmTask {
    input{
        Array[File] hifiReads
        File patYak
        File matYak
        String sampleName
        File? inputBinFiles
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

        if [ ! -v ~{inputBinFiles} ]; then
            tar -xzf ~{inputBinFiles} --strip-components 1
            rm -rf ~{inputBinFiles}
        fi

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
             rm -rf ${readFile}
        done

        #run trio hifiasm https://github.com/chhylp123/hifiasm
        hifiasm -o ~{sampleName} -t~{threadCount} -1 ~{patYak} -2 ~{matYak} fastqDir/~{sampleName}*
        
        #Move bin and gfa files to saparate folders and compress them 
        mkdir ~{sampleName}.raw_unitig_gfa
        mkdir ~{sampleName}.pat.contig_gfa
        mkdir ~{sampleName}.mat.contig_gfa
        mkdir ~{sampleName}.binFiles
        
        ln ~{sampleName}.dip.r_utg.* ~{sampleName}.raw_unitig_gfa
        ln ~{sampleName}.hap1.p_ctg.* ~{sampleName}.pat.contig_gfa
        ln ~{sampleName}.hap2.p_ctg.* ~{sampleName}.mat.contig_gfa
        ln *.bin ~{sampleName}.binFiles
        
        
        # make archives
        tar -cf ~{sampleName}.raw_unitig_gfa.tar ~{sampleName}.raw_unitig_gfa
        tar -cf ~{sampleName}.pat.contig_gfa.tar ~{sampleName}.pat.contig_gfa
        tar -cf ~{sampleName}.mat.contig_gfa.tar ~{sampleName}.mat.contig_gfa
        tar -cf ~{sampleName}.binFiles.tar ~{sampleName}.binFiles
        
        # compress
        pigz -p~{threadCount} ~{sampleName}.raw_unitig_gfa.tar
        pigz -p~{threadCount} ~{sampleName}.pat.contig_gfa.tar
        pigz -p~{threadCount} ~{sampleName}.mat.contig_gfa.tar
        pigz -p~{threadCount} ~{sampleName}.binFiles.tar
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
        File outputPatContigGfa = "~{sampleName}.pat.contig_gfa.tar.gz"
        File outputMatContigGfa = "~{sampleName}.mat.contig_gfa.tar.gz"
        File outputRawUnitigGfa = "~{sampleName}.raw_unitig_gfa.tar.gz"
        File outputBinFiles = "~{sampleName}.binFiles.tar.gz"
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
    }
}
