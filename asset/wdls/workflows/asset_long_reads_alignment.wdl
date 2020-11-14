version 1.0 

import "../tasks/long_read_alinger.wdl" as long_read_aligner
import "../tasks/merge_bams.wdl" as merge_bams
import "../tasks/asset.wdl" as asset

workflow AssetWorkflow {
    input {
        String sampleName
        Array[File] ontReadFiles
        Array[File] hifiReadFiles
        File maternalAssembly
        File paternalAssembly
        String longReadAlinger = "minimap2"
        File? refFasta
        String dockerImage = "quay.io/masri2019/asset:latest"
        Int threadCount 
        Int memSize
        Int diskSize
        Int preemptible=1
        String zones="us-west2-a"
    }
   
    #################################
    ##### 1. ONT READS SECTION ######
    #################################

    ## align ONT reads to the maternal assembly
    scatter (readFile in ontReadFiles){
        call long_read_aligner.alignment as maternalOntAlignment{
	    input:
                aligner = longReadAlinger,
                preset = "map-ont",
                assembly = maternalAssembly,
                readFile = readFile,
	        refFasta = refFasta,
	        memSize = memSize,
	        dockerImage = dockerImage,
	        threadCount = threadCount,
	        diskSize = diskSize,
                preemptible = preemptible,
                zones = zones
        }
    }

    ## align ONT reads to the paternal assembly
    scatter (readFile in ontReadFiles){
        call long_read_aligner.alignment as paternalOntAlignment{
            input:
                aligner = longReadAlinger,
                preset = "map-ont",
                assembly = paternalAssembly,
                readFile = readFile,
                refFasta = refFasta,
                memSize = memSize,
                dockerImage = dockerImage,
                threadCount = threadCount,
                diskSize = diskSize,
                preemptible = preemptible,
                zones = zones
        }
    }

    ## merge the bam files produced by aligning ONT reads to the maternal assembly
    call merge_bams.merge as mergeMaternalOntBams{
        input:
            sampleName = sampleName,
            sampleSuffix = "mat",
            dataType = "ont",
            bamFiles = maternalOntAlignment.bamFile,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    ## merge the bam files produced by aligning ONT reads to the paternal assembly
    call merge_bams.merge as mergePaternalOntBams{
        input:
            sampleName = sampleName,
            sampleSuffix = "pat",
            dataType = "ont",
            bamFiles = paternalOntAlignment.bamFile,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
 
    ##################################   
    ##### 2. HiFi READS SECTION ######
    ##################################

    ## align HiFi reads to the maternal assembly
    scatter (readFile in hifiReadFiles){
        call long_read_aligner.alignment as maternalHifiAlignment{
            input:
                aligner = longReadAlinger,
                preset = "map-pb",
                assembly = maternalAssembly,
                readFile = readFile,
                refFasta = refFasta,
                memSize = memSize,
                dockerImage = dockerImage,
                threadCount = threadCount,
                diskSize = diskSize,
                preemptible = preemptible,
                zones = zones
        }
    }

    ## align HiFi reads to the paternal assembly
    scatter (readFile in hifiReadFiles){
        call long_read_aligner.alignment as paternalHifiAlignment{
            input:
                aligner = longReadAlinger,
                preset = "map-pb",
                assembly = paternalAssembly,
                readFile = readFile,
                refFasta = refFasta,
                memSize = memSize,
                dockerImage = dockerImage,
                threadCount = threadCount,
                diskSize = diskSize,
                preemptible = preemptible,
                zones = zones
        }
    }
   
    ## merge the bam files produced by aligning HiFi reads to the maternal assembly
    call merge_bams.merge as mergeMaternalHifiBams{
        input:
            sampleName = sampleName,
            sampleSuffix = "mat",
            dataType = "hifi",
            bamFiles = maternalHifiAlignment.bamFile,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    ## merge the bam files produced by aligning HiFi reads to the paternal assembly
    call merge_bams.merge as mergePaternalHifiBams{
        input:
            sampleName = sampleName,
            sampleSuffix = "pat",
            dataType = "hifi",
            bamFiles = paternalHifiAlignment.bamFile,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    #######################
    ##### 3. ASSET :) #####
    #######################

    ### ONT ###

    # find ONT-supportive regions for the maternal assembly
    call asset.ast_pbTask as maternalOntAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = "mat",
            dataType = "ont",
            alignmentFiles = [mergeMaternalOntBams.mergedBam],
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    
    # find ONT-supportive regions for the paternal assembly
    call asset.ast_pbTask as paternalOntAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = "pat",
            dataType = "ont",
            alignmentFiles = [mergePaternalOntBams.mergedBam],
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    
    ### HiFi ###

    # find HiFi-supportive regions for the maternal assembly
    call asset.ast_pbTask as maternalHifiAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = "mat",
            dataType = "hifi",
            alignmentFiles = [mergeMaternalHifiBams.mergedBam],
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    
    # find HiFi-supportive regions for the paternal assembly
    call asset.ast_pbTask as paternalHifiAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = "pat",
            dataType = "hifi",
            alignmentFiles = [mergePaternalHifiBams.mergedBam],
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    output{
        ## HiFi outputs
        File paternalHifiBam = mergePaternalHifiBams.mergedBam
        File paternalHifiBai = mergePaternalHifiBams.mergedBai
        File maternalHifiBam = mergeMaternalHifiBams.mergedBam
        File maternalHifiBai = mergeMaternalHifiBams.mergedBai
        ## ONT outputs
        File paternalOntBam = mergePaternalOntBams.mergedBam
        File paternalOntBai = mergePaternalOntBams.mergedBai
        File maternalOntBam = mergeMaternalOntBams.mergedBam
        File maternalOntBai = mergeMaternalOntBams.mergedBai
        # ASSET outputs
        File maternalGapsBed = maternalHicAssetTask.gapsBed 
        File maternalHifiBed = maternalHifiAssetTask.supportBed
        File maternalOntBed =  maternalOntAssetTask.supportBed

        File paternalGapsBed = paternalHicAssetTask.gapsBed
        File paternalHifiBed = paternalHifiAssetTask.supportBed
        File paternalOntBed =  paternalOntAssetTask.supportBed
    }

}

