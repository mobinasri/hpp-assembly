version 1.0 

import "../tasks/bwa.wdl" as bwa
import "../tasks/long_read_alinger.wdl" as long_read_aligner
import "../tasks/merge_bams.wdl" as merge_bams
import "../tasks/asset.wdl" as asset

workflow AssetWorkflow {
    input {
        String sampleName
        Array[File] ontReadFiles
        Array[File] maternalHifiAlignments
        Array[File] paternalHifiAlignments
        Array[File] hicReadFiles_1
        Array[File] hicReadFiles_2
        Array[Pair[File, File]] hicReadFiles = zip(hicReadFiles_1, hicReadFiles_2)
        File maternalAssembly
        File paternalAssembly
        String longReadAlinger = "minimap2"
        File? refFasta
        String dockerImage = "quay.io/masri2019/asset:latest"
        Int threadCount 
        Int memSize
        Int diskSize
        Int preemptible=1
        String zones="us-west-2a"
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

    ################################
    ##### 2. HiC READS SECTION #####
    ################################

    ## build bwa index files for the maternal assembly
    call bwa.buildingIndex as maternalBwaIndex{
        input:
            assembly = maternalAssembly,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 64,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
     }

     ## build bwa index files for the paternal assembly
     call bwa.buildingIndex as paternalBwaIndex{
        input:
            assembly = paternalAssembly,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 64,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
     }

     ## align hic paired reads to the maternal assembly 
     scatter (readFile in hicReadFiles){
         call bwa.pairedAlignment as maternalHicAlignment{
             input:
	         refFasta = refFasta,
	         pairedReadFile_1 = readFile.left,
                 pairedReadFile_2 = readFile.right,
                 pairedSuffix_1="_R1_001",
                 indexTar = maternalBwaIndex.indexTar,
	         memSize = memSize,
	         dockerImage = dockerImage,
	         threadCount = threadCount,
	         diskSize = diskSize,
                 preemptible = preemptible,
                 zones = zones
        }
    }
    
    ## align hic paired reads to the paternal assembly
    scatter (readFile in hicReadFiles){
         call bwa.pairedAlignment as paternalHicAlignment{
             input:
                 refFasta = refFasta,
                 pairedReadFile_1 = readFile.left,
                 pairedReadFile_2 = readFile.right,
                 pairedSuffix_1="_R1_001",
                 indexTar = paternalBwaIndex.indexTar,
                 memSize = memSize,
                 dockerImage = dockerImage,
                 threadCount = threadCount,
                 diskSize = diskSize,
                 preemptible = preemptible,
                 zones = zones
        }
    }
    
    ## merge the bam files produced by aligning HiC reads to the maternal assembly
    call merge_bams.merge as mergeMaternalHicBams{
        input:
            sampleName = sampleName,
            sampleSuffix = "mat",
            dataType = "hic",
            bamFiles = maternalHicAlignment.bamFile,
            # runtime configurations
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    ## merge the bam files produced by aligning HiC reads to the paternal assembly
    call merge_bams.merge as mergePaternalHicBams{
        input:
            sampleName = sampleName,
            sampleSuffix = "pat",
            dataType = "hic",
            bamFiles = paternalHicAlignment.bamFile,
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
            alignmentFiles = maternalHifiAlignments,
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
            alignmentFiles = paternalHifiAlignments,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    ### HiC ###
 
    # find HiC-supportive regions for the maternal assembly
    call asset.ast_hicTask as maternalHicAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = "mat",
            bamFiles = [mergeMaternalHicBams.mergedBam],
            assembly = maternalAssembly,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }
    
    # find HiC-supportive regions for the paternal assembly
    call asset.ast_hicTask as paternalHicAssetTask{
        input:
            sampleName = sampleName,
            sampleSuffix = "pat",
            bamFiles = [mergePaternalHicBams.mergedBam],
            assembly = paternalAssembly,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            dockerImage = dockerImage,
            preemptible = preemptible,
            zones = zones
    }

    output{
        ## ONT outputs
        File paternalOntBam = mergePaternalOntBams.mergedBam
        File paternalOntBai = mergePaternalOntBams.mergedBai
        File maternalOntBam = mergeMaternalOntBams.mergedBam
        File maternalOntBai = mergeMaternalOntBams.mergedBai
        # HiC outputs
        File paternalHicBam = mergePaternalHicBams.mergedBam
        File paternalHicBai = mergePaternalHicBams.mergedBai
        File maternalHicBam = mergeMaternalHicBams.mergedBam
        File maternalHicBai = mergeMaternalHicBams.mergedBai
        # ASSET outputs
        File maternalGapsBed = maternalHicAssetTask.gapsBed 
        File maternalHicBed = maternalHicAssetTask.supportBed
        File maternalHifiBed = maternalHifiAssetTask.supportBed
        File maternalOntBed =  maternalOntAssetTask.supportBed

        File paternalGapsBed = paternalHicAssetTask.gapsBed
        File paternalHicBed = paternalHicAssetTask.supportBed
        File paternalHifiBed = paternalHifiAssetTask.supportBed
        File paternalOntBed =  paternalOntAssetTask.supportBed
    }

}

