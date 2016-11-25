/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
// Workflow to run GATK best practices and HaplotypeCaller
// before running you have to get the reference and and other tools
// installed on your system

// Extra notes:
// - The tiny dataset is the same as it is in used in the SciLifeLab/CAW project, its coordinates are -L "1:131941-141339"
//

import nextflow.extension.FilesEx

sampleTSVconfig = file(params.sample)

if (!params.sample) {
    text = Channel.from(
    "GATK-BP WGS WORKFLOW - Missing the sample TSV config file: ",
    "    Usage",
    "       nextflow run main.nf -c <file.config> --sample <sample.tsv>")
    text.subscribe { println "$it" }
    exit 1
}

def extractFastqFiles(tsvFile) {
/*
* Channeling the TSV file containing FASTQ
* The format is: "libPrep flowCell lane fastq1 fastq2"
*/
    fastqFiles = Channel
        .from(tsvFile.readLines())
        .map{ line ->
        list      = line.split()
        
        libPrep = list[0]
        flowCell = list[1]        
        lane = list[2]
        fastqFile1 = file(list[3])
        fastqFile2 = file(list[4])

        checkFileExistence(fastqFile1)
        checkFileExistence(fastqFile2)

      [ libPrep, flowCell, lane, fastqFile1, fastqFile2 ]
    }
  return fastqFiles
}

def checkFileExistence(fileToCheck) {
  try { assert file(fileToCheck).exists() }
  catch (AssertionError ae) {
    exit 1, "Missing file in TSV file: ${fileToCheck}, see --help for more information"
  }
}

referenceMap = [
  "genomeFile":  params.genome,       // genome reference
  "genomeIndex": params.genomeIndex,  // genome reference index
  "genomeDict":  params.genomeDict,   // genome reference dictionary
  "kgIndels":    params.kgIndels,     // 1000 Genomes SNPs
  "kgIndex":     params.kgIndex,      // 1000 Genomes SNPs index
  "dbsnp":       params.dbsnp,        // dbSNP
  "dbsnpIndex":  params.dbsnpIndex,   // dbSNP index
  "millsIndels": params.millsIndels,  // Mill's Golden set of SNPs
  "millsIndex":  params.millsIndex,   // Mill's Golden set index
  "cosmic41":    params.cosmic41,     // cosmic vcf file with VCF4.1 header
  "cosmic":      params.cosmic       // cosmic vcf file
]

fastqFiles = extractFastqFiles(file(params.sample))

process MapReads { 
    tag { params.projectID}
    time { params.runTime * task.attempt }
    publishDir "bwa_mem"

    input:
    set libPrep, flowCell, lane, file(fq1), file(fq2) from fastqFiles
    file referenceMap["genomeFile"]

    output:
    file("${libPrep}_${flowCell}_${lane}.bam") into bams

    script:
    readGroupString="\"@RG\\tID:${lane}\\tSM:${params.projectID}_${params.sampleID}\\tLB:${flowCell}_${libPrep}\\tPL:illumina\""
    """
    #!/bin/bash
    set -eo pipefail
    bwa mem -R ${readGroupString} -B 3 -t ${task.cpus} -M ${referenceMap["genomeFile"]} ${fq1} ${fq2} | \
    samtools view -bS -t ${referenceMap["genomeIndex"]} - | \
    samtools sort - > ${libPrep}_${flowCell}_${lane}.bam
    """
}

process MergeBAMs { 
    tag { params.projectID }
    publishDir "merged"

    time { params.runTime * task.attempt }

    input:
    file '*.bam' from bams.toList()

    output:
    file "${params.projectID}_${params.sampleID}.bam" into mergedBam

    script:
    """
    samtools merge ${params.projectID}_${params.sampleID}.bam *.bam 
    """
}

process MarkDuplicates { 
    tag {params.projectID}
    publishDir "md"

    input:
    file(merged) from mergedBam

    output:
    file "${params.projectID}_${params.sampleID}.md.bam" into mdBam
    file "${params.projectID}_${params.sampleID}.md.bai" into mdBamIdx
    file "${params.projectID}_${params.sampleID}.bam.metrics" into mdMetrics

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar ${params.picardHome}/MarkDuplicates.jar \
    INPUT=${merged} \
    METRICS_FILE=${merged}.metrics \
    TMP_DIR=. \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=TRUE \
    OUTPUT=${params.projectID}_${params.sampleID}.md.bam
    """
}

process CreateIntervals { 
    tag {params.projectID}
    publishDir "CreateIntervals"

    input:
    file md from mdBam
    file mdIdx from mdBamIdx
    file gf from file(referenceMap["genomeFile"])
    file gi from file(referenceMap["genomeIndex"])
    file gd from file(referenceMap["genomeDict"])
    file ki from file(referenceMap["kgIndels"])
    file kix from file(referenceMap["kgIndex"])
    file mi from file(referenceMap["millsIndels"])
    file mix from file(referenceMap["millsIndex"])

    output:
    file "${params.projectID}_${params.sampleID}.intervals" into intervals
    // We are actually duplicating the mdBam/Idx channels, but it have to be duplicated anyway
    file "${params.projectID}_${params.sampleID}.md.bam" into mdCIBam
    file "${params.projectID}_${params.sampleID}.md.bai" into mdCIBamIdx
    """
    java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I $md \
    -R $gf \
    -known $ki \
    -known $mi \
    -nt ${task.cpus} \
    -L "1:131941-141339" \
    -XL hs37d5 \
    -XL NC_007605 \
    -o ${params.projectID}_${params.sampleID}.intervals
    """
}

process IndelRealign {
    tag {params.projectID}
    publishDir "IndelRealigner"

    input:
    file cibam from mdCIBam
    file cibai from mdCIBamIdx
    file intervals from intervals
    file gf from file(referenceMap["genomeFile"])
    file gi from file(referenceMap["genomeIndex"])
    file gd from file(referenceMap["genomeDict"])
    file ki from file(referenceMap["kgIndels"])
    file kix from file(referenceMap["kgIndex"])
    file mi from file(referenceMap["millsIndels"])
    file mix from file(referenceMap["millsIndex"])

    output:
    file("*.md.real.bam") into realignedBam
    file("*.md.real.bai") into realignedBai

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -I $cibam \
    -R $gf \
    -targetIntervals $intervals \
    -known $ki \
    -known $mi \
    -XL hs37d5 \
    -XL NC_007605 \
    -L "1:131941-141339" \
    -o ${params.projectID}_${params.sampleID}.md.real.bam
    """     
}

process CreateRecalibrationTable { 
    tag {params.projectID}
    publishDir "CreateRecalibrationTable"

    input:
    file raBam from realignedBam
    file raBai from realignedBai
    file referenceMap["genomeFile"]
    file referenceMap["dbsnp"]
    file referenceMap["kgIndels"]
    file referenceMap["millsIndels"]

    output:
    file "${params.projectID}_${params.sampleID}.recal.table" into recalibrationTable
    file "${params.projectID}_${params.sampleID}.md.real.bam" into crRaBam
    file "${params.projectID}_${params.sampleID}.md.real.bai" into crRaBai
    //file "*.md.real.bam" into rcRealignedBam
    //file "*.md.real.bai" into rcRealignedBai

    script:
    """
    java -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir="/tmp" \
    -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${referenceMap["genomeFile"]} \
    -I $raBam \
    -knownSites ${referenceMap["dbsnp"]} \
    -knownSites ${referenceMap["kgIndels"]} \
    -knownSites ${referenceMap["millsIndels"]} \
    -nct ${task.cpus} \
    -XL hs37d5 \
    -XL NC_007605 \
    -l INFO \
    -L "1:131941-141339" \
    -o ${params.projectID}_${params.sampleID}.recal.table
    """
}

//process RecalibrateBam { 
//    input:
//    file realBambam), file(bai), recalibrationReport from recalibrationTable
//    file referenceMap["genomeFile"]
//
//  output:
//    set idPatient, gender, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam
//    set idPatient, gender, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV
//
//  script:
//  """
//  java -Xmx${task.memory.toGiga()}g \
//  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
//  -T PrintReads \
//  -R ${referenceMap["genomeFile"]} \
//  -nct ${task.cpus} \
//  -I $bam \
//  -XL hs37d5 \
//  -XL NC_007605 \
//  --BQSR $recalibrationReport \
//  -o ${idSample}.recal.bam
//  """}
//
//process HaplotypeCaller { }
//process RecalibrateVariants { }
//
