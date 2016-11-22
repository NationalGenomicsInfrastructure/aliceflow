/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/
// Workflow to run GATK best practices and HaplotypeCaller
// before running you have to get the reference and and other tools
// installed on your system

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

refs = [
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

    input:
    set libPrep, flowCell, lane, file(fq1), file(fq2) from fastqFiles
    file refs["genomeFile"]

    output:
    set libPrep, flowCell, lane, file("${libPrep}_${flowCell}_${lane}.bam") into bams

    script:
    readGroupString="\"@RG\\tID:${lane}\\tSM:${params.projectID}_${params.sampleID}\\tLB:${flowCell}_${libPrep}\\tPL:illumina\""
    """
    #!/bin/bash
    set -eo pipefail
    bwa mem -R ${readGroupString} -B 3 -t ${task.cpus} -M ${refs["genomeFile"]} ${fq1} ${fq2} | \
    samtools view -bS -t ${refs["genomeIndex"]} - | \
    samtools sort - > ${libPrep}_${flowCell}_${lane}.bam
    """
}

process MergeBAMs { 
  tag { params.projectID }

  //queue 'core'
  time { params.runTime * task.attempt }

  input:
  set libPrep, flowCell, lane, file(bam) from bams

  output:
  set libPrep, flowCell, lane, file("${params.projectID}_${params.sampleID}.bam") into mergedBam

  script:
  """
  #!/bin/bash
  java -Xmx6g -jar ${params.picardHome}/MergeSamFiles.jar \
  ASSUME_SORTED=true \
  INPUT=${bam} \
  OUTPUT=${params.projectID}_${params.sampleID}.bam



#  samtools merge ${params.projectID}_${params.sampleID}.bam ${bam}
  """
}

//process MarkDuplicates { }
//process CreateIntervals { }
//process IndelRealign { }
//process CreateRecalibrationTable { }
//process RecalibrateBam { }
//process HaplotypeCaller { }
//process RecalibrateVariants { }
//
