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

fastqFiles = extractFastqFiles(file(params.sample))
//process BWAMap { }
//process MergeBAMs { }
//process MarkDuplicates { }
//process CreateIntervals { }
//process IndelRealign { }
//process CreateRecalibrationTable { }
//process RecalibrateBam { }
//process HaplotypeCaller { }
//process RecalibrateVariants { }
//
