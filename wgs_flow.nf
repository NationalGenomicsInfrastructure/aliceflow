// Test workflow to run a simple run with the genalice server
// before running you have to get a reference and start the server:
// gaReference --input=/home/genome/build/GRCh37/GRCh37.fa --output=/store/nf/reference/GRCh37.grf
// gaServer --server_name=DEMO --reference=/store/nf/reference/GRCh37.grf

//fq1=file("/store/G15511_S1/tiny_1.fastq")
//fq2=file("/store/G15511_S1/tiny_2.fastq")

//fq1=file("/store/G15511_S1/C09DFACXX111207.1.TTGAGCCT_1.fastq.gz")
//fq2=file("/store/G15511_S1/C09DFACXX111207.1.TTGAGCCT_2.fastq.gz")


sampleTSVconfig = file(params.sample)

if (!params.sample) {
	text = Channel.from(
	"GENALICE WGS WORKFLOW - Missing the sample TSV config file: ",
	"    Usage",
	"       nextflow run wgs_flow.nf -c <file.config> --sample <sample.tsv>")
  	text.subscribe { println "$it" }
  	exit 1
}

rawFastqFiles = Channel
  .from(sampleTSVconfig.readLines())
  .map { line ->
    list        = line.split()
    idSample	= list[0]
    idType      = list[2]
    fastqFile1  = file(list[2])
    fastqFile2  = file(list[3])
    [ idSample, idType, fastqFile1, fastqFile2 ]
}

/*
	To considerread groups, Genalice mapping have to be set up in a certain way.
	Advice from Bas Tolhuis:

We have data set from two lanes (L001 and L002).
myData-L001.R1.fastq.gz
myData-L001.R2.fastq.gz
myData-L002.R1.fastq.gz
myData-L002.R2.fastq.gz

We create two directories:
mkdir myData.pair1
mkdir myData.pair2

We make symbolic links of the R1 files in the myData.pair1 folder
ln -s /full/path/to/myData-L001.R1.fastq.gz myData.pair1/myData-L001.R1.fastq.gz
ln -s /full/path/to/myData-L002.R1.fastq.gz myData.pair1/myData-L002.R1.fastq.gz

We make symbolic links of the R2 files in the myData.pair2 folder
ln -s /full/path/to/myData-L001.R2.fastq.gz myData.pair2/myData-L001.R2.fastq.gz
ln -s /full/path/to/myData-L002.R2.fastq.gz myData.pair2/myData-L002.R2.fastq.gz

We give directories as input to gaMap
Note that its a comma separated list and each directory ends with a "/" character
gaMap --server_name=myServer --input=myData.pair1/,myData.pair2/ --output=... etc

Now the fastq links in the myData.pair* directories will be mapped into one GAR file 
using alphabetical order. The order in pair1 and pair2 directories should be the same 
(i.e. first L001 and then L002).. It is not possible to give the individual fastq inputs 
different read group names. 
 */

(rawFastqFiles, idSample) = getIdSample(rawFastqFiles)

println "Sample ID: " + idSample

R1_dir = file("${idSample}.R1")
R2_dir = file("${idSample}.R2")
R1_dir.mkdir()
R2_dir.mkdir()


// create files containing full path names for R1 and R2 reads

R1File = file("R1.readgroups")
R2File = file("R2.readgroups")
R1Channel = Channel.create()
R2Channel = Channel.create()
Channel.from rawFastqFiles.separate(R1Channel,R2Channel) {x -> [x,x] }

// Now we have two channels for reads, make symlinks in the corresponding directories
// TODO: these lines are not checking for existence and are throwing an error if 
// files are already there
R1FileNames = R1Channel.map { x -> x.get(2)}
	.subscribe { it -> file(it).mklink(idSample +".R1/" + it.fileName) }
R2FileNames = R2Channel.map { x -> x.get(3)}
	.subscribe { it -> file(it).mklink(idSample +".R2/" + it.fileName) }

//process map {
//	input: 
//	file fq1 
//	file fq2
//
//	output: 
//	file "G15511_S1.gar" into result_gar
//	file "G15511_S1.report" into result_report
//
//	"""
//	gaMap --server_name=DEMO --input=${fq1},${fq2} --output=`pwd`/G15511_S1.gar --run_name=G15511_1 --report_file=`pwd`/G15511_S1.report --cmd_file=/store/params/human.map.conf
//	"""
//}
//process variant_call {
//	input:
//	file result_gar
//
//	output:
//	file "G15511_S1.vcf" into VCF
//	file "G15511_S1.variant.json" into json
//
//	"""
//	gaVariant --input=$result_gar --output=`pwd`/G15511_S1.vcf --run_name=G15511_1 --cmd_file=/store/params/human.variant.conf --output_statistics=true --statistics_file=`pwd`/G15511_S1.variant.json 
//	"""
//}
//
def getIdSample(aCh) {

    consCh = Channel.create()
    originalCh = Channel.create()

    // get the patient ID
    // duplicate channel to get sample name
    Channel.from aCh.separate(consCh,originalCh) {x -> [x,x]}

    // use the "consumed" channel to get it
    // we are assuming the first column is the same for the patient, as hoping
    // people do not want to compare samples from differnet patients
    idPatient = consCh.map { x -> [x.get(0)]}.unique().getVal()[0]
    // we have to close to make sure remainding items are not
    consCh.close()

    return [ originalCh, idPatient]
}
