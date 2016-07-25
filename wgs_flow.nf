// Test workflow to run a simple run with the genalice server
// before running you have to get a reference and start the server:
// gaReference --input=/home/genome/build/GRCh37/GRCh37.fa --output=/store/nf/reference/GRCh37.grf
// gaServer --server_name=DEMO --reference=/store/nf/reference/GRCh37.grf

//fq1=file("/store/G15511_S1/tiny_1.fastq")
//fq2=file("/store/G15511_S1/tiny_2.fastq")

fq1=file("/store/G15511_S1/C09DFACXX111207.1.TTGAGCCT_1.fastq.gz")
fq2=file("/store/G15511_S1/C09DFACXX111207.1.TTGAGCCT_2.fastq.gz")

process map {
	input: 
	file fq1 
	file fq2

	output: 
	file "G15511_S1.gar" into result_gar
	file "G15511_S1.report" into result_report

	"""
	gaMap --server_name=DEMO --input=${fq1},${fq2} --output=`pwd`/G15511_S1.gar --run_name=G15511_1 --report_file=`pwd`/G15511_S1.report --cmd_file=/store/params/human.map.conf
	"""
}
process variant_call {
	input:
	file result_gar

	output:
	file "G15511_S1.vcf" into VCF
	file "G15511_S1.variant.json" into json

	"""
	gaVariant --input=$result_gar --output=`pwd`/G15511_S1.vcf --run_name=G15511_1 --cmd_file=/store/params/human.variant.conf --output_statistics=true --statistics_file=`pwd`/G15511_S1.variant.json 
	"""
}
