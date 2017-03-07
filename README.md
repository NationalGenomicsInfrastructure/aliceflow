# aliceflow
NextFlow framework to streamline the Genalice variant call pipeline

Having two sets of calls, A and B:

calls only in B: 

    vcfintersect -r ref.fasta -v -i A.vcf B.vcf

calls only in A: 

    vcfintersect -r ref.fasta -v -i B.vcf A.vcf

calls both in A and B (intersect): 

    vcfintersect -r ref.fasta -i A.vcf B.vcf


union calls: 

    vcfintersect -r ref.fasta -u A.vcf B.vcf

    ps -eo pmem,pid,pcpu,rss,vsz,time,args | sort -k 1 -r| less -S

    perl -pi -e 's/chr//' PL.vcf

Generated calls that ar in NIST, but not in GATK

    vcfintersect -r ~/genome/human_g1k_v37_decoy.fasta -v -i PL.vcf NIST.vcf > NIST_not_GATK.vcf

then calls that are in NIST, and in GA (and not in GATK):

    vcfintersect -r ~/genome/human_g1k_v37_decoy.fasta -i NIST_not_GATK.vcf GA.vcf > GA_NIST_not_in_GATK.vcf

calls that are in GA only, neither in GATK or NIST:

    vcfintersect -r ~/genome/human_g1k_v37_decoy.fasta -u NIST.vcf PL.vcf > NIST_U_Plat.vcf # get the union (KILLED, can't do that)

calls in GA but not in NIST

    vcfintersect -r ~/genome/human_g1k_v37_decoy.fasta -v -i NIST.vcf GA.vcf > GA_compl_NIST.vcf

#n of calls:
   3642054 NIST.vcf
   3953641 PL.vcf
   4564558 GA.vcf
    172054 NIST_not_GATK.vcf
    141543 GA_NIST_not_in_GATK.vcf

