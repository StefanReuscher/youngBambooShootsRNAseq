#!/bin/sh
#$ -S /bin/sh
#$ -N TISSUE_tophat2_PH_refseq_maskScaf
fwdReads=$(ls -d ~/bambooRNAseq/TISSUE_fwd_P.fastq)
revReads=$(ls -d ~/bambooRNAseq/TISSUE_rev_P.fastq)
singletons=$(ls -dm ~/bambooRNAseq/TISSUE_*_S.fastq | tr -d "\n")
export PATH=/usr/local/pkg/bowtie2/2.2.8/:$PATH
/usr/local/pkg/tophat2/2.1.0/tophat2 -o tophat2_TISSUE ~/ref/bowtie2Index/PH_refseq_maskScaf $fwdReads $revReads $singletons
