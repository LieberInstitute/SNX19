#!/bin/sh

# qsub -t 1-738 -tc 100 -l mf=4G,h_vmem=5G,h_stack=256M cut_bam_files_DLPFC.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/SNX19/LIBD_Samples_SNX19_analysis_DLPFC.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam
CUTBAM=/dcl01/lieber/ajaffe/Brain/SNX19/BAM/DLPFC_PolyA_${ID}_SNX19.bam

# make samtools range: chr11:129718630-131718630
r=chr11:129718630-131718630
samtools view -b $BAM $r > $CUTBAM
samtools index $CUTBAM

####### run feature counts ######
a=/dcl01/lieber/ajaffe/Brain/SNX19/SNX19_Ensembl_v75_chrPrefix.gtf

oGene=/dcl01/lieber/ajaffe/Brain/SNX19/Counts/Gene/${ID}_SNX19_Ensembl75_genes_DLPFC.counts
if [ ! -e $oGene ] ; then

	mkdir -p /dcl01/lieber/ajaffe/Brain/SNX19/Counts

	### GENE
	mkdir -p /dcl01/lieber/ajaffe/Brain/SNX19/Counts/Gene

	featureCounts -a $a -o $oGene $CUTBAM

	### EXON
	mkdir -p /dcl01/lieber/ajaffe/Brain/SNX19/Counts/Exon
	oExon=/dcl01/lieber/ajaffe/Brain/SNX19/Counts/Exon/${ID}_SNX19_Ensembl75_exons_DLPFC.counts
	featureCounts -O -f -a $a -o $oExon $CUTBAM
fi 