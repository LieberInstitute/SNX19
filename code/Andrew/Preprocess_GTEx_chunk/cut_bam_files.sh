#!/bin/sh

# qsub -t 1-8208 -tc 150 -l mf=3G,h_vmem=4G,h_stack=256M cut_bam_files.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/SNX19/GTEX/GTEX_Samples_SNX19_analysis.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/TopHat/$ID/accepted_hits.bam
CUTBAM=/dcl01/lieber/ajaffe/Brain/SNX19/GTEX/BAM/${ID}_SNX19.bam

# make samtools range: chr11:129718630-131718630
r=chr11:129718630-131718630
samtools view -b $BAM $r > $CUTBAM
samtools index $CUTBAM

####### run feature counts ######
a=/dcl01/lieber/ajaffe/Brain/SNX19/SNX19_Ensembl_v75_chrPrefix.gtf

### GENE
oGene=/dcl01/lieber/ajaffe/Brain/SNX19/GTEX/Counts/Gene/${ID}_SNX19_Ensembl75_genes.counts
featureCounts -a $a -o $oGene $CUTBAM

### EXON
oExon=/dcl01/lieber/ajaffe/Brain/SNX19/GTEX/Counts/Exon/${ID}_SNX19_Ensembl75_exons.counts
featureCounts -O -f -a $a -o $oExon $CUTBAM
