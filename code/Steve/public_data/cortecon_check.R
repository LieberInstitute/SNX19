## Check correlation of junc8.10
load('/dcl01/ajaffe/data/lab/brainseq_phase1/preprocessed_data/rpkmCounts_DLPFC_polyA_BrainSeq_Phase1_LIBD_n959.rda')
pheno = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseI/LIBD_RNA-seq_ALL_jhs_04_09_2014_toALL_750.csv",as.is=TRUE)
pheno = pheno[,1:33]
pheno$RNum = paste0("R", pheno$RNum)
pheno$BRNum = paste0("Br", pheno$BRNum)
pheno = pheno[,1:8]
pd = cbind(pheno[match(metrics$RNum, pheno$RNum),], metrics[,-ncol(metrics)])
names(pd)[c(2,5)] = c("BrNum", "Age")
jaffelab::ss(colnames(jRpkm),"\\_")
jRpkm_hg38 = jRpkm
jMap_hg38=jMap
## load dlpfc polya
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")

##
x = data.table::fread('/dcl01/lieber/ajaffe/Steve/SNX19/cortecon/GSE31845_raw_and_processed_SNP.txt')
x1 = subset(x, Chr=="11")
x1[x1$Position> 130718380-10000 & x1$Position< 130718380+10000,]
chr11:130718380-130718880

rsnps::ld_search

rs10894376
rs10894377
rs11606643 
rs12288364 
rs2187574
rs2511775  
rs2513532  
rs2513533 
rs4300402  
rs7938687
rs7947194  
rs7952557  
rs7952585