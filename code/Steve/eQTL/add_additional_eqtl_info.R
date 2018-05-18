### Add additional eqtl infos
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.rda')

### Load exon ensembl IDs
load('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/exon_names_hg19_hg38_deduplicated.rda')
hg19_exons$chrpos = paste0(hg19_exons$Chr, ":", hg19_exons$hg19_Start, "-", hg19_exons$hg19_End )
table(duplicated(hg19_exons$chrpos) )

cisMerge$gencode_exonID = gsub("\\_.*","",hg19_exons[match(cisMerge[,'ucsc_feature'], hg19_exons$chrpos),'gencode_exonID'])
cisMerge$gencode_exonID[cisMerge$exprsType!="Exon"] = NA #Set non-exonic features to NA
cisMerge$hg19_eID = hg19_exons[match(cisMerge[,'ucsc_feature'], hg19_exons$chrpos),'hg19_eID']
cisMerge$hg19_eID[cisMerge$exprsType!="Exon"] = NA #Set non-exonic features to NA
save(cisMerge,file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
write.csv(cisMerge,file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/tables/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.csv',row.names=F)

#QC
table(is.na(cisMerge$gencode_exonID), cisMerge$exprsType)
table(is.na(cisMerge$hg19_eID), cisMerge$exprsType)

##
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
minus_log_PVals = -log10(cisMerge[,grep("pvalue", colnames(cisMerge) )])
colnames(minus_log_PVals) = paste0("-log10_",colnames(minus_log_PVals))
cisMerge = cbind(cisMerge, minus_log_PVals)
write.csv(cisMerge, '/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/tables/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.csv',row.names=F)