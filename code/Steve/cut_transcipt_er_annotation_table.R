# ERs + Tx
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")

regions_SNX19=as.data.frame(regions[which(regions$nearestSymbol=="SNX19")])
regions_SNX19$ER_id=rownames(regions_SNX19)
tMap_SNX19 = as.data.frame(tMap[which(tMap$gene_name=="SNX19")])


mapWrite=list(er_map_snx19=regions_SNX19,tx_map_snx19=tMap_SNX19)

openxlsx::write.xlsx(mapWrite,file='/dcl01/lieber/ajaffe/Steve/SNX19/map_objects/er_tx_SNX19_map.xlsx')

### Cut PGC SNP Objects
pgc69 = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/daner_pgcEU69_pgcAS.FOR_IPSYC.FOR_run2.gz.p4.clump.areator.sorted.1mhc.csv')
pgc69=pgc69[,1:9]
pgc69=tidyr::separate(pgc69,A1A2,c('A1','A2'), sep="/") 

pgc69_subset = pgc69[pgc69$CHR=="11" & pgc69$BP<=130718630+1e6 & pgc69$BP >130718630-1e6,]

write.csv(pgc69,file='/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/FUMA_cut/daner_pgcEU69_pgcAS.FOR_IPSYC.FOR_run2.gz.p4.clump.areator.sorted.1mhc_2mb_Cut.csv',row.names=FALSE)
### 
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/SCZ52_may13_9444230snps-scz2-snp-results-ckqny-scz2snpres-rall.rda')
pgc52_subset = PGC_52[PGC_52$CHR=="chr11" & PGC_52$BP<=130718630+1e6 & PGC_52$BP >130718630-1e6,]
write.csv(pgc52_subset,file='/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/FUMA_cut/SCZ52_may13_9444230snps-scz2-snp-results-ckqny-scz2snpres-rall_2mb_Cut_all_PGC_SNPs.csv',row.names=FALSE)

### Write as .txt files
#write.csv(pgc69,file='/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/FUMA_cut/daner_pgcEU69_pgcAS.FOR_IPSYC.FOR_run2.gz.p4.clump.areator.sorted.1mhc_2mb_Cut.csv',row.names=FALSE)
