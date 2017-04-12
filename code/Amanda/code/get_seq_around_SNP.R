library('BSgenome.Hsapiens.UCSC.hg19')
library(readxl)

table = read_excel("./Dropbox/SNX19/Junc8_10_SNX19_eQTL_results_updatedSNPs.xlsx")
table = table[which(table$Set!="na"),]
dim(table)
table = table[,c("updated_newRef","updated_newCount","updated_SNP","snpChrPos","UniqueID",
                 "snps","All_pvalue","PGC2014_P" ,"PGC2016_P","MAF_All","Set")]
table = table[which(table$updated_SNP!="rs67215852"),]
table[,c("updated_SNP","snpChrPos")]
ids = as.numeric(gsub("chr11:", "", table$snpChrPos, fixed=T))
ids

chr = 'chr11'
ref = table$updated_newRef
alt = table$updated_newCount
offset = 15

RefSeq = paste0(">Set",table$Set,"_", table$updated_SNP,"\n",
                getSeq(Hsapiens,chr,ids-offset,ids-1),
                ref,
                getSeq(Hsapiens,chr,ids+1,ids+offset))
AltSeq = paste0(">Set",table$Set,"_", table$updated_SNP,"\n",
                getSeq(Hsapiens,chr,ids-offset,ids-1),
                alt,
                getSeq(Hsapiens,chr,ids+1,ids+offset))
write.table(RefSeq, "./Dropbox/SNX19/ref.fa", quote = F, row.names=F, col.names=F)
write.table(AltSeq, "./Dropbox/SNX19/alt.fa", quote = F, row.names=F, col.names=F)
