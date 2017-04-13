library('BSgenome.Hsapiens.UCSC.hg19')
library(readxl)

table = read_excel("./Dropbox/SNX19/code/Amanda/supporting_files/SCZ52_updated_pgc_vals_marquee_junction_eqtl_040517_CUR.xlsx", sheet=2)
table = data.frame(table)

alleles = strsplit(table$observed, "/", fixed=T)
allele1 = allele2 = c()
for (i in 1:length(alleles)){allele1[i] = alleles[[i]][2]}
for (i in 1:length(alleles)){allele2[i] = alleles[[i]][1]}
pos = table$snpPos
chr = 'chr11'
offset = 10

RefSeq = paste0(">Extended_Hap_", table$name,"\n",
                getSeq(Hsapiens,chr,pos-offset,pos-1),
                allele1,
                getSeq(Hsapiens,chr,pos+1,pos+offset))
AltSeq = paste0(">Extended_Hap_",table$name,"\n",
                getSeq(Hsapiens,chr,pos-offset,pos-1),
                allele2,
                getSeq(Hsapiens,chr,pos+1,pos+offset))
write.table(RefSeq, "./Dropbox/SNX19/code/Amanda/supporting_files/ref_extended_Hap.fa", quote = F, row.names=F, col.names=F)
write.table(AltSeq, "./Dropbox/SNX19/code/Amanda/supporting_files/alt_extended_Hap.fa", quote = F, row.names=F, col.names=F)
