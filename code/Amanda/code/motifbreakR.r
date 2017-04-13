library("motifbreakR")
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(MotifDb)

snps = c("rs10791099", "rs1448360", "rs4937582", "rs10750450", "rs10791098")
snps = c("rs10791097","rs10791098","rs10791099","rs10750450")
snps = c("rs10894268","rs4366492","rs1944142","rs6590512","rs34529622",
         "rs10791097","rs10791098","rs10791099","rs10750450","rs67215852","rs4601795	rs3831404")

# load in SNPs
setwd("/dcl01/lieber/ajaffe/Amanda/ATAC/SNX19/")
snps.mb = snps.from.rsid(rsid = snps,
                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)
# find broken motifs
data(motifbreakR_motif)
motifbreakR_motif
results_hocomoco <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                                pwmList = hocomoco,
                                threshold = 1e-3,
                                method = "ic",
                                bkg = c(A=0.25, C=0.25, G=0.25, T=0.25) )
results_encode <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                              pwmList = encodemotif,
                              threshold = 1e-3,
                              method = "ic",
                              bkg = c(A=0.25, C=0.25, G=0.25, T=0.25) )
results_factorbook <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                                  pwmList = factorbook,
                                  threshold = 1e-3,
                                  method = "ic",
                                  bkg = c(A=0.25, C=0.25, G=0.25, T=0.25) )
results_homer <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                             pwmList = homer,
                             threshold = 1e-3,
                             method = "ic",
                             bkg = c(A=0.25, C=0.25, G=0.25, T=0.25) )
save(results_encode,results_hocomoco,results_factorbook,results_homer,file="./motifbreakR/results.0.001.rda")

# Calculate P Values
results_homer = calculatePvalue(results_homer)
save(results_homer, file="./motifbreakR/results_homer_0.001.rda")
results_factorbook = calculatePvalue(results_factorbook)
save(results_factorbook, file="./motifbreakR/results_factorbook_0.001.rda")
results_hocomoco = calculatePvalue(results_hocomoco)
save(results_hocomoco, file="./motifbreakR/results_hocomoco_0.001.rda")

rs10791099_encode <- results_encode[names(results_encode) %in% "rs10791099"]
rs10791099_encode <- calculatePvalue(rs10791099_encode)
save(rs10791099_encode, file="./motifbreakR/rs10791099_encode_0.001.rda")
rs10791098_encode <- results_encode[names(results_encode) %in% "rs10791098"]
rs10791098_encode <- calculatePvalue(rs10791098_encode)
save(rs10791098_encode, file="./motifbreakR/rs10791098_encode_0.001.rda")
rs10750450_encode <- results_encode[names(results_encode) %in% "rs10750450"]
rs10750450_encode <- calculatePvalue(rs10750450_encode)
save(rs10750450_encode, file="./motifbreakR/rs10750450_encode_0.001.rda")

rs10791097_encode <- results_encode[names(results_encode) %in% "rs10791097"]
rs10791097_encode1 = rs10791097_encode[1:35]
rs10791097_encode2 = rs10791097_encode[36:56]
rs10791097_encode3 = rs10791097_encode[57:67]
rs10791097_encode4 = rs10791097_encode[68:70]
rs10791097_encode5 = rs10791097_encode[71]
rs10791097_encode6 = rs10791097_encode[73:74]
rs10791097_encode7 = rs10791097_encode[72]
rs10791097_encode1 <- calculatePvalue(rs10791097_encode1)
save(rs10791097_encode1, file="./motifbreakR/rs10791097_encode1_0.001.rda")
rs10791097_encode2 <- calculatePvalue(rs10791097_encode2)
save(rs10791097_encode2, file="./motifbreakR/rs10791097_encode2_0.001.rda")
rs10791097_encode3 <- calculatePvalue(rs10791097_encode3)
save(rs10791097_encode3, file="./motifbreakR/rs10791097_encode3_0.001.rda")
rs10791097_encode4 <- calculatePvalue(rs10791097_encode4)
save(rs10791097_encode4, file="./motifbreakR/rs10791097_encode4_0.001.rda")
#rs10791097_encode5 <- calculatePvalue(rs10791097_encode5)
#save(rs10791097_encode5, file="./motifbreakR/rs10791097_encode5_0.001.rda")
rs10791097_encode6 <- calculatePvalue(rs10791097_encode6)
save(rs10791097_encode6, file="./motifbreakR/rs10791097_encode6_0.001.rda")
rs10791097_encode7 <- calculatePvalue(rs10791097_encode7)
save(rs10791097_encode7, file="./motifbreakR/rs10791097_encode7_0.001.rda")

# Combine results into table
load("./motifbreakR/rs10791097_encode1_0.001.rda")
load("./motifbreakR/rs10791097_encode2_0.001.rda")
load("./motifbreakR/rs10791097_encode3_0.001.rda")
load("./motifbreakR/rs10791097_encode4_0.001.rda")
load("./motifbreakR/rs10791097_encode6_0.001.rda")
load("./motifbreakR/rs10791097_encode7_0.001.rda")
rs10791097_encode1 <- as.data.frame(mcols(rs10791097_encode1))
rs10791097_encode2 <- as.data.frame(mcols(rs10791097_encode2))
rs10791097_encode3 <- as.data.frame(mcols(rs10791097_encode3))
rs10791097_encode4 <- as.data.frame(mcols(rs10791097_encode4))
rs10791097_encode5 <- as.data.frame(mcols(rs10791097_encode5))
rs10791097_encode6 <- as.data.frame(mcols(rs10791097_encode6))
rs10791097_encode7 <- as.data.frame(mcols(rs10791097_encode7))
rs10791097_encode = rbind(rs10791097_encode1,rs10791097_encode2,rs10791097_encode3,rs10791097_encode4,
                          rs10791097_encode5,rs10791097_encode6,rs10791097_encode7)
rs10791097_encode$SNP = "rs10791097"

load("./motifbreakR/rs10791098_encode_0.001.rda")
rs10791098_encode <- as.data.frame(mcols(rs10791098_encode))
rs10791098_encode$SNP <- "rs10791098"
load("./motifbreakR/rs10750450_encode_0.001.rda")
rs10750450_encode <- as.data.frame(mcols(rs10750450_encode))
rs10750450_encode$SNP <- "rs10750450"
load("./motifbreakR/rs10791099_encode_0.001.rda")
rs10791099_encode <- as.data.frame(mcols(rs10791099_encode))
rs10791099_encode$SNP <- "rs10791099"
res_encode = rbind(rs10750450_encode, rs10791098_encode, rs10791099_encode, rs10791097_encode)

load("./motifbreakR/results_homer_0.001.rda")
load("./motifbreakR/results_factorbook_0.001.rda")
load("./motifbreakR/results_hocomoco_0.001.rda")
res_hocomoco <- as.data.frame(mcols(results_hocomoco))
res_hocomoco$SNP <- names(results_hocomoco)
res_factorbook <- as.data.frame(mcols(results_factorbook))
res_factorbook$SNP <- names(results_factorbook)
res_homer <- as.data.frame(mcols(results_homer))
res_homer$SNP <- names(results_homer)
all_res <- rbind(res_hocomoco,res_encode,res_factorbook, res_homer)
all_res <- all_res[,c(ncol(all_res),1:(ncol(all_res)-1) )]
write.csv(all_res, file = "./motifbreakR/motifbreakR.output.pval.calc.0.001.csv", row.names=FALSE)





#extracting results for rs10791099
rs10791099_hocomoco <- results_hocomoco[names(results_hocomoco) %in% "rs10791099"]
rs10791099_hocomoco <- calculatePvalue(rs10791099_hocomoco)
rs10791099_hocomoco <- as.data.frame(mcols(rs10791099_hocomoco))
rs10791099_hocomoco$SNP <- "rs10791099"
rs10791099_homer <- results_homer[names(results_homer) %in% "rs10791099"]
rs10791099_homer <- calculatePvalue(rs10791099_homer)
rs10791099_homer <- as.data.frame(mcols(rs10791099_homer))
rs10791099_homer$SNP <- "rs10791099"
rs10791099_factorbook <- results_factorbook[names(results_factorbook) %in% "rs10791099"]
rs10791099_factorbook <- calculatePvalue(rs10791099_factorbook)
rs10791099_factorbook <- as.data.frame(mcols(rs10791099_factorbook))
rs10791099_factorbook$SNP <- "rs10791099"
rs10791099_encode <- results_encode[names(results_encode) %in% "rs10791099"]
length(rs10791099_encode)
rs10791099_encode
rs10791099_encode1 <- rs10791099_encode[1:60]
rs10791099_encode2 <- rs10791099_encode[61:120]
rs10791099_encode3 <- rs10791099_encode[121:180]
rs10791099_encode4 <- rs10791099_encode[181:240]
rs10791099_encode5 <- rs10791099_encode[241:287]
rs10791099_encode1 <- calculatePvalue(rs10791099_encode1)
rs10791099_encode2 <- calculatePvalue(rs10791099_encode2)
rs10791099_encode3 <- calculatePvalue(rs10791099_encode3)
rs10791099_encode4 <- calculatePvalue(rs10791099_encode4)
rs10791099_encode5 <- calculatePvalue(rs10791099_encode5)
rs10791099_encode1 <- as.data.frame(mcols(rs10791099_encode1))
rs10791099_encode2 <- as.data.frame(mcols(rs10791099_encode2))
rs10791099_encode3 <- as.data.frame(mcols(rs10791099_encode3))
rs10791099_encode4 <- as.data.frame(mcols(rs10791099_encode4))
rs10791099_encode5 <- as.data.frame(mcols(rs10791099_encode5))
rs10791099_encode = rbind(rs10791099_encode1,rs10791099_encode2,rs10791099_encode3,rs10791099_encode4,rs10791099_encode5)
rs10791099_encode$SNP <- "rs10791099"
rs10791099 = rbind(rs10791099_hocomoco, rs10791099_homer, rs10791099_factorbook, rs10791099_encode)
write.csv(rs10791099, file = "./MotifbreakR_rs10791099.csv", row.names=FALSE)

#extracting results for rs10750450
load("./results.0.01.rda")
rs10750450_hocomoco <- results_hocomoco[names(results_hocomoco) %in% "rs10750450"]
length(rs10750450_hocomoco)
rs10750450_hocomoco1 <- rs10750450_hocomoco[1:45]
rs10750450_hocomoco2 <- rs10750450_hocomoco[46:66]
rs10750450_hocomoco3 <- rs10750450_hocomoco[67:98]
rs10750450_hocomoco1 <- calculatePvalue(rs10750450_hocomoco1)
save(rs10750450_hocomoco1,file="./rs10750450_hocomoco1.rda")
rs10750450_hocomoco2 <- calculatePvalue(rs10750450_hocomoco2)
save(rs10750450_hocomoco2,file="./rs10750450_hocomoco2.rda")
rm("first_time","rs10750450_hocomoco","rs10750450_hocomoco1","snps","snps.mb","rs10750450_hocomoco2")
rs10750450_hocomoco3 <- calculatePvalue(rs10750450_hocomoco3)
save(rs10750450_hocomoco3,file="./rs10750450_hocomoco3.rda")
load("./rs10750450_hocomoco1.rda")
load("./rs10750450_hocomoco2.rda")
rs10750450_hocomoco1 <- as.data.frame(mcols(rs10750450_hocomoco1))
rs10750450_hocomoco2 <- as.data.frame(mcols(rs10750450_hocomoco2))
rs10750450_hocomoco3 <- as.data.frame(mcols(rs10750450_hocomoco3))

rs10750450_hocomoco = rbind(rs10750450_hocomoco1, rs10750450_hocomoco2,rs10750450_hocomoco3)
rs10750450_hocomoco$SNP <- "rs10750450"
write.csv(rs10750450_hocomoco, file = "./rs10750450_hocomoco.csv", row.names=FALSE)

rs10750450_encode <- results_encode[names(results_encode) %in% "rs10750450"]
length(rs10750450_encode)
rs10750450_encode1 <- rs10750450_encode[1:20]
rs10750450_encode2 <- rs10750450_encode[21:40]
rs10750450_encode3 <- rs10750450_encode[41:60]
rs10750450_encode4 <- rs10750450_encode[61:80]
rs10750450_encode5 <- rs10750450_encode[81:100]
rs10750450_encode6 <- rs10750450_encode[101:120]
rs10750450_encode7 <- rs10750450_encode[121:140]
rs10750450_encode8 <- rs10750450_encode[141:160]
rs10750450_encode9 <- rs10750450_encode[161:180]
rs10750450_encode10 <- rs10750450_encode[181:200]
rs10750450_encode11 <- rs10750450_encode[201:220]
rs10750450_encode12 <- rs10750450_encode[221:240]
rs10750450_encode13 <- rs10750450_encode[241:260]
rs10750450_encode14 <- rs10750450_encode[261:280]
rs10750450_encode15 <- rs10750450_encode[281:300]
rs10750450_encode16 <- rs10750450_encode[301:320]
rs10750450_encode17 <- rs10750450_encode[321:340]
rs10750450_encode18 <- rs10750450_encode[341:367]

rs10750450_encode1 <- calculatePvalue(rs10750450_encode1)
rs10750450_encode2 <- calculatePvalue(rs10750450_encode2)
rs10750450_encode3 <- calculatePvalue(rs10750450_encode3)
rs10750450_encode4 <- calculatePvalue(rs10750450_encode4)
rs10750450_encode5 <- calculatePvalue(rs10750450_encode5)
rs10750450_encode6 <- calculatePvalue(rs10750450_encode6)
rs10750450_encode7 <- calculatePvalue(rs10750450_encode7)
rs10750450_encode8 <- calculatePvalue(rs10750450_encode8)
rs10750450_encode9 <- calculatePvalue(rs10750450_encode9)
rs10750450_encode10 <- calculatePvalue(rs10750450_encode10)
rs10750450_encode11 <- calculatePvalue(rs10750450_encode11)
rs10750450_encode12 <- calculatePvalue(rs10750450_encode12)
rs10750450_encode13 <- calculatePvalue(rs10750450_encode13)
rs10750450_encode14 <- calculatePvalue(rs10750450_encode14)
rs10750450_encode15 <- calculatePvalue(rs10750450_encode15)
rs10750450_encode16 <- calculatePvalue(rs10750450_encode16)
rs10750450_encode17 <- calculatePvalue(rs10750450_encode17)
rs10750450_encode18 <- calculatePvalue(rs10750450_encode18)

rs10750450_encode1 <- as.data.frame(mcols(rs10750450_encode1))
rs10750450_encode2 <- as.data.frame(mcols(rs10750450_encode2))
rs10750450_encode3 <- as.data.frame(mcols(rs10750450_encode3))
rs10750450_encode4 <- as.data.frame(mcols(rs10750450_encode4))
rs10750450_encode5 <- as.data.frame(mcols(rs10750450_encode5))
rs10750450_encode6 <- as.data.frame(mcols(rs10750450_encode6))
rs10750450_encode = rbind(rs10750450_encode1,rs10750450_encode2,rs10750450_encode3,rs10750450_encode4,rs10750450_encode5,rs10750450_encode6)
rs10750450_encode$SNP <- "rs10750450"
write.csv(rs10750450_encode, file = "./rs10750450_encode.csv", row.names=FALSE)

rs10750450_factorbook <- results_factorbook[names(results_factorbook) %in% "rs10750450"]
length(rs10750450_factorbook)
rs10750450_factorbook <- calculatePvalue(rs10750450_factorbook)
rs10750450_factorbook <- as.data.frame(mcols(rs10750450_factorbook))
rs10750450_factorbook$SNP <- "rs10750450"
write.csv(rs10750450_factorbook, file = "./rs10750450_factorbook.csv", row.names=FALSE)

rs10750450_homer <- results_homer[names(results_homer) %in% "rs10750450"]
length(rs10750450_homer)
rs10750450_homer <- calculatePvalue(rs10750450_homer)
rs10750450_homer <- as.data.frame(mcols(rs10750450_homer))
rs10750450_homer$SNP <- "rs10750450"
write.csv(rs10750450_homer, file = "./rs10750450_homer.csv", row.names=FALSE)

rs10750450 = rbind(rs10750450_hocomoco, rs10750450_homer, rs10750450_factorbook, rs10750450_encode)
write.csv(rs10750450, file = "./MotifbreakR_rs10750450.csv", row.names=FALSE)

#extracting results for rs10791098
rs10791098_hocomoco <- results_hocomoco[names(results_hocomoco) %in% "rs10791098"]
length(rs10791098_hocomoco)
rs10791098_hocomoco <- calculatePvalue(rs10791098_hocomoco)
rs10791098_hocomoco <- as.data.frame(mcols(rs10791098_hocomoco))
rs10791098_hocomoco$SNP <- "rs10791098"
write.csv(rs10791098_hocomoco, file = "./rs10791098_hocomoco.csv", row.names=FALSE)

rs10791098_encode <- results_encode[names(results_encode) %in% "rs10791098"]
length(rs10791098_encode)
rs10791098_encode1 <- rs10791098_encode[1:80]
rs10791098_encode2 <- rs10791098_encode[81:160]
rs10791098_encode3 <- rs10791098_encode[161:240]
rs10791098_encode4 <- rs10791098_encode[241:320]
rs10791098_encode5 <- rs10791098_encode[321:383]
rs10791098_encode1 <- calculatePvalue(rs10791098_encode1)
rs10791098_encode2 <- calculatePvalue(rs10791098_encode2)
rs10791098_encode3 <- calculatePvalue(rs10791098_encode3)
rs10791098_encode4 <- calculatePvalue(rs10791098_encode4)
rs10791098_encode5 <- calculatePvalue(rs10791098_encode5)
rs10791098_encode1 <- as.data.frame(mcols(rs10791098_encode1))
rs10791098_encode2 <- as.data.frame(mcols(rs10791098_encode2))
rs10791098_encode3 <- as.data.frame(mcols(rs10791098_encode3))
rs10791098_encode4 <- as.data.frame(mcols(rs10791098_encode4))
rs10791098_encode5 <- as.data.frame(mcols(rs10791098_encode5))
rs10791098_encode = rbind(rs10791098_encode1,rs10791098_encode2,rs10791098_encode3,rs10791098_encode4,rs10791098_encode5)
rs10791098_encode$SNP <- "rs10791098"
write.csv(rs10791098_encode, file = "./rs10791098_encode.csv", row.names=FALSE)

rs10791098_factorbook <- results_factorbook[names(results_factorbook) %in% "rs10791098"]
length(rs10791098_factorbook)
rs10791098_factorbook <- calculatePvalue(rs10791098_factorbook)
rs10791098_factorbook <- as.data.frame(mcols(rs10791098_factorbook))
rs10791098_factorbook$SNP <- "rs10791098"
write.csv(rs10791098_factorbook, file = "./rs10791098_factorbook.csv", row.names=FALSE)

rs10791098_homer <- results_homer[names(results_homer) %in% "rs10791098"]
length(rs10791098_homer)
rs10791098_homer <- calculatePvalue(rs10791098_homer)
rs10791098_homer <- as.data.frame(mcols(rs10791098_homer))
rs10791098_homer$SNP <- "rs10791098"
write.csv(rs10791098_homer, file = "./rs10791098_homer.csv", row.names=FALSE)

rs10791098 = rbind(rs10791098_hocomoco, rs10791098_homer, rs10791098_factorbook, rs10791098_encode)
write.csv(rs10791098, file = "./MotifbreakR_rs10791098.csv", row.names=FALSE)