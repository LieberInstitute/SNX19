
##########################
#Author: Liang Ma, Qiang Chen
#SNX19
#conditional analysis for junc8.10 and junc8.8a
##########################

library(MatrixEQTL)
library(GenomicRanges)
library(sva)
library(XVector)
library(Biostrings)


# load data
#--------------------------------------------------------------------------------------
load("rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495.rda")
load("SNX19_DNAm_data_DLPFC.rda");





rownames(snpAll) = paste("snp_", rownames(snpAll), sep="");
rownames(methMat) = paste("meth_", rownames(methMat), sep="");


snpAll = data.frame(t(snpAll));
snpAll$BrNum = row.names(snpAll);

meths = t(methMat);
meths = data.frame(meths);
meths$BrNum = row.names(meths);


meths_pheno = merge(meths, methPheno, by="BrNum");


meths_snp = merge(meths_pheno, snpAll, by="BrNum");



## Screen samples
##------------------------------------------------------------------------------
aIndex = which(meths_snp$Age>13)     
methMat2 = as.matrix(meths_snp[aIndex, grep("meth_", names(meths_snp))])
snp2 = as.matrix(meths_snp[aIndex, grep("snp_", names(meths_snp))]) 
meths_snp = meths_snp[aIndex,]



## PCA
#--------------------------------------------------------------------------------------
pca = prcomp(methMat2)


meths_snp_pca = cbind(meths_snp, exprPC)


exprpc_all = meths_snp_pca[, grep("^PC", colnames(meths_snp_pca))];

exprpc_all = as.matrix(exprpc_all) 

snppc_all = meths_snp_pca[, grep("^snpPC", colnames(meths_snp_pca))];


snppc_all = as.matrix(snppc_all) 

methMat3 = as.matrix(meths_snp_pca[grep("meth_", names(meths_snp_pca))])
snp3 = as.matrix(meths_snp_pca[grep("snp_", names(meths_snp_pca))])


## covariates 
#--------------------------
mod = model.matrix(~ snppc_all + exprpc_all)
colnames(mod)[2:ncol(mod)] = c(paste0("snpPC",1:5), paste0("exprsPC",1:15))


covs = SlicedData$new(t(mod[,-1])) 

# snp position
#--------------------------
snpspos = snpMapAll[,c("SNP","CHR","POS")] 
snpspos$SNP = gsub(":", "\\.", snpspos$SNP); 
snpspos$CHR = paste0("chr",snpspos$CHR) 
colnames(snpspos) = c("snpid","chr","pos")
snpspos[,1] = paste("snp_", snpspos[,1], sep=""); 


# meth position
#--------------------------
methpos = cbind(rownames(methMap), methMap[,c("chr", "pos", "pos")]);
names(methpos) = c("geneid", "chr", "left", "right");
methpos[,1] = paste("meth_", methpos[,1], sep="");
methpos[1:5,]
dim(methpos) 

#snpSlice
#-----------------
snpSlice = SlicedData$new(t(snp2)) 
snpSlice$ResliceCombined(sliceSize = 5000)



#exprsSlice
#-----------------
exprsSlice = SlicedData$new(t(methMat2)) 
exprsSlice$ResliceCombined(sliceSize = 5000)


#all adult
#-----------------------------------------
meJoint = Matrix_eQTL_main(snps=snpSlice, gene = exprsSlice, 
                           cvrt = covs, output_file_name.cis =  ".ctxt" ,
                           pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                           snpspos = snpspos, genepos = methpos, 
                           useModel = modelLINEAR,	cisDist=5e6,
                           pvalue.hist = 100, min.pv.by.genesnp = T)	



# all CAUC 
#-------------------------------
caucIndex=which(meths_snp_pca$Race=="CAUC")


snpSliceCauc = SlicedData$new(t(snp3[caucIndex,]))
snpSliceCauc$ResliceCombined(sliceSize = 5000)

exprsCauc = SlicedData$new(t(methMat3[caucIndex,])) 
exprsCauc$ResliceCombined(sliceSize = 5000)



modCauc = mod[caucIndex,]

covsCauc = SlicedData$new(t(modCauc[,-1]))

meJointCauc = Matrix_eQTL_main(snps=snpSliceCauc, gene = exprsCauc, 
                               cvrt = covsCauc, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all aa 
#-------------------------------
aaIndex=which(meths_snp_pca$Race=="AA")


snpSliceaa = SlicedData$new(t(snp3[aaIndex,]))
snpSliceaa$ResliceCombined(sliceSize = 5000)

exprsaa = SlicedData$new(t(methMat3[aaIndex,])) 
exprsaa$ResliceCombined(sliceSize = 5000)

modaa = mod[aaIndex,]

covsaa = SlicedData$new(t(modaa[,-1]))

meJointaa = Matrix_eQTL_main(snps=snpSliceaa, gene = exprsaa, 
                               cvrt = covsaa, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all ct 
#-------------------------------
ctIndex=which(meths_snp_pca$Dx=="Control")


snpSlicect = SlicedData$new(t(snp3[ctIndex,]))
snpSlicect$ResliceCombined(sliceSize = 5000)

exprsct = SlicedData$new(t(methMat3[ctIndex,])) 
exprsct$ResliceCombined(sliceSize = 5000)

modct = mod[ctIndex,]

covsct = SlicedData$new(t(modct[,-1]))

meJointct = Matrix_eQTL_main(snps=snpSlicect, gene = exprsct, 
                               cvrt = covsct, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all sz 
#-------------------------------
szIndex=which(meths_snp_pca$Dx=="Schizo")


snpSlicesz = SlicedData$new(t(snp3[szIndex,]))
snpSlicesz$ResliceCombined(sliceSize = 5000)

exprssz = SlicedData$new(t(methMat3[szIndex,])) 
exprssz$ResliceCombined(sliceSize = 5000)

modsz = mod[szIndex,]

covssz = SlicedData$new(t(modsz[,-1]))

meJointsz = Matrix_eQTL_main(snps=snpSlicesz, gene = exprssz, 
                               cvrt = covssz, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all cauc.ct 
#-------------------------------
cauc.ctIndex=which((meths_snp_pca$Race=="CAUC") & (meths_snp_pca$Dx=="Control"))


snpSlicecauc.ct = SlicedData$new(t(snp3[cauc.ctIndex,]))
snpSlicecauc.ct$ResliceCombined(sliceSize = 5000)

exprscauc.ct = SlicedData$new(t(methMat3[cauc.ctIndex,])) 
exprscauc.ct$ResliceCombined(sliceSize = 5000)

modcauc.ct = mod[cauc.ctIndex,]

covscauc.ct = SlicedData$new(t(modcauc.ct[,-1]))

meJointcauc.ct = Matrix_eQTL_main(snps=snpSlicecauc.ct, gene = exprscauc.ct, 
                               cvrt = covscauc.ct, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all cauc.sz 
#-------------------------------
cauc.szIndex=which((meths_snp_pca$Race=="CAUC") & (meths_snp_pca$Dx=="Schizo"))


snpSlicecauc.sz = SlicedData$new(t(snp3[cauc.szIndex,]))
snpSlicecauc.sz$ResliceCombined(sliceSize = 5000)

exprscauc.sz = SlicedData$new(t(methMat3[cauc.szIndex,])) 
exprscauc.sz$ResliceCombined(sliceSize = 5000)

modcauc.sz = mod[cauc.szIndex,]

covscauc.sz = SlicedData$new(t(modcauc.sz[,-1]))

meJointcauc.sz = Matrix_eQTL_main(snps=snpSlicecauc.sz, gene = exprscauc.sz, 
                               cvrt = covscauc.sz, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all aa.ct 
#-------------------------------
aa.ctIndex=which((meths_snp_pca$Race=="AA") & (meths_snp_pca$Dx=="Control"))


snpSliceaa.ct = SlicedData$new(t(snp3[aa.ctIndex,]))
snpSliceaa.ct$ResliceCombined(sliceSize = 5000)

exprsaa.ct = SlicedData$new(t(methMat3[aa.ctIndex,])) 
exprsaa.ct$ResliceCombined(sliceSize = 5000)

modaa.ct = mod[aa.ctIndex,]

covsaa.ct = SlicedData$new(t(modaa.ct[,-1]))

meJointaa.ct = Matrix_eQTL_main(snps=snpSliceaa.ct, gene = exprsaa.ct, 
                               cvrt = covsaa.ct, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	

# all aa.sz
#-------------------------------
aa.szIndex=which((meths_snp_pca$Race=="AA") & (meths_snp_pca$Dx=="Schizo"))


snpSliceaa.sz = SlicedData$new(t(snp3[aa.szIndex,]))
snpSliceaa.sz$ResliceCombined(sliceSize = 5000)

exprsaa.sz = SlicedData$new(t(methMat3[aa.szIndex,])) 
exprsaa.sz$ResliceCombined(sliceSize = 5000)

modaa.sz = mod[aa.szIndex,]

covsaa.sz = SlicedData$new(t(modaa.sz[,-1]))

meJointaa.sz = Matrix_eQTL_main(snps=snpSliceaa.sz, gene = exprsaa.sz, 
                               cvrt = covsaa.sz, output_file_name.cis =  ".ctxt" ,
                               pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                               snpspos = snpspos, genepos = methpos, 
                               useModel = modelLINEAR,	cisDist=5e6,
                               pvalue.hist = 100,min.pv.by.genesnp = T)	



#save
#-------------------
save()




