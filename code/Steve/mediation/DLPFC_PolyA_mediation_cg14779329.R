############### DLPFC POLYA ################# 
## Load DLPFC PolyA Data
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")
pd_rna = pd

##
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/SNX19_DNAm_data_DLPFC.rda")
pd_meth = methPheno 

## subset to shared people and line everything up
shared_people=intersect(pd_meth$BrNum, pd_rna$BrNum)

pd_rna = pd_rna[match(shared_people, pd_rna$BrNum), ]
pd_meth = pd_meth[match(shared_people, pd_meth$BrNum), ]

snpAll = snpAll[,shared_people]
methMat = methMat[,shared_people]
yExprs_DlpfcPolyA = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1)[,pd_rna$RNum]

## subset phenotype
keepIndex = which(pd_rna$Race %in% c("CAUC", "AA") &  pd_rna$age > 13)

pd_rna = pd_rna[keepIndex, ]
pd_meth = pd_meth[keepIndex, ]

snpAll = snpAll[, keepIndex]
methMat = methMat[, keepIndex]
yExprs_DlpfcPolyA = yExprs_DlpfcPolyA[,keepIndex ]

## Load in features of interest
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')

## load statistics
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/tidyStats_FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
library(dplyr)
library(tidyr)

bestHits = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature & !is.na(pgc52_P) &pgc52_P<5e-8 ) %>% 
						 group_by(Region,Feature) %>%  
						 slice(which.min(pvalue)) %>% 
						 as.data.frame()				 
bestHits_DlpfcPolyA = bestHits[bestHits$Region=='DLPFC_PolyA',]

## merge objects
dat = cbind(pd_rna, t(yExprs_DlpfcPolyA[foi$Feature,]) )
colnames(dat)[match(foi$Feature, colnames(dat))] <- foi[match( foi$Feature, colnames(dat)[match(foi$Feature, colnames(dat))], ),'Tag']
## load PCs
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_polyA/rdas/all3Feature_PCs_DLPFC_polyA_age13_matchedUp.rda',verbose=T)
all3PCs = all3PCs[dat$RNum,]
load("/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/rdas/meth_PCs_DLPFC_polyA_age13_matchedUp.rda",verbose=T)
pca_meth = pca_meth[dat$BrNum,]
##
mod = model.matrix(~dat$snpPC1 + dat$snpPC2 + dat$snpPC3 + dat$snpPC4 + dat$snpPC5 )

getMediationStats = function(row_i) {

## subset to features of interest
snp_i = bestHits_DlpfcPolyA[row_i, 'SNP' ]
cpg_i = methMat['cg14779329', ]
feature_i = bestHits_DlpfcPolyA[row_i, 'Feature']

keepSNP = !is.na(snpAll[snp_i,])

##### Calculate a: effect of SNP on methylation
mod_a = cbind(mod, pca_meth, t(snpAll[snp_i,])  )
dat_a = data.frame(cg14779329 = cpg_i[keepSNP], mod_a[keepSNP,-1] )
colnames(dat_a)[ncol(dat_a)] <- c(snp_i)

fastResA=summary(lm(cg14779329~ ., data=dat_a) )
fastResA=data.frame(coef=rownames(fastResA$coefficients),fastResA$coefficients, stringsAsFactors=FALSE)
colnames(fastResA) <- c('coef','Estimate','StdErr','t.value','P')
fastResA$coef = gsub("`","",fastResA$coef)
#fastResA$SNP = snp_i
#fastResA$CpG = cpg_i
fastResA = fastResA[fastResA$coef==snp_i ,]

##### Calculate b: effect of methylation on expression
mod_b = cbind(mod, all3PCs, pca_meth, t(snpAll[snp_i,])  )
dat_b = data.frame(cg14779329 = cpg_i[keepSNP], mod_b[keepSNP,-1], t(yExprs_DlpfcPolyA[feature_i,keepSNP]) )
colnames(dat_b)[(ncol(dat_b)-1):ncol(dat_b)] <- c(snp_i, 'Feature')

fastResB = summary( lm(Feature~.,dat_b) )
fastResB=data.frame(coef=rownames(fastResB$coefficients),fastResB$coefficients, stringsAsFactors=FALSE)
colnames(fastResB) <- c('coef','Estimate','StdErr','t.value','P')
fastResB$coef = gsub("`","",fastResB$coef)
fastResCPrime = fastResB[fastResB$coef==snp_i,]
fastResB = fastResB[fastResB$coef=='cg14779329',]

##### Calculate c: full effect
mod_c = cbind(mod, all3PCs, t(snpAll[snp_i,])  )
dat_c = data.frame(mod_c[keepSNP,-1], t(yExprs_DlpfcPolyA[feature_i,keepSNP]) )
colnames(dat_c)[(ncol(dat_c)-1):ncol(dat_c)] <- c(snp_i, 'Feature')

fastResC = summary( lm(Feature~.,dat_c) )
fastResC=data.frame(coef=rownames(fastResC$coefficients),fastResC$coefficients, stringsAsFactors=FALSE)
colnames(fastResC) <- c('coef','Estimate','StdErr','t.value','P')
fastResC$coef = gsub("`","",fastResC$coef)
fastResC = fastResC[fastResC$coef==snp_i,]

##### Calculate b naive: effect of methylation on expression
mod_b0 = cbind(mod, pca_meth )
dat_b0 = data.frame(cg14779329 = cpg_i[keepSNP], mod_b0[keepSNP,-1], t(yExprs_DlpfcPolyA[feature_i,keepSNP]) )
colnames(dat_b0)[ncol(dat_b0)] <- c( 'Feature')

fastResB0 = summary( lm(Feature~.,dat_b0) )
fastResB0=data.frame(coef=rownames(fastResB0$coefficients),fastResB0$coefficients, stringsAsFactors=FALSE)
colnames(fastResB0) <- c('coef','Estimate','StdErr','t.value','P')
fastResB0$coef = gsub("`","",fastResB0$coef)
fastResB0 = fastResB0[fastResB0$coef=='cg14779329',]

##
res = data.frame(SNP=snp_i,
		CpG='cg14779329',
		Feature=foi[match(feature_i,foi$Feature),'Tag'] ,
		
		a.Estimate_SNP_x_CpG = fastResA$Estimate,
		b.Estimate_CpG_x_Exprs = fastResB0$Estimate,
		c.Estimate_SNP_x_Exprs = fastResC$Estimate, 		
		bPrime.Estimate_CpG_x_Exprs = fastResB$Estimate,
		cPrime.Estimate_SNP_x_Exprs = fastResCPrime$Estimate, 

		a.SE_SNP_x_CpG = fastResA$StdErr,
		b.SE_CpG_x_Exprs = fastResB0$StdErr,
		c.SE_SNP_x_Exprs = fastResC$StdErr,		
		bPrime.SE_CpG_x_Exprs = fastResB$StdErr,
		cPrime.SE_SNP_x_Exprs = fastResCPrime$StdErr,
		
		a.Pval_SNP_x_CpG = fastResA$P,
		b.Pval_CpG_x_Exprs = fastResB0$P,		
		c.Pval_SNP_x_Exprs = fastResC$P,
		bPrime.Pval_CpG_x_Exprs = fastResB$P,		
		cPrime.Pval_SNP_x_Exprs = fastResCPrime$P
		)
		
### Baron and Kenny analysis		
res$BK_Step1_to_3 = ifelse(res$a.Pval_SNP_x_CpG<0.01 & res$b.Pval_CpG_x_Exprs <0.01 & res$c.Pval_SNP_x_Exprs <0.01, "Mediation Possible", "Mediation Unlikely")
res$BK_Step4 = ifelse(res$BK_Step1_to_3=="Mediation Unlikely", "Steps 1-3 Failed", ifelse(res$bPrime.Pval_CpG_x_Exprs>0.01, "No Mediation", ifelse(res$cPrime.Pval_SNP_x_Exprs<0.01, "Partial Mediation", "Full Mediation")  ) )

### Sobel test for mediation
a=res$`a.Estimate_SNP_x_CpG`
sa=res$`a.SE_SNP_x_CpG`
b=res$`bPrime.Estimate_CpG_x_Exprs`
sb=res$`bPrime.SE_CpG_x_Exprs`

## Sobel Test for mediation
res$Sobel_Mediation_Z = a*b/sqrt(b^2*sa^2+a^2*sb^2)
res$Sobel_Mediation_P = pnorm(-abs(res$Sobel_Mediation_Z))*2

## Aroian Test for mediation
res$Aroian_Mediation_Z = a*b/sqrt(b^2*sa^2 + a^2*sb^2 + sa^2*sb^2)
res$Aroian_Mediation_P = pnorm(-abs(res$Aroian_Mediation_Z))*2
##
return(res)
 }

mediationStats = data.table::rbindlist(lapply(1:nrow(bestHits_DlpfcPolyA), getMediationStats))

write.csv(as.data.frame(mediationStats), file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/mediation_stats_10features_bestSNP_dlpfc_polya_cg14779329_mediator.csv',row.names=F )
######################## Now with rs10791097
getMediationStats = function(row_i) {

## subset to features of interest
snp_i = snpMapAll[which(snpMapAll$name=='rs10791097'),'SNP']
cpg_i = methMat['cg14779329', ]
feature_i = bestHits_DlpfcPolyA[row_i, 'Feature']

keepSNP = !is.na(snpAll[snp_i,])

##### Calculate a: effect of SNP on methylation
mod_a = cbind(mod, pca_meth, t(snpAll[snp_i,])  )
dat_a = data.frame(cg14779329 = cpg_i[keepSNP], mod_a[keepSNP,-1] )
colnames(dat_a)[ncol(dat_a)] <- c(snp_i)

fastResA=summary(lm(cg14779329~ ., data=dat_a) )
fastResA=data.frame(coef=rownames(fastResA$coefficients),fastResA$coefficients, stringsAsFactors=FALSE)
colnames(fastResA) <- c('coef','Estimate','StdErr','t.value','P')
fastResA$coef = gsub("`","",fastResA$coef)
#fastResA$SNP = snp_i
#fastResA$CpG = cpg_i
fastResA = fastResA[fastResA$coef==snp_i ,]

##### Calculate b: effect of methylation on expression
mod_b = cbind(mod, all3PCs, pca_meth, t(snpAll[snp_i,])  )
dat_b = data.frame(cg14779329 = cpg_i[keepSNP], mod_b[keepSNP,-1], t(yExprs_DlpfcPolyA[feature_i,keepSNP]) )
colnames(dat_b)[(ncol(dat_b)-1):ncol(dat_b)] <- c(snp_i, 'Feature')

fastResB = summary( lm(Feature~.,dat_b) )
fastResB=data.frame(coef=rownames(fastResB$coefficients),fastResB$coefficients, stringsAsFactors=FALSE)
colnames(fastResB) <- c('coef','Estimate','StdErr','t.value','P')
fastResB$coef = gsub("`","",fastResB$coef)
fastResCPrime = fastResB[fastResB$coef==snp_i,]
fastResB = fastResB[fastResB$coef=='cg14779329',]

##### Calculate c: full effect
mod_c = cbind(mod, all3PCs, t(snpAll[snp_i,])  )
dat_c = data.frame(mod_c[keepSNP,-1], t(yExprs_DlpfcPolyA[feature_i,keepSNP]) )
colnames(dat_c)[(ncol(dat_c)-1):ncol(dat_c)] <- c(snp_i, 'Feature')

fastResC = summary( lm(Feature~.,dat_c) )
fastResC=data.frame(coef=rownames(fastResC$coefficients),fastResC$coefficients, stringsAsFactors=FALSE)
colnames(fastResC) <- c('coef','Estimate','StdErr','t.value','P')
fastResC$coef = gsub("`","",fastResC$coef)
fastResC = fastResC[fastResC$coef==snp_i,]

##### Calculate b naive: effect of methylation on expression
mod_b0 = cbind(mod, pca_meth )
dat_b0 = data.frame(cg14779329 = cpg_i[keepSNP], mod_b0[keepSNP,-1], t(yExprs_DlpfcPolyA[feature_i,keepSNP]) )
colnames(dat_b0)[ncol(dat_b0)] <- c( 'Feature')

fastResB0 = summary( lm(Feature~.,dat_b0) )
fastResB0=data.frame(coef=rownames(fastResB0$coefficients),fastResB0$coefficients, stringsAsFactors=FALSE)
colnames(fastResB0) <- c('coef','Estimate','StdErr','t.value','P')
fastResB0$coef = gsub("`","",fastResB0$coef)
fastResB0 = fastResB0[fastResB0$coef=='cg14779329',]

##
res = data.frame(SNP=snp_i,
		CpG='cg14779329',
		Feature=foi[match(feature_i,foi$Feature),'Tag'] ,
		
		a.Estimate_SNP_x_CpG = fastResA$Estimate,
		b.Estimate_CpG_x_Exprs = fastResB0$Estimate,
		c.Estimate_SNP_x_Exprs = fastResC$Estimate, 		
		bPrime.Estimate_CpG_x_Exprs = fastResB$Estimate,
		cPrime.Estimate_SNP_x_Exprs = fastResCPrime$Estimate, 

		a.SE_SNP_x_CpG = fastResA$StdErr,
		b.SE_CpG_x_Exprs = fastResB0$StdErr,
		c.SE_SNP_x_Exprs = fastResC$StdErr,		
		bPrime.SE_CpG_x_Exprs = fastResB$StdErr,
		cPrime.SE_SNP_x_Exprs = fastResCPrime$StdErr,
		
		a.Pval_SNP_x_CpG = fastResA$P,
		b.Pval_CpG_x_Exprs = fastResB0$P,		
		c.Pval_SNP_x_Exprs = fastResC$P,
		bPrime.Pval_CpG_x_Exprs = fastResB$P,		
		cPrime.Pval_SNP_x_Exprs = fastResCPrime$P
		)
		
### Baron and Kenny analysis		
res$BK_Step1_to_3 = ifelse(res$a.Pval_SNP_x_CpG<0.01 & res$b.Pval_CpG_x_Exprs <0.01 & res$c.Pval_SNP_x_Exprs <0.01, "Mediation Possible", "Mediation Unlikely")
res$BK_Step4 = ifelse(res$BK_Step1_to_3=="Mediation Unlikely", "Steps 1-3 Failed", ifelse(res$bPrime.Pval_CpG_x_Exprs>0.01, "No Mediation", ifelse(res$cPrime.Pval_SNP_x_Exprs<0.01, "Partial Mediation", "Full Mediation")  ) )

### Sobel test for mediation
a=res$`a.Estimate_SNP_x_CpG`
sa=res$`a.SE_SNP_x_CpG`
b=res$`bPrime.Estimate_CpG_x_Exprs`
sb=res$`bPrime.SE_CpG_x_Exprs`

## Sobel Test for mediation
res$Sobel_Mediation_Z = a*b/sqrt(b^2*sa^2+a^2*sb^2)
res$Sobel_Mediation_P = pnorm(-abs(res$Sobel_Mediation_Z))*2

## Aroian Test for mediation
res$Aroian_Mediation_Z = a*b/sqrt(b^2*sa^2 + a^2*sb^2 + sa^2*sb^2)
res$Aroian_Mediation_P = pnorm(-abs(res$Aroian_Mediation_Z))*2
##
return(res)
 }

mediationStats = data.table::rbindlist(lapply(1:nrow(bestHits_DlpfcPolyA), getMediationStats))
write.csv(as.data.frame(mediationStats), file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/mediation_stats_10features_rs10791097_dlpfc_polya_cg14779329_mediator.csv',row.names=F )