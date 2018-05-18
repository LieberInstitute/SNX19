## Adding additional snx19_info to big table 03/19/2018 ##

#load template
template = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/fixes/template-Table1_v2_FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417_minusLogP_cur4.xlsx')

# load the data
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild.rda')

######## Adding tag
template$Tag <- NA
template$Tag = foi[match(template$Feature, foi$Feature),'Tag']

########### Adding directionality #############

#### directionality
dataSets =c("DlpfcPolya","DlpfcRibo","Caudate", "Hippo", "CMC")
eqtlTraits =c("All",
			  "Control_AA",
			  "Control_CAUC",
			  "Schizo_AA",
			  "Schizo_CAUC",
			  "CAUC",
			  "AA")
length_eqtl_traits=length(eqtlTraits)
length_gwas_traits = length(dataSets)	
## set up matrics
snpMatchMat = eqtlSignMat = riskSignMat = matrix(0, 
	ncol = length_gwas_traits*length_eqtl_traits+1, nrow = nrow(cisMerge))
colnames(snpMatchMat) = c(paste0(rep(dataSets,
	each=length_eqtl_traits), "_", rep(eqtlTraits, times=length_gwas_traits)), "Caudate_Bipolar_CAUC")
colnames(snpMatchMat) <- c(colnames(snpMatchMat)[1:21], "Caudate_Bipolar_CAUC",colnames(snpMatchMat)[22:35])
colnames(eqtlSignMat) = colnames(riskSignMat) = colnames(snpMatchMat)

# does snp match?
snpMatchMat[cisMerge$pgc52_A2 == cisMerge$snpCounted,] = 1
snpMatchMat[grepl("I", cisMerge$pgc52_A2) & cisMerge$snpType == "Insertion",] = 1
snpMatchMat[cisMerge$pgc52_A2 == "D" & cisMerge$snpType == "Deletion",] = 1
snpMatchMat[cisMerge$pgc52_A1 == cisMerge$snpCounted,] = -1
snpMatchMat[grepl("I", cisMerge$pgc52_A2)& cisMerge$snpType == "Deletion",] = -1
snpMatchMat[cisMerge$pgc52_A2 == "D" & cisMerge$snpType == "Insertion",] = -1

# risk allele
riskSignMat[cisMerge$pgc52_A2 == cisMerge$pgc52_riskAllele, 1:36] = 1
riskSignMat[cisMerge$pgc52_A1 == cisMerge$pgc52_riskAllele, 1:36] = -1

# risk allele

#DLPFC PolyA
eqtlSignMat[, 1] =sign(cisMerge$DLPFC_PolyA_Linear_All_statistic)
eqtlSignMat[, 2] =sign(cisMerge$DLPFC_PolyA_Linear_Control_AA_statistic)
eqtlSignMat[, 3] =sign(cisMerge$DLPFC_PolyA_Linear_Control_CAUC_statistic)
eqtlSignMat[, 4] =sign(cisMerge$DLPFC_PolyA_Linear_Schizo_AA_statistic)
eqtlSignMat[, 5] =sign(cisMerge$DLPFC_PolyA_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 6] =sign(cisMerge$DLPFC_PolyA_Linear_CAUC_statistic)
eqtlSignMat[, 7] =sign(cisMerge$DLPFC_PolyA_Linear_AA_statistic)

#DLPFC RiboZero
eqtlSignMat[, 8] =sign(cisMerge$DLPFC_RiboZero_Linear_All_statistic)
eqtlSignMat[, 9] =sign(cisMerge$DLPFC_RiboZero_Linear_Control_AA_statistic)
eqtlSignMat[, 10] =sign(cisMerge$DLPFC_RiboZero_Linear_Control_CAUC_statistic)
eqtlSignMat[, 11] =sign(cisMerge$DLPFC_RiboZero_Linear_Schizo_AA_statistic)
eqtlSignMat[, 12] =sign(cisMerge$DLPFC_RiboZero_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 13] =sign(cisMerge$DLPFC_RiboZero_Linear_CAUC_statistic)
eqtlSignMat[, 14] =sign(cisMerge$DLPFC_RiboZero_Linear_AA_statistic)

# Caudate
eqtlSignMat[, 15] =sign(cisMerge$Caudate_Linear_All_statistic)
eqtlSignMat[, 16] =sign(cisMerge$Caudate_Linear_Control_AA_statistic)
eqtlSignMat[, 17] =sign(cisMerge$Caudate_Linear_Control_CAUC_statistic)
eqtlSignMat[, 18] =sign(cisMerge$Caudate_Linear_Schizo_AA_statistic)
eqtlSignMat[, 19] =sign(cisMerge$Caudate_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 20] =sign(cisMerge$Caudate_Linear_CAUC_statistic)
eqtlSignMat[, 21] =sign(cisMerge$Caudate_Linear_AA_statistic)
eqtlSignMat[, 22] =sign(cisMerge$Caudate_Linear_Bipolar_CAUC_statistic)

# Hippo
eqtlSignMat[, 23] =sign(cisMerge$Hippo_Linear_All_statistic)
eqtlSignMat[, 24] =sign(cisMerge$Hippo_Linear_Control_AA_statistic)
eqtlSignMat[, 25] =sign(cisMerge$Hippo_Linear_Control_CAUC_statistic)
eqtlSignMat[, 26] =sign(cisMerge$Hippo_Linear_Schizo_AA_statistic)
eqtlSignMat[, 27] =sign(cisMerge$Hippo_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 28] =sign(cisMerge$Hippo_Linear_CAUC_statistic)
eqtlSignMat[, 29] =sign(cisMerge$Hippo_Linear_AA_statistic)

# CMC DLPFC
eqtlSignMat[, 30] =sign(cisMerge$CMC_DLPFC_Linear_All_statistic)
eqtlSignMat[, 31] =sign(cisMerge$CMC_DLPFC_Linear_Control_AA_statistic)
eqtlSignMat[, 32] =sign(cisMerge$CMC_DLPFC_Linear_Control_CAUC_statistic)
eqtlSignMat[, 33] =sign(cisMerge$CMC_DLPFC_Linear_Schizo_AA_statistic)
eqtlSignMat[, 34] =sign(cisMerge$CMC_DLPFC_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 35] =sign(cisMerge$CMC_DLPFC_Linear_CAUC_statistic)
eqtlSignMat[, 36] =sign(cisMerge$CMC_DLPFC_Linear_AA_statistic)

# combine	
directionalityMat = snpMatchMat*riskSignMat*eqtlSignMat
directionalityMat[directionalityMat==1] = "U"
directionalityMat[directionalityMat==-1] = "D"
directionalityMat[is.na(directionalityMat)] = "M"

# make strings
dirString = apply(directionalityMat,1, function(x) {
	paste(paste(x[1:7], collapse=""), 
		paste(x[8:14], collapse=""),  
		paste(x[15:22], collapse=""),
		paste(x[23:29], collapse=""),
		paste(x[30:36], collapse=""),
		sep="_")
})

dirString = paste(ss(dirString,"_" ,1),  ss(dirString,"_" ,2), ss(dirString,"_" ,5), ss(dirString,"_" ,4), ss(dirString,"_" ,3), sep="_" )
dirString = dirString[match(cisMerge$UniqueID,template$UniqueID)]
template$dirString=dirString
write.csv(template, '/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/fixes/ssAddedCols_template-Table1_v2_FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417_minusLogP_cur4.csv',row.names=F)