#compare original coexpression to new coexpression

#load original coexpression
old_coexp = read.csv( gzfile('/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/Multivariate_Coexpression_SNX19_FullStats.csv.gz') )
#load new coexpression
load('/dcl01/lieber/ajaffe/Steve/SNX19/coexpression/DLPFC_polyA/rdas/Multivariate_Coexpression_SNX19_PolyA_FullStats.rda') #outstats

head(outStats)

summary(abs(outStats[,"Control_CAUC.junc8.10.tstat"]- old_coexp[,"CAUC_CONT.D9.tstat"]))
cor(outStats[,"Control_CAUC.junc8.10.tstat"],old_coexp[,"CAUC_CONT.D9.tstat"])

max(abs(-log10(outStats[,"Control_CAUC.junc8.10.pval"])- -log10(old_coexp[,"CAUC_CONT.D9.pval"]) ) )
plot(-log10(outStats[,"Control_CAUC.junc8.10.pval"]), -log10(old_coexp[,"CAUC_CONT.D9.pval"])  )


outStats$"featureID" == old_coexp$"featureID"

