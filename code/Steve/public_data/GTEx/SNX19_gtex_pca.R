#qsub -l bluejay -l mf=30G,h_vmem=50G,h_stack=256M -M stephen.semick@gmail.com -cwd -b y R CMD BATCH SNX19_gtex_pca.R
#load objects
load("/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/SNX19_allFeature_EqtlGtex.rda")

#Throw out failed datasets
SNX19_allFeature_EqtlGtex <- SNX19_allFeature_EqtlGtex[sapply(SNX19_allFeature_EqtlGtex,is.data.frame)]

#define matching reference
ref_UniqueID=SNX19_allFeature_EqtlGtex[[1]]$UniqueID

#
SNX19_allFeature_EqtlGtex_pca_prep <- lapply (names(SNX19_allFeature_EqtlGtex), function(x) { 

#Appending name of study to column names
y = SNX19_allFeature_EqtlGtex[[x]] 
colnames(y)[colnames(y)%in%c("All_statistic", "All_pvalue", "All_FDR", "All_beta")] <- paste0(x,"_",colnames(y)[colnames(y)%in%c("All_statistic", "All_pvalue", "All_FDR", "All_beta")] )

#Lining up all of eQTL Results 
y = y[match(ref_UniqueID,y$UniqueID), ]
return(y)
} )
names(SNX19_allFeature_EqtlGtex_pca_prep) <-  names(SNX19_allFeature_EqtlGtex)


###
reoriented_eqtl_res <- lapply(SNX19_allFeature_EqtlGtex_pca_prep[2:length(SNX19_allFeature_EqtlGtex_pca_prep)], function(x) {
x[,!colnames(x)%in%c("UniqueID", "snps", "feature", "All_beta")]
})

#Creating a matrix for pca
SNX19_allFeature_EqtlGtex_beta_pca_prep <- sapply(names(SNX19_allFeature_EqtlGtex_pca_prep), function(x) { 
y = SNX19_allFeature_EqtlGtex_pca_prep[[x]]
y = y[,grepl("All_beta",colnames(y))] 
return(y) } )

#principle component analysis and saving
gtex_snx19_pca <- prcomp(t(SNX19_allFeature_EqtlGtex_beta_pca_prep))
save(gtex_snx19_pca, file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/gtex_eqtl_beta_pca.rda')