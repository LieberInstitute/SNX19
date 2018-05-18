
load("/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/SNX19_allFeature_EqtlGtex.rda")
SNX19_allFeature_EqtlGtex <- SNX19_allFeature_EqtlGtex[sapply(SNX19_allFeature_EqtlGtex,is.data.frame)]

snp_of_interest <- c('rs10791097',
'rs34529622',
'rs3794133',
'rs4337054',
'rs73028899',
'rs78968418',
'rs10791114')

lapply (SNX19_allFeature_EqtlGtex, function(x) { 
SNX19_allFeature_EqtlGtex[[x]] 
} )

ref_UniqueID=SNX19_allFeature_EqtlGtex[[1]]$UniqueID

SNX19_allFeature_EqtlGtex_pca_prep <- lapply (names(SNX19_allFeature_EqtlGtex), function(x) { 

#Appending name of study to column names
y = SNX19_allFeature_EqtlGtex[[x]] 
colnames(y)[colnames(y)%in%c("All_statistic", "All_pvalue", "All_FDR", "All_beta")] <- paste0(x,"_",colnames(y)[colnames(y)%in%c("All_statistic", "All_pvalue", "All_FDR", "All_beta")] )

#Lining up all of eQTL Results 
y = y[match(ref_UniqueID,y$UniqueID), ]
return(y)
} )
names(SNX19_allFeature_EqtlGtex_pca_prep) <-  names(SNX19_allFeature_EqtlGtex)

#Creating a matrix
SNX19_allFeature_EqtlGtex_beta_pca_prep <- sapply(names(SNX19_allFeature_EqtlGtex_pca_prep), function(x) { 
y = SNX19_allFeature_EqtlGtex_pca_prep[[x]]
y = y[,grepl("All_beta",colnames(y))] 
return(y) } )

gtex_snx19_pca <- prcomp(t(SNX19_allFeature_EqtlGtex_beta_pca_prep))