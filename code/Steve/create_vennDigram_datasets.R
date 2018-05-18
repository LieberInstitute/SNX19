library(dplyr)
library(tidyr)
## Load features of interest
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')

## venn diagram for Gwa positive SNPs 
#load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.rda')
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/tidyStats_FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
######## Gwas sig best hits
bestHitsGwasSig = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature & !is.na(pgc52_P) & pgc52_P<5e-8 )
bestHitsGwasSig$Tag <- foi$Tag[match(bestHitsGwasSig$Feature,foi$Feature)]
bestHitsGwasSig = bestHitsGwasSig %>% group_by(Region,Tag) %>% as.data.frame()
bestHitsGwasSig = bestHitsGwasSig[order(match(bestHitsGwasSig$Region, unique(bestHitsGwasSig$Region) ),match(bestHitsGwasSig$Tag, foi$Tag)),]

openxlsx::write.xlsx(list(`p_1e5`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-5, ],
	 `p_1e6`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-6, ],
	 `p_1e7`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-7, ],
	 `p_1e8`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-8, ],
	 `p_1e9`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-9, ],
	 `p_1e10`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-10, ],
	 `p_1e11`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-11, ],
	 `p_1e12`=bestHitsGwasSig[bestHitsGwasSig$pvalue<1e-12, ]), 
	 file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/gwas_SNP_sig_eqtl_5datasets_10features.xlsx',row.names=FALSE)

##
colOrder = paste0( rep(c("DLPFC_PolyA ","DLPFC_RiboZero ", "CMC_DLPFC ", "Hippo ","Caudate " ),each=10 ), rep(foi$Tag,5)  )
colOrder = c('id', colOrder)  

bestHitsGwasSig = unite(bestHitsGwasSig,RegionFeature, Region, Tag, sep=' ') %>% 
				  mutate(id=1:n())
				  
p_1e5 =bestHitsGwasSig %>% filter(pvalue<1e-5) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA,drop=FALSE)
p_1e5[,colOrder[!colOrder%in%colnames(p_1e5)]  ] <- NA
  
p_1e6 =bestHitsGwasSig %>% filter(pvalue<1e-6) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)
p_1e6[,colOrder[!colOrder%in%colnames(p_1e6)]  ] <- NA
  
  
p_1e7 =bestHitsGwasSig %>% filter(pvalue<1e-7) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)
p_1e7[,colOrder[!colOrder%in%colnames(p_1e7)]  ] <- NA
  
p_1e8 =bestHitsGwasSig %>% filter(pvalue<1e-8) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)
p_1e8[,colOrder[!colOrder%in%colnames(p_1e8)]  ] <- NA  
  
p_1e9 =bestHitsGwasSig %>% filter(pvalue<1e-9) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)
p_1e9[,colOrder[!colOrder%in%colnames(p_1e9)]  ] <- NA  
  
p_1e10 =bestHitsGwasSig %>% filter(pvalue<1e-10) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)  
p_1e10[,colOrder[!colOrder%in%colnames(p_1e10)]  ] <- NA  
  
p_1e11 =bestHitsGwasSig %>% filter(pvalue<1e-11) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)
p_1e11[,colOrder[!colOrder%in%colnames(p_1e11)]  ] <- NA  
    
p_1e12 =bestHitsGwasSig %>% filter(pvalue<1e-12) %>% select( c('SNP','RegionFeature')  )%>% filter(!is.na(SNP)) %>% # remove NA values
  unique %>%mutate(id=SNP) %>%  # remove duplicated rows
  spread(RegionFeature,SNP,fill=NA)   
p_1e12[,colOrder[!colOrder%in%colnames(p_1e12)]  ] <- NA  
  
########

openxlsx::write.xlsx(list(`p_1e5`=p_1e5[,colOrder],
	 `p_1e6`=p_1e6[,colOrder],
	 `p_1e7`=p_1e7[,colOrder],
	 `p_1e8`=p_1e8[,colOrder],
	 `p_1e9`=p_1e9[,colOrder],
	 `p_1e10`=p_1e10[,colOrder],
	 `p_1e11`=p_1e11[,colOrder],
	 `p_1e12`=p_1e12[,colOrder]), 
	 file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/gwas_SNP_sig_eqtl_5datasets_10features_spread.xlsx',row.names=FALSE)

##  