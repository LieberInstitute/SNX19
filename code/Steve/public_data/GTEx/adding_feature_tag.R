## Fixing random GTEX SNX19 things
gtex = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/for-steve-72799-GTEX.xlsx')

## foi 
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')

gtex$FeatureName =NA
gtex$FeatureName <- foi[match(gtex$feature, foi$Feature),'Tag']

write.csv(gtex,file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/gtex_featurename_added.csv',row.names=F)