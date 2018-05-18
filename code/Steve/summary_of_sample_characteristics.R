#Table One
library(tableone)
setwd('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_polyA')
#load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 
library(LIBDpheno)
library(ggplot2)
library(reshape2)
#####################
### DLPFC PolyA ###
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")
keepIndex= which(pd$age >= 13) #Indices for subjects above the age of 13 
pd_DLPFC_polyA= pd[keepIndex,] #Keeping phenotypic data information only for the subjects above the age of 13
pd_DLPFC_polyA$Region = "DLPFC_PolyA"

### DLPFC RiboZero ###
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_RiboZero_n485_updateMap.rda")
keepIndex= which(pd$age >= 13 & pd$Race %in% c("AA","CAUC") ) #Indices for subjects above the age of 13 
pd_DLPFC_RiboZero= pd[keepIndex,] #Keeping phenotypic data information only for the subjects above the age of 13
pd_DLPFC_RiboZero$Region = "DLPFC_RiboZero"

### Caudate RiboZero ###
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/Caudate_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_Caudate_n468_updateMap.rda")
keepIndex= which(pd$age >= 13 & pd$Race %in% c("AA","CAUC") ) #Indices for subjects above the age of 13 
pd_Caudate= pd[keepIndex,] #Keeping phenotypic data information only for the subjects above the age of 13
pd_Caudate$Region = "Caudate"

### Hippo RiboZero ###
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/HIPPO_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_HIPPO_RiboZero_n452_updateMap.rda")
pd = plyr::rename(pd, c('Age'='age') )
keepIndex= which(pd$age >= 13 & pd$Race %in% c("AA","CAUC") ) #Indices for subjects above the age of 13 
pd_Hippo = pd[keepIndex,] #Keeping phenotypic data information only for the subjects above the age of 13
pd_Hippo$Region = "Hippo"
### CMC DLPFC ###
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/CMC/rawAndRpkmCounts_plusGenotype_SNX19_CMC_n547.rda")
pd$Age_of_Death[pd$Age_of_Death=="90+"] <- NA
pd <- plyr::rename(pd, c("Age_of_Death"="age"))
pd <- plyr::rename(pd, c("Gender"="Sex"))
pd$Sex = plyr::mapvalues(pd$Sex, c("Male", "Female"), c("M","F") )

pd <- plyr::rename(pd, c("DLPFC_RNA_isolation_RIN"="RIN"))
pd$age <- as.numeric(pd$age)
pd$Race = plyr::mapvalues(pd$Ethnicity, c("African-American", "Caucasian"), c("AA","CAUC") )
pd$Dx = plyr::mapvalues(pd$Dx, c("SCZ"), c("Schizo") )

keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
pd_CMC = pd[keepIndex,]
pd_CMC$Region="CMC DLPFC"

### Combined pd object ###
cols_interest = c('age','Race','Sex',"Region","Dx", "RIN" )
pd = rbind(pd_DLPFC_polyA[ ,cols_interest],
		   pd_DLPFC_RiboZero[ ,cols_interest],
		   pd_Caudate[ ,cols_interest], 
		   pd_Hippo[ ,cols_interest],
		   pd_CMC[ ,cols_interest])

table1 <- CreateTableOne(vars = c("age","Race","Sex", "RIN"), strata = c("Dx","Region"), data = pd)
table1 <- print(table1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, pDigits=100)
## Save to a CSV file
write.csv(table1, file = "/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/TableI.csv", row.names=TRUE)