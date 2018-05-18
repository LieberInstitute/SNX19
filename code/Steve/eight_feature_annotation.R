#Getting annotattion information SNX19 8features of interest 
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.rda')

### SCRIPT UPDATEDS 03 22 18 to add new transcript feautres
## load features of interest
#foi = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_100517.csv',stringsAsFactors=FALSE)
foi = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_03222018.csv', stringsAsFactors=FALSE)

foi$Tag[foi$Tag=="junc8c.8d\xa0"] = "junc8c.8d"#Fixing Tag which is misread
foi = cbind(foi, cisMerge[match(foi$Feature, cisMerge$Feature ),c('exprsType','exprsClass','ucsc_feature','NumTx','WhichTx','exprsChr','exprsStart','exprsEnd') ] )

##
foi$MergeName <- foi$Feature
foi$MergeName[foi$exprsType=="Exon"] <- foi[foi$exprsType=="Exon",'Tag']

### Lift SNX19 junctions
library('GenomicRanges')
library(AnnotationHub)
library(rtracklayer)

ahub <- AnnotationHub()
ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
query(ahub.chain, c("hg38", "hg19"))
chain <- ahub.chain[ahub.chain$title == "hg19ToHg38.over.chain.gz"]
chain <- chain[[1]]

## Lift end position of hg19 junctions to hg38
hg19_jxn_start = makeGRangesFromDataFrame (data.frame(start=foi$exprsStart[foi$exprsType=="Junction"], end=foi$exprsStart[foi$exprsType=="Junction"], chr='chr11'))
names(hg19_jxn_start) =foi$MergeName[foi$exprsType=="Junction"]
## Lift start position of hg19 junctions to hg38
hg19_jxn_stop = makeGRangesFromDataFrame (data.frame(start=foi$exprsEnd[foi$exprsType=="Junction"], end=foi$exprsEnd[foi$exprsType=="Junction"], chr='chr11'))
names(hg19_jxn_stop) =foi$MergeName[foi$exprsType=="Junction"]

## hg19 junctions lifted to hg38
hg19_lifted_to_hg38_start = liftOver(hg19_jxn_start, chain)
hg19_lifted_to_hg38_stop = liftOver(hg19_jxn_stop, chain)

#Drop any junctions in which start or end position multimaps
idx = lengths(hg19_lifted_to_hg38_start)==1 & lengths(hg19_lifted_to_hg38_stop)==1
liftedStart = unlist(hg19_lifted_to_hg38_start[idx] )
liftedStop = unlist(hg19_lifted_to_hg38_stop[idx] )
##
df1 = data.frame(start=start(liftedStart), 
				 end = start(liftedStop),
				 strand = strand(liftedStart),
				 seqnames=seqnames(liftedStart),
				 Name = names(hg19_lifted_to_hg38_start[idx])  )				 
df1$hg19Liftedhg38_chrPos = paste0(df1$seqnames,":",df1$start,"-",df1$end, "(", "*",")")
colnames(df1)[-6] = paste0(colnames(df1),"_hg19")[-6]
##
foi$MergeName[foi$exprsType=="Junction"] = df1[match(df1$Name_hg19, foi$MergeName[foi$exprsType=="Junction"]),'hg19Liftedhg38_chrPos']
##
#save(foi,file='/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild.rda')
rownames(foi)=foi$Tag
foi = foi[c('junc8.10',
'junc8.8a',
'junc8c.9',
'ENSE00001757032.2',
'ENSE00002148491.1',
'junc1.2_long',
'junc2.3_short',
'junc4.6',
'junc8.9',
'SNX19 Gene'),]
rownames(foi) = NULL
save(foi,file='/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')

write.csv(foi,file='/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_10Features_of_Interest_crossBuild_03_22_18.csv',row.names=F)
foi_bed =  foi[,c('exprsChr','exprsStart','exprsEnd','Tag','Feature')]
colnames(foi_bed) <- c("chrom",'chromStart','chromEnd','name','strand')
foi_bed$score = 0
foi_bed = foi_bed[,c("chrom",'chromStart','chromEnd','name','score','strand')]
foi_bed$strand <- gsub(".*\\(","",foi_bed$strand)
foi_bed$strand <-gsub("\\)","",foi_bed$strand)
foi_bed$strand[foi_bed$strand!="+"&foi_bed$strand!="-"] <- "."
write.table(foi_bed,file='/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_10Features_of_Interest_crossBuild_03_22_18.bed',row.names=F,sep='\t',quote=F,col.names=F)
