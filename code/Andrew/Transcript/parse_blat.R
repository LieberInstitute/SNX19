library(jaffelab)
library(GenomicRanges)
library(dplyr)

## read in text version
hits = read.delim("SNX19_BLAT-PSL-output-style-15-transcripts-021517-cur.txt",
	as.is=TRUE, header=FALSE, row.names=1)
	
## colnames from http://useast.ensembl.org/info/website/upload/psl.html
colnames(hits) = c("matches", "misMatches", "repMatches",
	"nCount", "qNumInsert", "qBaseInsert","tNumInsert",
	"tBaseInsert", "strand", "qName", "qSize",
	"qStart", "qEnd", "tName", "tSize", "tStart",
	"tEnd", "blockCount", "blockSizes", "qStarts",
	"tStarts")

# get sizes and starts
sizeList = lapply(strsplit(hits$blockSizes, ","),as.numeric)
startList = lapply(strsplit(hits$tStarts, ","), as.numeric)

# create ends
endList = mapply(function(x,d) x+d, startList,sizeList)

## make GRanges
grList = mapply(function(s,e) GRanges("chr11",
	IRanges(start=s,end = e), strand="-"), startList,endList)
names(grList) = hits$qName
grList = GRangesList(grList)
save(grList, file="snx19_transcripts.rda")

## test as bed file into browser
write.table(as.data.frame(grList[[1]])[,1:3],sep=" ",
	row.names=FALSE, col.names=FALSE,quote=FALSE)