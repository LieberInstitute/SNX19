## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(RColorBrewer)

## SNX19 Features over times
load('/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/SNX19/stemcell_SNX19_features.rda')

snxRpkm = log2(snxRpkm+1)

## 
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5)) )
pdf("/dcl01/lieber/ajaffe/Steve/SNX19/plots/StemCell_8features_singleGenomes_timecourse.pdf",h=8.5,w=11)
for (i in unique(pd$Donor) ) {
dat = cbind(as.data.frame(pd),t(snxRpkm) )
dat = dat[dat$Donor %in% i, ]
dat2=tidyr::gather(dat, Feature,Expression, (ncol(dat)-7):ncol(dat) )
dat2$Feature = foi[match(dat2$Feature,rownames(foi)),'Tag']

a <- ggplot(data=dat2, aes(x=DAY, y=Expression, col = Feature)) +
		geom_point(size=2.5) + geom_smooth(aes(group=Feature), method="loess", size=1.5, se=F) +
		theme(plot.title = element_text(hjust = 0.5),
			  legend.position="right") +
		scale_colour_brewer(palette="Dark2") +
		scale_fill_brewer(palette="Dark2") +
		labs( x = "Day", y="log2(Exprs+1)", title = paste0("Donor ", i ))
print(a)		
}
dev.off()

# Colors by Donor
lineIndexes = splitit(paste0(pd$LINE,":",pd$REP))

pal = c(colorRampPalette(brewer.pal(10,"Paired")[1:2])(5), 
		colorRampPalette(brewer.pal(10,"Paired")[3:4])(3), 
		brewer.pal(10,"Paired")[9:10], brewer.pal(10,"Paired")[5:8])
pd$lineCol = pal[as.numeric(factor(pd$LINE))]

#### Expression of NANOG gene across timecourse
snxRpkm

## Day as factor
pdf("/dcl01/lieber/ajaffe/Steve/SNX19/plots/StemCell_8features_allGenomes_timecourse.pdf",h=6,w=8)
for (i in rownames(snxRpkm) ) {

par(mar=c(6,6,4,2),cex.axis=2,cex.lab=2,cex.main=2)
plot(as.numeric(pd$DAY), snxRpkm[i,], 
	pch = 21, bg=pd$lineCol,
	cex=2, xlab="Day",
	ylab="log2(Exprs + 1)",
	main = foi[match(i,rownames(foi)),'Tag'], xaxt="n")
axis(1, at=1:length(levels(pd$DAY)), labels = levels(as.factor(pd$DAY)))
#lines(as.numeric(pd$DAY)[lineInd], snxRpkm[i,lineInd], lwd=1.5, col=pal[5])
legend("topright", levels(factor(pd$LINE)),
	pch = 15, col = pal, cex=1)
	
}	
dev.off()

