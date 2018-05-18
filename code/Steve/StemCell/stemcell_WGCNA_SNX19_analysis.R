libary(corrplot)
load('dcl01/lieber/ajaffe/lab/libd_stem_timecourse/SNX19/stemcell_SNX19_features.rda')

#cor(datME,t(snxRpkm) )

res = Hmisc::rcorr(as.matrix(cbind(datME,t(snxRpkm) )))
cormat=res$r[colnames(datME),rownames(snxRpkm)]
pmat = p.mat = res$P[colnames(datME),rownames(snxRpkm)]
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_stemcell_WGCNA_11Modules.pdf',height=11,width=8.5)
corrplot(cormat ,p.mat=pmat , sig.level = .05)
pmat[cormat<0] =1
cormat[cormat<0]=0
corrplot(corr=cormat, p.mat=pmat, sig.level = .05)

dev.off()
