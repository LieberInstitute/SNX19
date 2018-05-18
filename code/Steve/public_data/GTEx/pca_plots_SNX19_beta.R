load('/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/gtex_eqtl_beta_pca.rda')
library(jaffelab)
library(ggrepel)
library(ggplot2)
library(ggplot2)
theme_set(theme_bw(base_size=12) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5, face='italic') ))
				 
				 
PC_Vars <- getPcaVars(gtex_snx19_pca)
gtex_snx19_pca$x

dat <- data.frame(gtex_snx19_pca$x)
dat$Tissue <- rownames(dat)
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/plots/GTEX_PCA_Plot.pdf')
a <- ggplot(data= dat, aes(x=PC1, y=PC2) ) +geom_point() + geom_text_repel(data = subset(dat,abs(PC1)>1000 | abs(PC2)>1000 ), aes(label=Tissue)) + labs(x=paste0("PC1: ",PC_Vars[1],"%"),y=paste0("PC2:",PC_Vars[2],"%")  )
print(a)

b <- ggplot(data= dat, aes(x=PC2, y=PC3) ) +geom_point() + geom_text_repel(data = subset(dat,abs(PC2)>1000 | abs(PC3)>200 ), aes(label=Tissue)) + labs(x=paste0("PC2: ",PC_Vars[2],"%"),y=paste0("PC3:",PC_Vars[3],"%")  )
print(b)

c <- ggplot(data= dat, aes(x=PC3, y=PC4) ) +geom_point() + geom_text_repel(data = subset(dat,abs(PC3)>1000 | abs(PC4)>200 ), aes(label=Tissue)) + labs(x=paste0("PC3: ",PC_Vars[3],"%"),y=paste0("PC4:",PC_Vars[4],"%")  )
print(c)

d <- ggplot(data= dat, aes(x=PC4, y=PC5) ) +geom_point() + geom_text_repel(data = subset(dat,abs(PC4)>200 | abs(PC5)>200 ), aes(label=Tissue)) + labs(x=paste0("PC4: ",PC_Vars[4],"%"),y=paste0("PC5:",PC_Vars[5],"%")  )
print(d)
dev.off()