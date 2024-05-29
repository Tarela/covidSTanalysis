library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(harmony)
library(gplots)
library(pheatmap)
#S1A_image <- Read10X_Image(
#  image.dir = "/project/zanglab_project/sh8tv/SunLab/Data/visium/mouseLungCD8/10x_analysis_6433-HN/Sample_6433-HN-S1-A-GEX_TCAGAACT-GTGGACTG",
#  image.name = "tissue_lowres_image.png",
#  filter.matrix = TRUE
#)
covid1 <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/humanLung_spaceranger/covid1")
covid2 <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/humanLung_spaceranger/covid2")
covid3 <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/humanLung_spaceranger/covid3")
ctrl1  <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/humanLung_spaceranger/ctrl1")
ctrl2  <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/humanLung_spaceranger/ctrl2")
ctrl3  <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/humanLung_spaceranger/ctrl3")


sampleName <- c("covid1","covid2","covid3","ctrl1","ctrl2","ctrl3")
sampleName2 <- c("21-035","BB-22001","LIBB-0093","cc0005-19m3","cc003-20","cc07_18")
sampleCol <- c("#FF5555","#D43F00","#8B0000","#6495ED","#4169E1","#003366")

pdf(file="sample_legend.pdf")
plot(1,1,type="n",axes=F,xlab="",ylab="")
legend("center",legend=paste0(sampleName2, ":", sampleName),col=sampleCol,bty="n",cex=2,lwd=5)
dev.off()


### QC

pdf(file="sampleQC_outlier.pdf",width=12,height=4)
par(mfrow=c(1,3),mar=c(4,4,2,2))
barplot(c(length(covid1$nCount_Spatial), 
	    		length(covid2$nCount_Spatial),
	    		length(covid3$nCount_Spatial),
	    		length(ctrl1$nCount_Spatial),
	    		length(ctrl2$nCount_Spatial),
	    		length(ctrl3$nCount_Spatial)),
	    names=sampleName,ylab="#spot",main="#spot per spot",col=sampleCol)
boxplot(log10(as.numeric(covid1$nCount_Spatial)+1), 
	    	log10(as.numeric(covid2$nCount_Spatial)+1),
	    	log10(as.numeric(covid3$nCount_Spatial)+1),
	    	log10(as.numeric(ctrl1$nCount_Spatial)+1),
	    	log10(as.numeric(ctrl2$nCount_Spatial)+1),
	    	log10(as.numeric(ctrl3$nCount_Spatial)+1),
	    names=sampleName,ylab="log10 reads count",main="#reads per spot",col=sampleCol, outline=T)

boxplot(log10(as.numeric(covid1$nFeature_Spatial)+1), 
	      log10(as.numeric(covid2$nFeature_Spatial)+1),
	      log10(as.numeric(covid3$nFeature_Spatial)+1),
	      log10(as.numeric(ctrl1$nFeature_Spatial)+1),
	      log10(as.numeric(ctrl2$nFeature_Spatial)+1),
	      log10(as.numeric(ctrl3$nFeature_Spatial)+1),
	    names=sampleName,ylab="log10 covered gene",main="#genes per spot",col=sampleCol, outline=T)
abline(h=log10(2000),lwd=2)
dev.off()


#pdf(file="sampleQC.pdf",width=12,height=4)
#par(mfrow=c(1,3),mar=c(4,4,2,2))
#barplot(c(length(S1A$nCount_Spatial), 
#	    length(S1D$nCount_Spatial),
#	    length(S2A$nCount_Spatial),
#	    length(S2D$nCount_Spatial)),
#	    names=sampleName,ylab="#spot",main="#spot per spot",col=sampleCol)
#boxplot(log10(as.numeric(S1A$nCount_Spatial)), 
#	    log10(as.numeric(S1D$nCount_Spatial)),
#	    log10(as.numeric(S2A$nCount_Spatial)),
#	    log10(as.numeric(S2D$nCount_Spatial)),
#	    names=sampleName,ylab="log10 reads count",main="#reads per spot",col=sampleCol, outline=F)
#
#boxplot(log10(as.numeric(S1A$nFeature_Spatial)), 
#	    log10(as.numeric(S1D$nFeature_Spatial)),
#	    log10(as.numeric(S2A$nFeature_Spatial)),
#	    log10(as.numeric(S2D$nFeature_Spatial)),
#	    names=sampleName,ylab="log10 covered gene",main="#genes per spot",col=sampleCol, outline=F)
#dev.off()

### QC filtering
# nFeature >= 2000
covid1_highQ <- covid1[,which(covid1$nFeature_Spatial>=2000)]
covid2_highQ <- covid2[,which(covid2$nFeature_Spatial>=2000)]
covid3_highQ <- covid3[,which(covid3$nFeature_Spatial>=2000)]
ctrl1_highQ <- ctrl1[,which(ctrl1$nFeature_Spatial>=2000)]
ctrl2_highQ <- ctrl2[,which(ctrl2$nFeature_Spatial>=2000)]
ctrl3_highQ <- ctrl3[,which(ctrl3$nFeature_Spatial>=2000)]

#saveRDS(covid1_highQ, file="covid1_highQ.rds")
#saveRDS(covid2_highQ, file="covid2_highQ.rds")
#saveRDS(covid3_highQ, file="covid3_highQ.rds")
#saveRDS(ctrl1_highQ, file="ctrl1_highQ.rds")
#saveRDS(ctrl2_highQ, file="ctrl2_highQ.rds")
#saveRDS(ctrl3_highQ, file="ctrl3_highQ.rds")
#
#covid1_highQ <- readRDS(file="covid1_highQ.rds")
#covid2_highQ <- readRDS(file="covid2_highQ.rds")
#covid3_highQ <- readRDS(file="covid3_highQ.rds")
#ctrl1_highQ <- readRDS(file="ctrl1_highQ.rds")
#ctrl2_highQ <- readRDS(file="ctrl2_highQ.rds")
#ctrl3_highQ <- readRDS(file="ctrl3_highQ.rds")

pdf(file="spotNum_highQ.pdf")
#par(mfrow=c(1,3),mar=c(4,4,2,2))
barplot(c(length(covid1_highQ$nCount_Spatial), 
	    		length(covid2_highQ$nCount_Spatial),
	    		length(covid3_highQ$nCount_Spatial),
	    		length(ctrl1_highQ$nCount_Spatial),
	    		length(ctrl2_highQ$nCount_Spatial),
	    		length(ctrl3_highQ$nCount_Spatial)),
	    names=sampleName,ylab="#spot",main="highQ #spot per spot",col=sampleCol)
dev.off()


summary(covid1_highQ@images$slice1@coordinates[,2:3])
summary(covid2_highQ@images$slice1@coordinates[,2:3])
summary(covid3_highQ@images$slice1@coordinates[,2:3])
summary(ctrl1_highQ@images$slice1@coordinates[,2:3])
summary(ctrl2_highQ@images$slice1@coordinates[,2:3])
summary(ctrl3_highQ@images$slice1@coordinates[,2:3])



# test position
XLYL <- c(0,130)
useimage <- matrix(rep(NA,130*130),nrow=130)


norm01 <- function(indata){
	zmax <-quantile(indata,0.99)
	zmin <- quantile(indata,0.01)
	indata[indata > zmax] <- zmax#quantile(indata,0.99)
	indata[indata < zmin] <- zmin
	outdata0 <- (indata - min(indata)) / (max(indata) - min(indata))
	outdata <- round(outdata0*100) + 1
	return(outdata)
}
norm01fix <- function(indata,zmin,zmax){
	#zmax <-quantile(indata,0.99)
	#zmin <- quantile(indata,0.01)
	indata[indata > zmax] <- zmax#quantile(indata,0.99)
	indata[indata < zmin] <- zmin
	outdata0 <- (indata - min(indata)) / (max(indata) - min(indata))
	outdata <- round(outdata0*100) + 1
	return(outdata)
}

spotFeature_plot <- function(Coordinate, Feature, log_scale, usecolor, M){
	XLYL <- c(0,130)

	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	if(log_scale == 1){
		feature <- log10(Feature[spot_names]+1)
	}else{
		feature <- Feature[spot_names]
	}
	norm01_feature <- norm01(feature)
	ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
	color_feature <- ColorRamp[norm01_feature]
	data_range <- paste0(min(Feature)," - ",max(Feature))
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab=data_range, ylab="")
	#return(list(x_pos, y_pos, color_feature))
}

spotFeature_plot_fix <- function(Coordinate, Feature, log_scale, usecolor, M,zmin,zmax){
	XLYL <- c(0,130)

	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	if(log_scale == 1){
		feature <- log10(Feature[spot_names]+1)
	}else{
		feature <- Feature[spot_names]
	}
	norm01_feature <- norm01fix(feature,zmin,zmax)
	ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
	color_feature <- ColorRamp[norm01_feature]
	data_range <- paste0(zmin," - ",zmax)
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab=data_range, ylab="")
	#return(list(x_pos, y_pos, color_feature))
}

pdf(file="spotFeature_nCount.pdf",height=9,width=13.5)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
spotFeature_plot(covid1_highQ@images$slice1@coordinates, covid1_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="covid1 log10 nCount")
spotFeature_plot(covid2_highQ@images$slice1@coordinates, covid2_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="covid2 log10 nCount")
spotFeature_plot(covid3_highQ@images$slice1@coordinates, covid3_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="covid3 log10 nCount")
spotFeature_plot(ctrl1_highQ@images$slice1@coordinates, ctrl1_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="ctrl1 log10 nCount")
spotFeature_plot(ctrl2_highQ@images$slice1@coordinates, ctrl2_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="ctrl2 log10 nCount")
spotFeature_plot(ctrl3_highQ@images$slice1@coordinates, ctrl3_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="ctrl3 log10 nCount")
dev.off()


pdf(file="spotFeature_nGene.pdf",height=9,width=13.5)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
spotFeature_plot(covid1_highQ@images$slice1@coordinates, covid1_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="covid1 log10 nGene")
spotFeature_plot(covid2_highQ@images$slice1@coordinates, covid2_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="covid2 log10 nGene")
spotFeature_plot(covid3_highQ@images$slice1@coordinates, covid3_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="covid3 log10 nGene")
spotFeature_plot(ctrl1_highQ@images$slice1@coordinates, ctrl1_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="ctrl1 log10 nGene")
spotFeature_plot(ctrl2_highQ@images$slice1@coordinates, ctrl2_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="ctrl2 log10 nGene")
spotFeature_plot(ctrl3_highQ@images$slice1@coordinates, ctrl3_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="ctrl3 log10 nGene")
dev.off()


#### highQ SCT

covid1_highQ_SCT <- SCTransform(covid1_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
covid2_highQ_SCT <- SCTransform(covid2_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
covid3_highQ_SCT <- SCTransform(covid3_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
ctrl1_highQ_SCT <- SCTransform(ctrl1_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
ctrl2_highQ_SCT <- SCTransform(ctrl2_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
ctrl3_highQ_SCT <- SCTransform(ctrl3_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)


#saveRDS(covid1_highQ_SCT, file="covid1_highQ_SCT.rds")
#saveRDS(covid2_highQ_SCT, file="covid2_highQ_SCT.rds")
#saveRDS(covid3_highQ_SCT, file="covid3_highQ_SCT.rds")
#saveRDS(ctrl1_highQ_SCT, file="ctrl1_highQ_SCT.rds")
#saveRDS(ctrl2_highQ_SCT, file="ctrl2_highQ_SCT.rds")
#saveRDS(ctrl3_highQ_SCT, file="ctrl3_highQ_SCT.rds")
#
#covid1_highQ_SCT <- readRDS(file="covid1_highQ_SCT.rds")
#covid2_highQ_SCT <- readRDS(file="covid2_highQ_SCT.rds")
#covid3_highQ_SCT <- readRDS(file="covid3_highQ_SCT.rds")
#ctrl1_highQ_SCT <- readRDS(file="ctrl1_highQ_SCT.rds")
#ctrl2_highQ_SCT <- readRDS(file="ctrl2_highQ_SCT.rds")
#ctrl3_highQ_SCT <- readRDS(file="ctrl3_highQ_SCT.rds")

allgenes <- intersect(intersect(intersect(intersect(intersect(rownames(covid1_highQ_SCT@assays$SCT@scale.data),rownames(covid2_highQ_SCT@assays$SCT@scale.data)),
																										rownames(covid3_highQ_SCT@assays$SCT@scale.data)), 
																					rownames(ctrl1_highQ_SCT@assays$SCT@scale.data)), 
																rownames(ctrl2_highQ_SCT@assays$SCT@scale.data)), 
											rownames(ctrl3_highQ_SCT@assays$SCT@scale.data))

markerGeneList_spotFeatures <- function(Glist, usename){
	dirname <- paste0("spotFeatures_keygenes/",usename)
	if (!file.exists(dirname)){
	    dir.create(dirname)
	}

	for(this_keygene_raw in unique(Glist) ){
		this_keygene <- toupper(this_keygene_raw)
		#if(this_keygene %in% allgenes){
		print(paste0(this_keygene))
		pdf(file=paste0(dirname,"/spotFeature_scaleExp_",this_keygene,".pdf"),height=9,width=14.5)
		par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
		if (this_keygene %in% rownames(covid1_highQ_SCT@assays$SCT@scale.data)){
			spotFeature_plot(covid1_highQ_SCT@images$slice1@coordinates, covid1_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("covid1 scaleExp ",this_keygene))
		}else{
			plot(1,1,type="n",axes=F,xlab="",ylab="")
		}
		if (this_keygene %in% rownames(covid2_highQ_SCT@assays$SCT@scale.data)){
			spotFeature_plot(covid2_highQ_SCT@images$slice1@coordinates, covid2_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("covid2 scaleExp ",this_keygene))
		}else{
			plot(1,1,type="n",axes=F,xlab="",ylab="")
		}
		if (this_keygene %in% rownames(covid3_highQ_SCT@assays$SCT@scale.data)){
			spotFeature_plot(covid3_highQ_SCT@images$slice1@coordinates, covid3_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("covid3 scaleExp ",this_keygene))
		}else{
			plot(1,1,type="n",axes=F,xlab="",ylab="")
		}
		if (this_keygene %in% rownames(ctrl1_highQ_SCT@assays$SCT@scale.data)){
			spotFeature_plot(ctrl1_highQ_SCT@images$slice1@coordinates, ctrl1_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("ctrl1 scaleExp ",this_keygene))
		}else{
			plot(1,1,type="n",axes=F,xlab="",ylab="")
		}
		if (this_keygene %in% rownames(ctrl2_highQ_SCT@assays$SCT@scale.data)){
			spotFeature_plot(ctrl2_highQ_SCT@images$slice1@coordinates, ctrl2_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("ctrl2 scaleExp ",this_keygene))
		}else{
			plot(1,1,type="n",axes=F,xlab="",ylab="")
		}
		if (this_keygene %in% rownames(ctrl3_highQ_SCT@assays$SCT@scale.data)){
			spotFeature_plot(ctrl3_highQ_SCT@images$slice1@coordinates, ctrl3_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("ctrl3 scaleExp ",this_keygene))
		}else{
			plot(1,1,type="n",axes=F,xlab="",ylab="")
		}
		dev.off()
		#}else{
	#		print(paste0(this_keygene," NOT in human"))
	#	}
	}
}

markerGeneList_AVEspotFeatures <- function(Glist, usename,zmin,zmax){

		this_keygene <- intersect(toupper(unique(Glist) ),allgenes)
		pdf(file=paste0("spotFeature_scaleExpAVE_",usename,".pdf"),height=9,width=14.5)
		par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
		spotFeature_plot_fix(covid1_highQ_SCT@images$slice1@coordinates, apply(covid1_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid1"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(covid2_highQ_SCT@images$slice1@coordinates, apply(covid2_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid2"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(covid3_highQ_SCT@images$slice1@coordinates, apply(covid3_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid3"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(ctrl1_highQ_SCT@images$slice1@coordinates, apply(ctrl1_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl1"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(ctrl2_highQ_SCT@images$slice1@coordinates, apply(ctrl2_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl2"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(ctrl3_highQ_SCT@images$slice1@coordinates, apply(ctrl3_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl3"),zmin=zmin,zmax=zmax)
		dev.off()
		#}else{
	#		print(paste0(this_keygene," NOT in human"))
	#	}
}

markerGeneList_AVEspotFeatures_x <- function(Glist, usename,zmin,zmax){

		this_keygene <- intersect(toupper(unique(Glist) ),allgenes)
		spotFeature_plot_fix(covid1_highQ_SCT@images$slice1@coordinates, apply(covid1_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid1"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(covid2_highQ_SCT@images$slice1@coordinates, apply(covid2_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid2"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(covid3_highQ_SCT@images$slice1@coordinates, apply(covid3_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid3"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(ctrl1_highQ_SCT@images$slice1@coordinates, apply(ctrl1_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl1"),zmin=zmin,zmax=zmax)
		spotFeature_plot_fix(ctrl2_highQ_SCT@images$slice1@coordinates, apply(ctrl2_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl2"),zmin=zmin,zmax=zmax)
		#spotFeature_plot_fix(ctrl3_highQ_SCT@images$slice1@coordinates, apply(ctrl3_highQ_SCT@assays$SCT@scale.data[this_keygene,],2,mean), usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl3"),zmin=zmin,zmax=zmax)
#		dev.off()
		#}else{
	#		print(paste0(this_keygene," NOT in human"))
	#	}
}



### new human cell type marker
CD8_glist <- c("CD8A", "CCL5", "GZMH", "IL32", "CD8B", "NKG7", "GZMA", "GZMK", "CCL4", "IFNG", "TRBC2", "GZMB", "GNLY", "PRF1")
CD4_glist <- c("IL7R", "CD2", "CD3E", "IL32", "LTB", "CD3D", "SPOCK2", "RORA", "TRAC", "CD96", "TRBC2", "ICOS", "CD3G", "CD40LG", "ITK", "CTLA4")
Bcell_glist <- c("MS4A1", "CD79A", "BANK1", "LTB", "CD79B", "IGHM", "VPREB3", "RALGPS2", "TNFRSF13C", "MEF2C", "CD37", "LINC00926", "SMIM14", "ARHGAP24", "IGHD", "POU2F2")
AM_glist <- c("IFI27", "FABP4", "SERPING1", "RBP4", "INHBA", "GPD1", "PCOLCE2", "CES1", "CD52", "IGFBP2", "HP", "LDHB", "NUPR1", "IGF1", "MME", "CCL24", "PPARG", "APOC1", "CFD", "ITIH5", "SCD")
MDM_glist <- c("SPP1", "FCGR2B", "TYMP", "FGL2", "TMEM176A", "CD14", "LGMN", "IFITM3", "ZFP36L1", "C15orf48", "TIMP1", "MS4A6A", "TMEM176B", "CORO1A", "S100A8", "FCN1", "VCAN", "CCR2", "RNASE1", "EMP1", "CD14")
cDC_glist <- c("CCL17", "HLA-DPB1", "HLA-DQA1", "HLA-DPA1", "HLA-DRA", "HLA-DQB1", "HLA-DRB1", "CD74", "HLA-DRB5", "CST3", "HLA-DQA2", "SGK1", "MS4A6A", "S100B", "GPR183", "HLA-DMA")
Neutrophils_glist <- c("S100A8", "S100A9", "CXCL8", "S100A12", "G0S2", "NAMPT", "IFITM2", "FCGR3B", "SLC25A37", "CSF3R", "S100P", "SOD2", "BASP1", "FPR1", "MNDA", "MXD1")
NK_cells_glist <- c("GNLY", "NKG7", "PRF1", "GZMB", "CCL3", "KLRD1", "CCL4", "FGFBP2", "SPON2", "CTSW", "CTSW", "CD7", "KLRF1", "CLIC3", "CST7", "KLRB1")
pDCs_glist <- c("PLD4", "LILRA4", "TSPAN13", "CLEC4C", "EGLN3", "PPP1R14A", "SPIB", "MAP1A", "TPM2", "SCT", "LRRC26", "PTCRA", "DNASE1L3", "SEC11C")
Plasma_cells_glist <- c("IGKC", "JCHAIN", "IGHA1", "IGKV4-1", "IGHG1", "IGKV1-12", "IGHG3", "IGLV3-1", "MZB1", "IGHA2", "IGHG4", "IGLL5", "IGLV6-57", "SERPINF1")
Basal_cells_glist <- c("KRT17", "KRT15", "S100A2", "KRT5", "MMP1", "MT1X", "AQP3", "IGFBP6", "IL33", "KRT19", "IGFBP2", "LGALS7B", "LGALS7", "PERP")
Ciliated_glist <- c("CAPS", "C20orf85", "C9orf24", "TPPP3", "TMEM190", "FAM183A", "C1orf194", "RSPH1", "SNTN", "C5orf49", "CETN2", "PIFO", "C9orf116", "ODF3B", "MORN2")
Goblet_glist <- c("BPIFB1", "MSMB", "SCGB3A1", "LCN2", "TFF3", "BPIFA1", "MUC5B", "SAA1", "SAA2", "RARRES1", "SLPI", "TSPAN8", "VMO1")
Serous_glist <- c("LYZ", "PRR4", "C6orf58", "LTF", "AC020656.1", "PRB3", "PIP", "CRISP3", "S100A1", "AC006518.7", "S100B")
Transitional_glist <- c("PRSS2", "SPINK1", "TAGLN", "FN1", "IL32", "COL1A1", "CALD1", "FRZB", "MMP7", "ITGB6")
AT1_glist <- c("AGER", "CAV1", "EMP2", "CLDN18", "RTKN2", "MYL9", "NTM", "CCND2", "SCEL", "LMO7")
AT2_glist <- c("SFTPA1", "SFTPC", "SFTPA2", "PGC", "SFTPD", "NAPSA", "LAMP3", "C11orf96", "MFSD2A", "NPC2", "LRRK2", "CTSH", "PEBP4", "HHIP", "ABCA3")
summary3_glist <- c("Hopx","Clic5","Cldn18","Sftpa1","Abca3","CCL5","GZMK","SPP1","VCAN","Ces1","FABP4","CD52")
### new human gene lists
# Assuming your function is defined like this:

# Create a list of gene lists
# Create a list of gene lists with modified names
all_gene_lists <- list(CD8_humanNew = CD8_glist, CD4_humanNew = CD4_glist, 
                       Bcell_humanNew = Bcell_glist, AM_humanNew = AM_glist, 
                       MDM_humanNew = MDM_glist, cDC_humanNew = cDC_glist, 
                       Neutrophils_humanNew = Neutrophils_glist, 
                       NK_cells_humanNew = NK_cells_glist, pDCs_humanNew = pDCs_glist, 
                       Plasma_cells_humanNew = Plasma_cells_glist, 
                       Basal_cells_humanNew = Basal_cells_glist, 
                       Ciliated_humanNew = Ciliated_glist, Goblet_humanNew = Goblet_glist, 
                       Serous_humanNew = Serous_glist, Transitional_humanNew = Transitional_glist, 
                       AT1_humanNew = AT1_glist, AT2_humanNew = AT2_glist,newAdd_humanNew=summary3_glist)

# Use lapply to apply the function to each gene list
results <- lapply(names(all_gene_lists), function(x) {
  markerGeneList_spotFeatures(all_gene_lists[[x]],  x)
})

markerGeneList_spotFeatures(all_gene_lists[["newAdd_humanNew"]],  "newAdd_humanNew")


############ marker gene plots
CD8_glist <- c("Cd8a","Cd8b","Trbc1","Nkg7","Cd3d","Il7r","Itgae","Ccl5")
NK_glist <- c("Klrd1","Klra7","Klrc1","Klrk1")
CD4_glist <- c("Cd3d","Cd4","Trdc","Trbc1","Trbc2","Cd44","Tbx21","Ifng","Gata3","Stat4","Il4","Il5","Rorc","Rora","Il17a","Foxp3","Bcl6","Pdcd1","Izumo1r","Cxcr5")
B_glist <- c("Ptprc","Cd19","Apex1","Mzb1","Ly6a","Prdm1","Irf4","Ighm","Ighd")
Mono_glist <- c("Spn","Nr4a1","Cx3cr1","Cd36","Cd300e","Cd300a")
#AM_glist <- c("Fabp4","Pparg","Siglecf","Ccl6","Atp6v0d2","Mrc1","Fabp1","Ctsd","Abcg1") #alveolar macrophages
AM_glist <- c("Sparc","cd36","Flt1","Fabp5","Ca4") #alveolar macrophages
MDM_glist <- c("Junb","Cd14","Stat1","Il1rn","Ccl4","S100a6","Apoe","Mafb")
IM_glist <- c("C1qb","C1qc","C1qa","Ms4a7","Pf4","Aif1","Ccl12","Cd163","CLEC10A") #interstitial
DC_glist <- c("Irf4","Ccr7","Itgae","Xcr1","Zbtb46")
epi_glist <- c("Epcam","Sftpc","Sftpd","Pecam1")
Krt_glist <- c("Krt5","Krt8","Krt17","Tp53","Cldn4")
#Krt8_glist <- c("Krt8","Cldn4","Fosl1","Areg","Lgals3","Clu","Krt18","Sprr1a","Krt7","Cryab","Krt19")
#Krt5_glist <- c("Krt5","Trp63","Krt17","Krt14","Aqp3","Ngfr","Ly6d")
AT1_glist <- c("Ager","Akap5","Clic5","Rtkn2","Vegfa","Cldn18","Emp2","Spock2","Timp3","Lmo7","Col4a3","IGFBP2","Hopx")
AT2_glist <- c("Sftpa1","Slc34a2","Sftpb","Hc","Chi3l1","Lamp3","Sfta2","S100g","Sftpd")
activatedAT2_glist <- c("Lcn2","Sftpd","Chi3l1","Scd","Cxcl15","Lyz","Ager","Napsa","Lrg1","Cystm1","Rgcc")
summary_glist <- c("Cldn18","Clic5","Emp2","Vegfa","Sftpc","Lamp3","Hc","S100g","Slc34a2","Krt14","Aqp3","Krt17","Krt5","Cldn4","Krt8","Krt19","Cd8a","Cd8b","itgae","Cd3g","C1qa","C1qb","Ms4a7","Cx3cr1","Il1b","Il6","Ndrg1","Ifng","Tnf","Lum","ACTA2","Col1a1","Scgb3a2","scgb1a1","SCGB1A1","Foxj1","Pecam1","Plvap","Gpihbp1","Vwf","Tp63","ITGB4","Aqp5","Krt15","Csf1","Ccl2","Csf1r")
summary2_glist <- c("CD8B","GZMB","GZMK","RUNX3","ZNF683","CD69","ITGAE","IL7R","IL21","PRDM1","IGHA1","Tox")
#IL1bSig_glist <- c("Ikbkb","Il1a","Il1b","Il1r1","Il1rap","Il6","Irak1","Irak2","Irak3","Irak4","Mapk3","Myd88")
#inflammasome_glist <- c("Casp1","Ddx3x","Dhx33","Gsdmd","Aim2","Nlrp3","Pycard","Nlrp1a","Nlrp1b")
inflammasome_glist <- c("Casp1","Aim2","Pycard","Casp4")
posRegTCR_glist <- c("Ada","Bcl10","Card11","Ccr7","Cd81","Cd226","Cyld","Ikbkg","Kcnn4","Lipa","Nectin2","Prkd2","Rab29","Rela","Rps3","Tespa1","Trat1","Usp12","Usp46")
inflamComplex_glist <- c("Aim2","Casp1","Casp4","Ddx3x","Dhx33","Gsdmd","Mefv","Naip","Nlrc4","Nlrp1","Nlrp3","Nlrp6","Nlrp9","Pycard")
#AE_glist <- c("Cldn18","Clic5","Emp2","Vegfa","Sftpc","Lamp3","Hc","S100g","Slc34a2")
AE_glist <- c("Sftpc","Clic5","Emp2","Hopx","Spock2","Sftpa1","Pecam1","Cldn18")
TGFb_glist <- c("Gzmk","Eomes","Nkg7","Tox","Ccl2","Ccr2","Cxcl13","Lta","Ltb","Hif1a","Ifng","Ifna","Ifnb1","Cxcl16")
IL1bSig_glist <- c("Irak4","Il1rap","Il1b","Irak2","Il1r1","Myd88")

Aberrant_basaloid_glist <- c("Epcam","Chd1","Vim","Fn1","Col1a1","Chd2","Tnc","Vcan","Pcp4","Cux2","Spink1","Prss2","Cpa6","Ctse","Mmp7","Mdk","Gdf15","Ptgs2","Slco2a1","Ephb2","Itgb8","Itgav","Itgb6","Tgfbl","Kcnn4","Kcnq5","Kcns3","Cdkn1a","Cdkn2a","Cdkn2b","Ccnd1","Ccnd2","Mdm2","Hmga2","Ptchd4","Ociad2","Tp63","Krt17","Lab3","Lamc2")
Profibrotic_glist <- c("Csf1","Adamts9","Cdkn1a","Timp3","Col4a1","Col4a2","Slco2a1","Plvap","Mmp14","Fbn1","Fn1","Col5a1","Col1a2","Col1a1","Col3a1","Sparc","Sparcl1","Col6a1","Cd74","Mfap5","Spp1","Chil1")
ADI_DATP_PATS_glist <- c("Anxa2","Bcl2l1","Itga3","Jun","Krt8","Krt19","F3","Klf6","Msn","Tapbp","Cdkn1a","Icam1","Timp3","Anxa1","Krt18","Tnfrsf1a","Sdc4","Cd63","Npc2","Myh9","Actg1","B2m","Jund","F11r")
Wound_healing_glist <- c("Actg1","Adra2a","Ano6","Anxa1","Apoh","Arfgef1","Ccl2","Ccn4","Cd36","Cldn1","Cldn3","Cldn4","Cldn13","Clec7a","Cpb2","Cxcr4","Ddr2","Dmtn","Duox1","Duox2","Emilin2","Enpp4","F2","F2r","F3","F7","F12","Fermt1","Fermt2","Foxc2","Hbegf","Hif1a","Hmgb1","Hpse","Hras","Hrg","Insl3","Itgb1","Kank1","Klrh1","Mtor","Mylk","Nfe2l2","Plat","Plau","Plg","Prdx2","Prkce","Ptger4","Ptk2","Reg3a","Reg3g","Rreb1","S100a9","Serpine1","Serpinf2","Smoc2","St3gal4","Tbxa2r","Thbs1","Vegfb","Vtn","Vwf","Xbp1")
Epi_proliferation_glist <- c("B4galt1","Cldn1","Cxadr","Eppk1","Fzd7","Jaml","Lrg1","Mmp12","Wnt7a")
Fibrosis_glist <- c("Stn1","Dsp","Bmp7","Calca","Cebpb","Ccr3","Ccr2","Csf2","Csf3","Edn1","Egf","Eln","Fgf1","Fgf2","Fgf7","Ccn2","Hgf","Hmox1","Igf1","Il12b","Il13","Il1b","Il4","Il5","Il6","Smad7","Cma1","Mecp2","Mmp2","Mmp9","Mt2","Nfe2l2","Pdgfa","Pdgfb","Plau","Ptx3","Ccl11","Ccl2","Ccl3","Ccl4","Ccl5","Cxcl15","Cxcl2","Sftpa1","Sftpc","Skil","Spp1","Tert","Tgfa","Tgfb1","Timp1","Tnf","Dpp9","Elmod2","Rtel1","Atp11a","Fam13a","Cysltr2","Parn","Muc5b")

MDM_original <- c("VCAN","FCN1","S100A8","MS4A6A","IFITM3","CD14","IFITM2","S100A9","TNFRSF1B","CD99","TGFBI","ITGAM","HIF1A","CD36","S100A10","CLEC5A","TNFSF10","JUNB","CCR1","IL7R","CD44","NFKBIA","TREM2","ITGB2","ITGA4","SLC43A3","IL10RA","IL17RA","FN1","PRCP","JAK1","STAT1","IFNGR2","GBP1")
AM_Original <- c("CD46","CLDN7","FABP5","LYZ","CD59","SPARC","ACOT2","BHLHE41","ABHD5","PEBP1","LY6E","PPARG","TREM1","S100A13","FABP4")
Profibrotic_Original <- c("CCL2","CCL4","CD14","CD44","CSF1R","FN1","HIF1A","HLA-A","HLA-DQA2","IL1RN","ITGAX","MERTK","S100A8","S100A9","STAT1","TGFBI","ZNF385A")
MDM_Knowledge <- c("Ly6c2","Ccr2","Cd14","Fn1","Vcan","IL1b","S100a6","S100a4","Apoe","Mafb","C5ar1","Ms4a4c","Csf1r")
AM_Knowledge <- c("Chil3","Fabp5","Fabp4","Flt1","Pparg","Siglecf","Car4","Ear1","Krt79","Cd169","Tcf7l2","Abcg1","Bach2","Bhlhe41")
IM_Knowledge <- c("C1qa","C1qb","C1qc","Pf4","C5ar1","Apoe","Cd14","Csf1r","Mafb","Cxc3r1","Cd163","Siglec1","Calr","Clec4n","Mpeg1","Lpl")


### human new gene lists
CD8_Hglist <- c("CD8A","CCL5","GZMH","CD3D","IL32","NKG7","GZMA","CD3G","CD3E","GZMK","CD2","CCL4","IFNG","TRBC2","TRGC2")
CD4_Hglist <- c("IL7R","CD2","CD3E","IL32","LTB","CD3D","SPOCK2","RORA","TRAC","CD96","TRBC2","ICOS","CD3G","CD40LG","ITK","CTLA4")
B_Hglist <- c("MS4A1","CD79A","BANK1","LTB","CD79B","IGHM","VPREB3","RALGPS2","TNFRSF13C","MEF2C","CD37",,"SMIM14","ARHGAP24","IGHD","POU2F2")
AM_Hglist <- c("IFI27","FABP4","SERPING1","RBP4","INHBA","GPD1","PCOLCE2","CES1","CD52","IGFBP2","HP","LDHB","NUPR1","IGF1","MME","CCL24")
MDM_Hglist <- c("SPP1","FCGR2B","TYMP","FGL2","TMEM176A","CD14","LGMN","IFITM3","ZFP36L1","TIMP1","MS4A6A","TMEM176B","CORO1A","S100A8","VCAN","CCL2","RNASE1")
AT1_Hglist <- c("AGER","CAV1","EMP2","CLDN18","RTKN2","MYL9","NTM","CLIC3","CAV2","TNNC1","ANXA3","AQP4","ANKRD29","HOPX","SCEL","LMO7")
AT2_Hglist <- c("SFTPA1","SFTPC","SFTPA2","PGC","SFTPD","NAPSA1","LAMP3","LRRK2","CTSH","ABCA3","SERPINA1","PEBP41","HHIP","MFSD2A","NPC2")

cDC_Hglist <- c("CCL17","HLA-DPB11","HLA-DQA11","HLA-DPA1","HLA-DRA1","HLA-DQB11","HLA-DRB1","CD741","HLA-DRB5","CST3","HLA-DQA21","SGK1","MS4A6A","S100B","GPR1831","HLA-DMA1")
Neutrophils_Hglist <- c("S100A8","S100A9","CXCL8","S100A12","G0S2","NAMPT","IFITM2","FCGR3B","SLC25A37","CSF3R","S100P","SOD21","BASP1","FPR1","MNDA","MXD1")
NK_Hglist <- c("GNLY1","NKG71","PRF11","GZMB1","CCL3","KLRD11","CCL41","FGFBP21","SPON21","CD2471","CTSW1","CD73","KLRF1","CLIC31","CST71","KLRB11")
pDC_Hglist <- c("PLD41","LILRA4","TSPAN132","CLEC4C","EGLN3","PPP1R14A1","SPIB2","MAP1A","TPM2","SMPD3","SCT","LRRC26","PTPRS","LAMP5","DNASE1L3","PTCRA")
Plasma_Hglist <- c("IGKC1","JCHAIN1","IGHA1","IGKV4-1","IGHG1","IGKV1-12","IGHG3","IGLV3-1","MZB11","IGHG2","IGHA2","IGHG4","IGLL5","DERL31","IGLV6-57","SEC11C")
Basal_Hglist <- c("KRT17","KRT15","S100A2","KRT5","MMP1","MT1X","AQP3","IGFBP6","IL33","DST1","KRT191","IGFBP2","LGALS7B","LGALS7","PERP","SERPINF1")
Ciliated_Hglist <- c("CAPS","C20orf85","C9orf24","TPPP3","TMEM190","FAM183A","C1orf194","RSPH1","SNTN","C11orf88","C5orf49","CETN2","PIFO","C9orf116","MORN2","ODF3B")
Goblet_Hglist <- c("BPIFB1","MSMB","SCGB3A11","LCN2","TFF3","BPIFA1","SCGB1A11","MUC5B","SAA11","WFDC2","SAA21","S100P1","RARRES11","SLPI","TSPAN8","VMO1")
Serous_Hglist <- c("LYZ","PRR4","C6orf58","LTF1","AC020656.1","PRB3","ZG16B","PIP1","PRH2","TCN1","CRISP3","S100A1","LPO","FRZB","AC006518.7","S100B")
Transitional_Hglist <- c("PRSS2","SPINK1","TAGLN1","FN11","IL321","COL1A1","CALD11","PCP4","PMEPA11","PTGS22","TMSB102","S100A102","ITGB62","MMP72","KRT71","TM4SF12","CTSE3","IFI272")
###


CD8_refine <- c("CD8A","CCL5","GZMH","GZMA","GZMK","GZMB","GNLY","PRF1","TRBC2","NKG7")
AM_refine <- c("CES1","GPD1","INHBA","MME","APOC1","CD52","IFI27","PCOLCE2","PPARG","SCD","FABP4")
MDM_refine <- c("CD14","FCGR2B","RNASE1","S100A8","SPP1","MS4A6A","TYMP","VCAN")
Krt_refine <- c("KRT8","KRT5","KRT17","CALD1","COL1A1","FN1","MMP7","ITGB6","AQP3","IL33","KRT19","MMP1","S100A2")
AE_refine <- c("AGER","CA1","CLDN18","ABCA3","LAMP3","SFTPA1","SFTPD","SFTPC","CLIC5","AQP4")


#### new normalization
covid1_countMat <- covid1_highQ@assays$Spatial@counts
covid2_countMat <- covid2_highQ@assays$Spatial@counts
covid3_countMat <- covid3_highQ@assays$Spatial@counts
ctrl1_countMat <- ctrl1_highQ@assays$Spatial@counts
ctrl2_countMat <- ctrl2_highQ@assays$Spatial@counts


covid1_countMatSCT <- covid1_highQ_SCT@assays$SCT@scale.data
covid2_countMatSCT <- covid2_highQ_SCT@assays$SCT@scale.data
covid3_countMatSCT <- covid3_highQ_SCT@assays$SCT@scale.data
ctrl1_countMatSCT <- ctrl1_highQ_SCT@assays$SCT@scale.data
ctrl2_countMatSCT <- ctrl2_highQ_SCT@assays$SCT@scale.data

idx <- Reduce(intersect, list(rownames(covid1_countMatSCT), 
	                          rownames(covid2_countMatSCT),
	                          rownames(covid3_countMatSCT),
	                          rownames(ctrl1_countMatSCT),
	                          rownames(ctrl2_countMatSCT)))
colnames(covid1_countMatSCT)<-paste0("covid1_",colnames(covid1_countMatSCT))
colnames(covid2_countMatSCT)<-paste0("covid2_",colnames(covid2_countMatSCT))
colnames(covid3_countMatSCT)<-paste0("covid3_",colnames(covid3_countMatSCT))
colnames(ctrl1_countMatSCT)<-paste0("ctrl1_",colnames(ctrl1_countMatSCT))
colnames(ctrl2_countMatSCT)<-paste0("ctrl2_",colnames(ctrl2_countMatSCT))

cbSCT <- cbind(covid1_countMatSCT[idx,],
			   covid2_countMatSCT[idx,],
			   covid3_countMatSCT[idx,],
			   ctrl1_countMatSCT[idx,],
			   ctrl2_countMatSCT[idx,]
	           	             )

usegenes <- c(CD8_refine,
AM_refine ,
MDM_refine,
Krt_refine,
AE_refine )
cor_usegenes <- cor(t(cbSCT[usegenes,]))
library(gplots)
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(from = -1, to = 1, length.out = 101)

pdf(file="cor_keygenes_SCT_heatmap.pdf",width=10,height=10)
heatmap.2(cor_usegenes, trace="none",margins=c(5,5),main="combine",col=my_palette,breaks=breaks,scale="none",dendrogram="none")
dev.off()

pdf(file="cor_CD8refine_SCT_heatmap.pdf",width=10,height=10)
heatmap.2(cor(t(cbSCT[CD8_refine,])), trace="none",margins=c(5,5),col=my_palette,breaks=breaks,scale="none",dendrogram="none",main="CD8")
dev.off()

pdf(file="cor_AMrefine_SCT_heatmap.pdf",width=10,height=10)
heatmap.2(cor(t(cbSCT[AM_refine,])), trace="none",margins=c(5,5),col=my_palette,breaks=breaks,scale="none",dendrogram="none",main="AM")
dev.off()

pdf(file="cor_MDMrefine_SCT_heatmap.pdf",width=10,height=10)
heatmap.2(cor(t(cbSCT[MDM_refine,])), trace="none",margins=c(5,5),col=my_palette,breaks=breaks,scale="none",dendrogram="none",main="MDM")
dev.off()

pdf(file="cor_Krtrefine_SCT_heatmap.pdf",width=10,height=10)
heatmap.2(cor(t(cbSCT[Krt_refine,])), trace="none",margins=c(5,5),col=my_palette,breaks=breaks,scale="none",dendrogram="none",main="Krt")
dev.off()

pdf(file="cor_AErefine_SCT_heatmap.pdf",width=10,height=10)
heatmap.2(cor(t(cbSCT[AE_refine,])), trace="none",margins=c(5,5),col=my_palette,breaks=breaks,scale="none",dendrogram="none",main="AE")
dev.off()

covid1_spotsumSCT <- apply(covid1_countMatSCT,2,sum)
covid2_spotsumSCT <- apply(covid2_countMatSCT,2,sum)
covid3_spotsumSCT <- apply(covid3_countMatSCT,2,sum)
ctrl1_spotsumSCT<-  apply(ctrl1_countMatSCT,2,sum)
ctrl2_spotsumSCT<-  apply(ctrl2_countMatSCT,2,sum)

boxplot(covid1_spotsumSCT, covid2_spotsumSCT, covid3_spotsumSCT, ctrl1_spotsumSCT, ctrl2_spotsumSCT,
	    names=c("covid1","covid2","covid3","ctrl1","ctrl2"),ylab="spotSum",main="SCT")



CD8_refine2 <- c("CD8A","CCL5","GZMH","GZMA","GZMK","GZMB","GNLY","PRF1","TRBC2","NKG7")
AM_refine2 <- c("GPD1","INHBA","MME","APOC1","CD52","PCOLCE2","PPARG","SCD","FABP4")
MDM_refine2 <- c("CD14","FCGR2B","RNASE1","S100A8","SPP1","MS4A6A","TYMP")
Krt_refine2 <- c("KRT8","KRT5","KRT17","MMP7","ITGB6","AQP3","KRT19","MMP1","S100A2")
AE_refine2 <- c("AGER","CLDN18","ABCA3","LAMP3","SFTPA1","SFTPD","SFTPC","CLIC5","AQP4")
CTkey_glist <- c("SFTPC","CD8A","KRT5","KRT8","KRT17","SFTPC","AGER","CD14","CX3CR1","CD68")
scaleMat <- cbSCT
TGFb_pathway <- intersect(rownames(scaleMat),c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP"))
IL1R_pathway <- intersect(rownames(scaleMat),c("IL1A","IL1B","IL1RN","TNF","IL6","TGFB1","IL1R1","IRAK1","IRAK1","TGFB2","IRAK3","RELA","NFKB1","IL1RAP","IL1RAP","IL1RAP","IL1RAP","MYD88","MYD88","MYD88","MYD88","IKBKB","IKBKB","RELA","RELA","CHUK","MAPK14","IKBKB","IRAK1","IRAK2","IFNB1","IL1RAP","JUN","MYD88","MAPK8","MAP2K3","MAP2K6","MAP3K7","TGFB2","TGFB3","MAP3K14","NFKB1","TRAF6","MAP3K1","TAB1","IRAK3","NFKBIA","RELA","IFNA1","IL1RAP","MAPK14","MAPK14","MAPK14","MAPK8"))
Inflammasome_pathway <- intersect(rownames(scaleMat),c("NFKB2","P2RX7","NLRC4","NLRP1","HSP90AB1","HMOX1","MEFV","PYCARD","NFKB1","PANX1","TXN","CASP1","PSTPIP1","APP","NLRP3","AIM2","SUGT1","BCL2L1","BCL2","RELA","TXNIP"))






#covid1_spotsum_norm <- (covid1_spotsum / median(covid1_spotsum))#/sd(covid1_spotsum)
#covid2_spotsum_norm <- (covid2_spotsum / median(covid2_spotsum))#/sd(covid2_spotsum)
#covid3_spotsum_norm <- (covid3_spotsum / median(covid3_spotsum))#/sd(covid3_spotsum)
#ctrl1_spotsum_norm <- (ctrl1_spotsum / median(ctrl1_spotsum))#/sd(ctrl1_spotsum)
#ctrl2_spotsum_norm <- (ctrl2_spotsum / median(ctrl2_spotsum))#/sd(ctrl2_spotsum)
#
#
#pdf(file="spotSumReadCount_boxplot.pdf",width=10,height=5)
#par(mfrow=c(1,2),mar=c(4,4,2,2))
#boxplot(covid1_spotsum, covid2_spotsum, covid3_spotsum, ctrl1_spotsum, ctrl2_spotsum,
#	    names=c("covid1","covid2","covid3","ctrl1","ctrl2"),ylab="spotSum",main="before Norm")
#boxplot(covid1_spotsum_norm, covid2_spotsum_norm, covid3_spotsum_norm, ctrl1_spotsum_norm, ctrl2_spotsum_norm,
#	    names=c("covid1","covid2","covid3","ctrl1","ctrl2"),ylab="spotSum",main="after Norm")
#dev.off()
#
#covid1_countMatNorm <- t(t(as.matrix(covid1_countMat)) / covid1_spotsum_norm)
#covid2_countMatNorm <- t(t(as.matrix(covid2_countMat)) / covid2_spotsum_norm)
#covid3_countMatNorm <- t(t(as.matrix(covid3_countMat)) / covid3_spotsum_norm)
#ctrl1_countMatNorm <- t(t(as.matrix(ctrl1_countMat)) / ctrl1_spotsum_norm)
#ctrl2_countMatNorm <- t(t(as.matrix(ctrl2_countMat)) / ctrl2_spotsum_norm)
#
#
#
#VlnPlot(cbdata, 
#        features = "nCount_SCT", 
#        group.by = "sample", 
#        pt.size = 0,
#        ncol = 2)


markerGeneList_spotFeatures(c(CD8_glist),"CD8")
markerGeneList_spotFeatures(c(MDM_glist),"MDM")
markerGeneList_spotFeatures(c(CD4_glist),"CD4")
markerGeneList_spotFeatures(c(B_glist),"B")
markerGeneList_spotFeatures(c(Mono_glist),"Mono")
markerGeneList_spotFeatures(c(AM_glist),"AM")
markerGeneList_spotFeatures(c(IM_glist),"IM")
markerGeneList_spotFeatures(c(DC_glist),"DC")
markerGeneList_spotFeatures(c(epi_glist),"epi")
markerGeneList_spotFeatures(c(Krt_glist),"Krt")

markerGeneList_spotFeatures(c(AT1_glist),"AT1")
markerGeneList_spotFeatures(c(AT2_glist),"AT2")
markerGeneList_spotFeatures(c(activatedAT2_glist),"activatedAT2")
markerGeneList_spotFeatures(c(summary_glist,summary2_glist),"summary")
markerGeneList_spotFeatures(c(IL1bSig_glist),"IL1bSig")
markerGeneList_spotFeatures(c(AT2_glist),"AT2")
markerGeneList_spotFeatures(c(inflammasome_glist),"inflammasome")
markerGeneList_spotFeatures(c(posRegTCR_glist),"posRegTCR")
markerGeneList_spotFeatures(c(inflamComplex_glist),"inflamComplex")
markerGeneList_spotFeatures(c(AE_glist),"AE")
markerGeneList_spotFeatures(c(TGFb_glist),"TGFb")


### human glist
markerGeneList_spotFeatures(c(CD8_Hglist),"CD8_human")
markerGeneList_spotFeatures(c(CD4_Hglist),"CD4_human")
markerGeneList_spotFeatures(c(B_Hglist),"B_human")
markerGeneList_spotFeatures(c(AM_Hglist),"AM_human")
markerGeneList_spotFeatures(c(MDM_Hglist),"MDM_human")
markerGeneList_spotFeatures(c(AT1_Hglist),"AT1_human")
markerGeneList_spotFeatures(c(AT2_Hglist),"AT2_human")

markerGeneList_spotFeatures(c(cDC_Hglist),"cDC_human")
markerGeneList_spotFeatures(c(Neutrophils_Hglist),"Neutrophils_human")
markerGeneList_spotFeatures(c(NK_Hglist),"NK_human")
markerGeneList_spotFeatures(c(pDC_Hglist),"pDC_human")
markerGeneList_spotFeatures(c(Plasma_Hglist),"Plasma_human")
markerGeneList_spotFeatures(c(Basal_Hglist),"Basal_human")
markerGeneList_spotFeatures(c(Ciliated_Hglist),"Ciliated_human")
markerGeneList_spotFeatures(c(Goblet_Hglist),"Goblet_human")
markerGeneList_spotFeatures(c(Serous_Hglist),"Serous_human")
markerGeneList_spotFeatures(c(Transitional_Hglist),"Transitional_human")


markerGeneList_AVEspotFeatures(c(CD8_Hglist),"CD8_human")
markerGeneList_AVEspotFeatures(c(CD4_Hglist),"CD4_human")
markerGeneList_AVEspotFeatures(c(B_Hglist),"B_human")
markerGeneList_AVEspotFeatures(c(AM_Hglist),"AM_human")
markerGeneList_AVEspotFeatures(c(MDM_Hglist),"MDM_human")
markerGeneList_AVEspotFeatures(c(AT1_Hglist),"AT1_human")
markerGeneList_AVEspotFeatures(c(AT2_Hglist),"AT2_human")
markerGeneList_AVEspotFeatures(c(AT1_Hglist,AT2_Hglist),"AT_human")
markerGeneList_AVEspotFeatures(c(Krt_glist),"Krt")
markerGeneList_AVEspotFeatures(c(AE_glist),"AE")



markerGeneList_AVEspotFeatures(c(Krt_refine2),"Krtrefine2",-1.7,4)

pdf(file=paste0("spotFeature_scaleExpAVE_combine.pdf"),height=20,width=24)
par(bg="grey",mfrow=c(5,6),mar=c(4,4,2,2))
markerGeneList_AVEspotFeatures_x(c(CD8_Hglist),"CD8_human",-0.8,2)
markerGeneList_AVEspotFeatures_x(c(Krt_glist),"Krt",-1.7,4)
markerGeneList_AVEspotFeatures_x(c(AE_glist),"AE",-2,5)
markerGeneList_AVEspotFeatures_x(c(MDM_Hglist),"MDM_human",-1,3)
markerGeneList_AVEspotFeatures_x(c(AT1_Hglist,AT2_Hglist),"AT_human",-2,5)
dev.off()


cbdata <- merge(covid1_highQ, y=c(covid2_highQ,covid3_highQ,ctrl1_highQ,ctrl2_highQ),add.cell.ids=c("covid1","covid2","covid3","ctrl1","ctrl2"),project="Spatial")

usesample <- unlist(lapply(names(cbdata$orig.ident), trimsamplename))
cbdata$sample <- usesample
disease <- rep("covid", length(cbdata$sample))
disease[which(cbdata$sample %in% c("ctrl1","ctrl2"))] <- "ctrl"
cbdata$disease <- disease

cbdata <- SCTransform(cbdata, assay = "Spatial", verbose = FALSE)
DefaultAssay(cbdata) <- "SCT"
set.seed(1)
cbdata <- RunPCA(cbdata, npcs = 30, verbose = FALSE, assay="SCT")
cbdata <- FindNeighbors(cbdata, dims = 1:30)
cbdata <- FindClusters(cbdata, verbose = FALSE)
cbdata <- RunUMAP(cbdata, dims = 1:30)



#saveRDS(cbdata, file="cbdata.rds")
#cbdata <- readRDS(file="cbdata.rds")


pdf(file="human_cbdata_SCTraw_UMAP_features.pdf",width=12,height=4)
p1 <- DimPlot(cbdata, label = TRUE,group.by="sample", reduction = "umap") + ggtitle("sample")
p2 <- DimPlot(cbdata, label = TRUE,group.by="disease", reduction = "umap") + ggtitle("disease")
p3 <- DimPlot(cbdata, label = TRUE,group.by="seurat_clusters", reduction = "umap") + ggtitle("cluster")
#p3 <- DimPlot(cbdata.combined, reduction = "umap", group.by = "predicted.id")
p1+p2+p3
dev.off()


#### disease wise DEG
DefaultAssay(cbdata) <- "SCT"
condition_DEG <- FindMarkers(cbdata,group.by="disease", ident.1 = "covid", ident.2 = "ctrl")
top10_condition_DEG <- condition_DEG[order(condition_DEG[,"avg_log2FC"]),][c(1:10, (nrow(condition_DEG)-9):nrow(condition_DEG)),]
pdf(file="human_UMAP_conditionDEGtop10.pdf",width=12,height=15)
FeaturePlot(cbdata, features=rownames(top10_condition_DEG), reduction="umap")
dev.off()

write.table(condition_DEG, file="human_condition_DEG.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(rownames(condition_DEG)[which(condition_DEG[,"avg_log2FC"]>0)], file="human_condition_covidHigh_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(rownames(condition_DEG)[which(condition_DEG[,"avg_log2FC"]<0)], file="human_condition_ctrlHigh_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)

condition_DEG <- read.table("human_condition_DEG.txt",row.names=1,header=T)
condition_DEG_filter <- condition_DEG[which(abs(condition_DEG[,"avg_log2FC"])>=log2(1.5) & condition_DEG[,"p_val_adj"] < 0.01),]
write.table(rownames(condition_DEG_filter)[which(condition_DEG_filter[,"avg_log2FC"]>0)], file="human_condition_covidHigh_filter_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(rownames(condition_DEG_filter)[which(condition_DEG_filter[,"avg_log2FC"]<0)], file="human_condition_ctrlHigh_filter_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)


condition_DEG <- FindMarkers(cbdata,group.by="disease", ident.1 = "covid", ident.2 = "ctrl")


saveRDS(cbdata, file="cbdata_new.rds")

#allspots_covid_vs_ctrl_DEG<- FindMarkers(cbdata, ident.1="covid",ident.2="ctrl",group.by="disease")
#allspots_covid_vs_ctrl_DEG_filter <- allspots_covid_vs_ctrl_DEG[which(abs(allspots_covid_vs_ctrl_DEG[,"avg_log2FC"])>=log2(1.5) & allspots_covid_vs_ctrl_DEG[,"p_val_adj"] < 0.01),]
#write.table(allspots_covid_vs_ctrl_DEG_filter[order(allspots_covid_vs_ctrl_DEG_filter[,"avg_log2FC"]),], file=paste0("allspots_covidVSctrl_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)
#write.table(rownames(allspots_covid_vs_ctrl_DEG_filter)[which(allspots_covid_vs_ctrl_DEG_filter[,"avg_log2FC"]>0)],file="allspots_covidVSctrl_DEG_covidHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
#write.table(rownames(allspots_covid_vs_ctrl_DEG_filter)[which(allspots_covid_vs_ctrl_DEG_filter[,"avg_log2FC"]<0)],file="allspots_covidVSctrl_DEG_ctrlHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)





scaleMat <- cbSCT#cbdataHAR2_scale_allgenes@assays$SCT@scale.data

Krt_score  <- apply(scaleMat[toupper(Krt_refine2),],2,mean)
Cd8_score <- apply(scaleMat[toupper(CD8_refine2),],2,mean)
Ae_score   <- apply(scaleMat[toupper(AE_refine2),],2,mean)
Mdm_score <- apply(scaleMat[toupper(MDM_refine2),],2,mean)
Am_score <- apply(scaleMat[toupper(AM_refine2),],2,mean)

TGFb_score <- apply(scaleMat[toupper(TGFb_pathway),],2,mean)
IL1R_score <- apply(scaleMat[toupper(IL1R_pathway),],2,mean)
Inflammasome_score <- apply(scaleMat[toupper(Inflammasome_pathway),],2,mean)
IFNG_score <- apply(scaleMat[toupper(IFNG_pathway),],2,mean)
TNF_score <- apply(scaleMat[toupper(TNF_pathway),],2,mean)
IFNGTNF_score <- apply(scaleMat[toupper(c(IFNG_pathway,TNF_pathway)),],2,mean)
#Cd8_score2 <- apply(scaleMat[toupper(c("Cd8a","Cx3cr1","Trbc1","Nkg7","Cd3d","Il7r")),],2,mean)
#Cd8_score3 <- apply(scaleMat[intersect(CD8MKG_top10, rownames(scaleMat)),],2,mean)
#Cd8_score4 <- apply(scaleMat[intersect(CD8MKG_top20, rownames(scaleMat)),],2,mean)
#Cd8_score5 <- apply(scaleMat[intersect(CD8MKG_FCge2, rownames(scaleMat)),],2,mean)

pdf(file="celltypeScore_hist_new.pdf",width=12,height=8)
par(mfrow=c(2,3),mar=c(4,4,2,2))
#hist(rawMat["Cd8a",],n=200, main="Cd8a",xlab="normalized exp")
#hist(rawMat["Krt8",],n=200, main="Krt8",xlab="normalized exp")
hist(Krt_score ,n=200, main="Krt_score", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(Ae_score  ,n=200, main="Ae_score",  xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(Cd8_score,n=200, main="Cd8_score",xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(Mdm_score,n=200, main="Mdm_score",xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(Am_score,n=200, main="Am_score",xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
dev.off()

CD8A_exp <- scaleMat["CD8A",]
KRT5_exp <- scaleMat["KRT5",]
KRT8_exp <- scaleMat["KRT8",]
KRT17_exp <- scaleMat["KRT17",]
SFTPC_exp <- scaleMat["SFTPC",]
AGER_exp <- scaleMat["AGER",]
CD14_exp <- scaleMat["CD14",]
CX3CR1_exp <- scaleMat["CX3CR1",]
CD68_exp <- scaleMat["CD68",]

pdf(file="keygeneScore_hist_new.pdf",width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,2))
hist(CD8A_exp ,n=200, main="CD8A", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(KRT5_exp ,n=200, main="KRT5", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(KRT8_exp ,n=200, main="KRT8", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(KRT17_exp ,n=200, main="KRT17", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(CD14_exp ,n=200, main="CD14", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(CX3CR1_exp ,n=200, main="CX3CR1", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(CD68_exp ,n=200, main="CD68", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(SFTPC_exp ,n=200, main="SFTPC", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
hist(AGER_exp ,n=200, main="AGER", xlab="SCT normalized exp",xlim=c(-3,5))
abline(v=c(0,2),col="red",lwd=2,lty=2)
dev.off()


Cd8Krt8Mix_spotFeature <- function(Coordinate, Cd8cell, Krt8cell, AEcell, M){
	XLYL <- c(0,130)
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[Cd8cell] <- "red"
	color_feature[Krt8cell] <- "blue"
	color_feature[intersect(Cd8cell, Krt8cell)] <- "purple"
	color_feature[AEcell] <- "#AAAA00"
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}

Cd8Krt8Mix_spotFeature2col <- function(Coordinate, Cd8cell, Krt8cell, M){
	XLYL <- c(0,130)
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[Cd8cell] <- "red"
	color_feature[Krt8cell] <- "blue"
	color_feature[intersect(Cd8cell, Krt8cell)] <- "purple"
	#color_feature[AEcell] <- "#AAAA00"
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}
Cd8Krt8Mix_spotFeature_single <- function(Coordinate, targetcell, M, usecol){
	XLYL <- c(0,130)
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[targetcell] <- usecol
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}



trimsamplename <- function(inname){
  return(strsplit(inname, "_")[[1]][1])
}
trimcellname <- function(inname){
  return(strsplit(inname, "_")[[1]][2])
}

fetch_sample_cell <- function(cell_list, samplename){
	sample_names <- unlist(lapply(cell_list, trimsamplename))
	cell_names <- unlist(lapply(cell_list, trimcellname))
	return(cell_names[which(sample_names == samplename)])
}
spotFeature_single <- function(Coordinate, targetcell, M, usecol){
	XLYL <- c(0,130)
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[targetcell] <- usecol
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}
ONspotPlot <- function(target_spots, outname, USECOL = "red"){
	target_spots_covid1 <- fetch_sample_cell(target_spots, "covid1")
	target_spots_covid2 <- fetch_sample_cell(target_spots, "covid2")
	target_spots_covid3 <- fetch_sample_cell(target_spots, "covid3")
	target_spots_ctrl1 <-  fetch_sample_cell(target_spots, "ctrl1")
	target_spots_ctrl2 <-  fetch_sample_cell(target_spots, "ctrl2")
	spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, target_spots_covid1, M=paste0(outname, " covid1"), usecol=USECOL)
	spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, target_spots_covid2, M=paste0(outname, " covid2"), usecol=USECOL)
	spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, target_spots_covid3, M=paste0(outname, " covid3"), usecol=USECOL)
	spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  target_spots_ctrl1,  M=paste0(outname, " ctrl1" ), usecol=USECOL)
	spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  target_spots_ctrl2,  M=paste0(outname, " ctrl2" ), usecol=USECOL)
}



spotFeature_plot_score <- function(Coordinate, Feature, log_scale, usecolor, M,zmin,zmax){
	XLYL <- c(0,130)
#	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[,"row"]
	y_pos <- Coordinate[,"col"]
	if(log_scale == 1){
		feature <- log10(Feature+1)
	}else{
		feature <- Feature
	}
	norm01_feature <- norm01fix(feature,zmin,zmax)
	ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
	color_feature <- ColorRamp[norm01_feature]
	data_range <- paste0(zmin," - ",zmax)
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab=data_range, ylab="")
	#return(list(x_pos, y_pos, color_feature))
}

spatialMap_score <- function(scores, usename,zmin,zmax){
	spotFeature_plot_score(covid1_highQ_SCT@images$slice1@coordinates, scores[paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))], usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid1"),zmin=zmin,zmax=zmax)
	spotFeature_plot_score(covid2_highQ_SCT@images$slice1@coordinates, scores[paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))], usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid2"),zmin=zmin,zmax=zmax)
	spotFeature_plot_score(covid3_highQ_SCT@images$slice1@coordinates, scores[paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))], usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " covid3"),zmin=zmin,zmax=zmax)
	spotFeature_plot_score(ctrl1_highQ_SCT@images$slice1@coordinates, scores[paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))], usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl1"),zmin=zmin,zmax=zmax)
	spotFeature_plot_score(ctrl2_highQ_SCT@images$slice1@coordinates, scores[paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))], usecolor = c("blue","white","red"), log_scale = 0, M=paste0(usename, " ctrl2"),zmin=zmin,zmax=zmax)
}


pdf(file=paste0("spotFeature_CTscoreSeparate_white0.pdf"),height=14*5,width=4*5)
par(bg="grey",mfrow=c(14,5),mar=c(4,4,2,2))
spatialMap_score(c(Cd8_score),"Cd8 score",-2,2)
spatialMap_score(c(Krt_score),"Krt score",-2,2)
spatialMap_score(c(Ae_score),"Ae score",-2,2)
spatialMap_score(c(Mdm_score),"Mdm score",-2,2)
spatialMap_score(c(Am_score),"Am score",-2,2)
spatialMap_score(c(CD8A_exp),"CD8A exp",-2,2)
spatialMap_score(c(KRT5_exp),"KRT5 exp",-2,2)
spatialMap_score(c(KRT8_exp),"KRT8 exp",-2,2)
spatialMap_score(c(KRT17_exp),"KRT17 exp",-2,2)
spatialMap_score(c(CD14_exp),"CD814 exp",-2,2)
spatialMap_score(c(CX3CR1_exp),"CX3CR1 exp",-2,2)
spatialMap_score(c(CD68_exp),"CD68 exp",-2,2)
spatialMap_score(c(SFTPC_exp),"AGER exp",-2,2)
dev.off()

pdf(file=paste0("spotFeature_CTscoreSeparate.pdf"),height=14*5,width=4*5)
par(bg="grey",mfrow=c(14,5),mar=c(4,4,2,2))
spatialMap_score(c(Cd8_score),"Cd8 score",-2,4)
spatialMap_score(c(Krt_score),"Krt score",-2,4)
spatialMap_score(c(Ae_score),"Ae score",-2,4)
spatialMap_score(c(Mdm_score),"Mdm score",-2,4)
spatialMap_score(c(Am_score),"Am score",-2,4)
spatialMap_score(c(CD8A_exp),"CD8A exp",-2,4)
spatialMap_score(c(KRT5_exp),"KRT5 exp",-2,4)
spatialMap_score(c(KRT8_exp),"KRT8 exp",-2,4)
spatialMap_score(c(KRT17_exp),"KRT17 exp",-2,4)
spatialMap_score(c(CD14_exp),"CD814 exp",-2,4)
spatialMap_score(c(CX3CR1_exp),"CX3CR1 exp",-2,4)
spatialMap_score(c(CD68_exp),"CD68 exp",-2,4)
spatialMap_score(c(SFTPC_exp),"AGER exp",-2,4)
dev.off()

#### given a cutoff, get ON cell with different genes/genesets
#ONcutoff <- 1.0
for(ONcutoff in seq(0, 2, 0.2)){
	print(ONcutoff)
	Krt_spots <- colnames(scaleMat)[which(Krt_score > ONcutoff)]
	Ae_spots <- colnames(scaleMat)[which(Ae_score > ONcutoff)]
	Cd8_spots <- colnames(scaleMat)[which(Cd8_score > ONcutoff)]
	Mdm_spots <- colnames(scaleMat)[which(Mdm_score > ONcutoff)]
	Am_spots <- colnames(scaleMat)[which(Am_score > ONcutoff)]
	CD8A_spots <- colnames(scaleMat)[which(CD8A_exp > ONcutoff)]
	KRT5_spots <- colnames(scaleMat)[which(KRT5_exp > ONcutoff)]
	KRT8_spots <- colnames(scaleMat)[which(KRT8_exp > ONcutoff)]
	KRT17_spots <- colnames(scaleMat)[which(KRT17_exp > ONcutoff)]
	SFTPC_spots <- colnames(scaleMat)[which(SFTPC_exp > ONcutoff)]
	AGER_spots <- colnames(scaleMat)[which(AGER_exp > ONcutoff)]
	CD14_spots <- colnames(scaleMat)[which(CD14_exp > ONcutoff)]
	CX3CR1_spots <- colnames(scaleMat)[which(CX3CR1_exp > ONcutoff)]
	CD68_spots <- colnames(scaleMat)[which(CD68_exp > ONcutoff)]
	
	pdf(file=paste0("spotFeature_CTscoreOnSeparate_cutoff",ONcutoff,".pdf"),height=14*5,width=4*5)
	par(bg="grey",mfrow=c(14,5),mar=c(4,4,2,2))
	ONspotPlot(Cd8_spots, paste0("cut=",ONcutoff," Cd8"),USECOL="red")
	ONspotPlot(Krt_spots, paste0("cut=",ONcutoff," Krt"),USECOL="blue")
	ONspotPlot(Ae_spots, paste0("cut=",ONcutoff," Ae"),USECOL="#AAAA00")
	ONspotPlot(Mdm_spots, paste0("cut=",ONcutoff," Mdm"),USECOL="green")
	ONspotPlot(Am_spots, paste0("cut=",ONcutoff," Am"),USECOL="brown")
	ONspotPlot(CD8A_spots, paste0("cut=",ONcutoff," CD8A"))
	ONspotPlot(KRT5_spots, paste0("cut=",ONcutoff," KRT5"))
	ONspotPlot(KRT8_spots, paste0("cut=",ONcutoff," KRT8"))
	ONspotPlot(KRT17_spots, paste0("cut=",ONcutoff," KRT17"))
	ONspotPlot(CD14_spots, paste0("cut=",ONcutoff," CD14"))
	ONspotPlot(CX3CR1_spots, paste0("cut=",ONcutoff," CX3CR1"))
	ONspotPlot(CD68_spots, paste0("cut=",ONcutoff," CD68"))
	ONspotPlot(SFTPC_spots, paste0("cut=",ONcutoff," SFTPC"))
	ONspotPlot(AGER_spots, paste0("cut=",ONcutoff," AGER"))
	dev.off()	
}




covid_ctrl_groupBox <- function(scores, usename){
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	boxplot(scores[covid1spots],
			scores[covid2spots],
			scores[covid3spots],
			scores[ctrl1spots],
			scores[ctrl2spots],
			names=c("covid1","covid2","covid3","ctrl1","ctrl2"),outline=F,
			col=c("red","red","red","blue","blue"),main=usename,ylab=usename)
	boxplot(scores[covidspots],
			scores[ctrlspots],
			names=c("covid","ctrl"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[covidspots],scores[ctrlspots],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
}

pdf(file=paste0("scoreBox_covidVSctrl.pdf"),height=20*4,width=4*2)
par(mfrow=c(20,2),mar=c(4,4,2,2))
covid_ctrl_groupBox(Cd8_score,"Cd8 score")
covid_ctrl_groupBox(Krt_score,"Krt score")
covid_ctrl_groupBox(Ae_score,"Ae score")
covid_ctrl_groupBox(Mdm_score,"Mdm score")
covid_ctrl_groupBox(Am_score,"Am score")
covid_ctrl_groupBox(CD8A_exp,"CD8A exp")
covid_ctrl_groupBox(KRT5_exp,"KRT5 exp")
covid_ctrl_groupBox(KRT8_exp,"KRT8 exp")
covid_ctrl_groupBox(KRT17_exp,"KRT17 exp")
covid_ctrl_groupBox(CD14_exp,"CD14 exp")
covid_ctrl_groupBox(CX3CR1_exp,"CX3CR1 exp")
covid_ctrl_groupBox(CD68_exp,"CD68 exp")
covid_ctrl_groupBox(SFTPC_exp,"SFTPC exp")
covid_ctrl_groupBox(AGER_exp,"AGER exp")
covid_ctrl_groupBox(TGFb_score,"TGFb_score")
covid_ctrl_groupBox(IL1R_score,"IL1R_score")
covid_ctrl_groupBox(Inflammasome_score,"Inflammasome_score")
covid_ctrl_groupBox(IFNG_score,"IFNG_score")
covid_ctrl_groupBox(TNF_score,"TNF_score")
covid_ctrl_groupBox(IFNGTNF_score,"IFNGTNF_score")
dev.off()




percent_ON <- function(values, ONcutoff){
	return(length(which(values > ONcutoff))/length(values) )
}

covid_ctrl_groupONpercent <- function(scores, usename, ONcutoff){
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	PV <- t.test(c(percent_ON(scores[covid1spots], ONcutoff),
			       percent_ON(scores[covid2spots], ONcutoff),
			       percent_ON(scores[covid3spots], ONcutoff)),
	             c(percent_ON(scores[ctrl1spots], ONcutoff),
			       percent_ON(scores[ctrl2spots], ONcutoff)))$p.val/2
	barplot(c(percent_ON(scores[covid1spots], ONcutoff),
			percent_ON(scores[covid2spots], ONcutoff),
			percent_ON(scores[covid3spots], ONcutoff),
			percent_ON(scores[ctrl1spots] , ONcutoff),
			percent_ON(scores[ctrl2spots] , ONcutoff)),
			names=c("covid1","covid2","covid3","ctrl1","ctrl2"),
			col=c("red","red","red","blue","blue"),main=usename,ylab="ONspots%",xlab=PV)
	barplot(c(percent_ON(scores[covidspots], ONcutoff),
			percent_ON(scores[ctrlspots], ONcutoff)),
			names=c("covid","ctrl"),col=c("red","blue"),main=usename,ylab="ONspots%")

	return(c(usename,ONcutoff,
		    percent_ON(scores[covidspots], ONcutoff),
			percent_ON(scores[ctrlspots], ONcutoff),
			percent_ON(scores[covid1spots], ONcutoff),
			percent_ON(scores[covid2spots], ONcutoff),
			percent_ON(scores[covid3spots], ONcutoff),
			percent_ON(scores[ctrl1spots] , ONcutoff),
			percent_ON(scores[ctrl2spots] , ONcutoff),
			PV
			))
}

covidCtrlSummary <- c()
for(ONcutoff in seq(0,2,0.2)){
	pdf(file=paste0("scoreONpercentBar_covidVSctrl_ONcutoff",ONcutoff,".pdf"),height=20*4,width=4*2)
	par(mfrow=c(20,2),mar=c(4,4,2,2))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(CD8A_exp,"CD8A exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(KRT5_exp,"KRT5 exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(KRT8_exp,"KRT8 exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(KRT17_exp,"KRT17 exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(CD14_exp,"CD14 exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(CX3CR1_exp,"CX3CR1 exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(CD68_exp,"CD68 exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(SFTPC_exp,"SFTPC exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(AGER_exp,"AGER exp",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(Cd8_score,"Cd8 score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(Krt_score,"Krt score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(Ae_score,"Ae score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(Mdm_score,"Mdm score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(Am_score,"Am score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(TGFb_score,"TGFb_score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(IL1R_score,"IL1R_score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(Inflammasome_score,"Inflammasome_score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(IFNG_score,"IFNG_score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(TNF_score,"TNF_score",ONcutoff))
	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(IFNGTNF_score,"IFNGTNF_score",ONcutoff))
	dev.off()
}
#rownames(covidCtrlSummary) <- c("CD8A","KRT5","KRt8","KRT17","CD14","CX3CR1","CD68","SFTPC","AGER",
#	          
#                  "CD8score","KRTscore","AEscore","MDMscore","AMscore","TGFb","IL1R","Inflammasome","IFNG","TNF","IFNGTNF")
colnames(covidCtrlSummary) <- c("term","ONcutoff","covid","ctrl","covid1","covid2","covid3","ctrl1","ctrl2","pval")

write.table(covidCtrlSummary, file="covidCtrl_ONpercent_Summary.txt",row.names=F,col.names=T,sep="\t",quote=F)



percent_ON <- function(values, ONcutoff){
	return(length(which(values > ONcutoff))/length(values) )
}

covid_ctrl_groupONpercent_CB_singleGene <- function( gname, ONcutoff){
	scores  <- scaleMat[toupper(gname),]
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	PV <- t.test(c(percent_ON(scores[covid1spots], ONcutoff),
			       percent_ON(scores[covid2spots], ONcutoff),
			       percent_ON(scores[covid3spots], ONcutoff)),
	             c(percent_ON(scores[ctrl1spots], ONcutoff),
			       percent_ON(scores[ctrl2spots], ONcutoff)))$p.val/2
	barplot(c(percent_ON(scores[covid1spots], ONcutoff),
			percent_ON(scores[covid2spots], ONcutoff),
			percent_ON(scores[covid3spots], ONcutoff),
			percent_ON(scores[ctrl1spots] , ONcutoff),
			percent_ON(scores[ctrl2spots] , ONcutoff)),
			names=c("covid1","covid2","covid3","ctrl1","ctrl2"),
			col=c("red","red","red","blue","blue"),main=gname,ylab="ONspots%",xlab=PV)
	barplot(c(percent_ON(scores[covidspots], ONcutoff),
			percent_ON(scores[ctrlspots], ONcutoff)),
			names=c("covid","ctrl"),col=c("red","blue"),main=gname,ylab="ONspots%",xlab="")
			
	return(c(percent_ON(scores[covidspots], ONcutoff),
			percent_ON(scores[ctrlspots], ONcutoff),
			percent_ON(scores[covid1spots], ONcutoff),
			percent_ON(scores[covid2spots], ONcutoff),
			percent_ON(scores[covid3spots], ONcutoff),
			percent_ON(scores[ctrl1spots] , ONcutoff),
			percent_ON(scores[ctrl2spots] , ONcutoff),
			PV
			))
}


CTkey_glist <- c("SFTPC","CD8A","KRT5","KRT8","KRT17","SFTPC","AGER","CD14","CX3CR1","CD68")

pdf(file="CTkeyGene_covidctrlONpercent.pdf",width=6,height=30)
par(mfrow=c(10,2),mar=c(4,4,2,2))
CTkey_covidctrlONpercent <- c()
for(G in c(CTkey_glist)){
	CTkey_covidctrlONpercent <- rbind(CTkey_covidctrlONpercent, 
		                              covid_ctrl_groupONpercent_CB_singleGene(G, 1.2))
}
dev.off()
rownames(CTkey_covidctrlONpercent) <- CTkey_glist
colnames(CTkey_covidctrlONpercent) <- c("covid","ctrl","covid1","covid2","covid3","ctrl1","ctrl2","pval")

write.table(CTkey_covidctrlONpercent, file="CTkey_covidctrlONpercent",row.names=T,col.names=T,sep="\t",quote=F)


pdf(file="newAddGene_covidctrlONpercent.pdf",width=6,height=33)
par(mfrow=c(10,2),mar=c(4,4,2,2))
newAddGene_covidctrlONpercent <- c()
for(G in c(summary3_glist)){
	newAddGene_covidctrlONpercent <- rbind(newAddGene_covidctrlONpercent, 
		                              covid_ctrl_groupONpercent_CB_singleGene(G, 1.2))
}
dev.off()
rownames(newAddGene_covidctrlONpercent) <- summary3_glist
colnames(newAddGene_covidctrlONpercent) <- c("covid","ctrl","covid1","covid2","covid3","ctrl1","ctrl2","pval")

write.table(newAddGene_covidctrlONpercent, file="newAddGene_covidctrlONpercent.txt",row.names=T,col.names=T,sep="\t",quote=F)


ONcutoff=1.2
for(gname in c(CTkey_glist,summary3_glist)){
	scores  <- scaleMat[toupper(gname),]
	target_spots <- colnames(scaleMat)[which(scores > ONcutoff)]
	pdf(file=paste0("spotFeatures_keygenes/keygeneOnSeparate/",gname,"_cutoff",ONcutoff,".pdf"),width=18,height=9)
	par(mfrow=c(2,3),mar=c(4,4,2,2),bg="grey")
	ONspotPlot(target_spots, paste0(gname," cut=",ONcutoff),USECOL="red")
	dev.off()
}


#SFTPC_spots <- colnames(scaleMat)[which(SFTPC_exp > ONcutoff)]
#	ONspotPlot(Am_spots, paste0("cut=",ONcutoff," Am"),USECOL="brown")
#
#summary3_glist
#pdf(file=paste0("scoreONpercentBar_covidVSctrl_ONcutoff",ONcutoff,".pdf"),height=17*4,width=4*2)
#par(mfrow=c(17,2),mar=c(4,4,2,2))
#covid_ctrl_groupONpercent(Cd8_score,"Cd8 score",ONcutoff)
#covid_ctrl_groupONpercent(Krt_score,"Krt score",ONcutoff)
#	covid_ctrl_groupONpercent(CD68_exp,"CD68 exp",ONcutoff)








#g3 <- rownames(scaleMat)[which(scaleMat[,"KRT17"] >= 2 & scaleMat[,"AGER"]  >= 2 )]

KRT17_AGER_groupBox <- function(scores, usename){	
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	KRT17high_spots <- colnames(scaleMat)[which(scaleMat["KRT17",] >= 2 & scaleMat["AGER",] < 2 )]
	AGERhigh_spots <- colnames(scaleMat)[which(scaleMat["KRT17",] < 2 & scaleMat["AGER",]  >= 2 )]
	KRT17high_covid <- intersect(KRT17high_spots, covidspots)
	AGERhigh_covid <- intersect(AGERhigh_spots, covidspots)
	KRT17high_covid1 <- intersect(KRT17high_spots, covid1spots)
	AGERhigh_covid1 <- intersect(AGERhigh_spots, covid1spots)
	KRT17high_covid2 <- intersect(KRT17high_spots, covid2spots)
	AGERhigh_covid2 <- intersect(AGERhigh_spots, covid2spots)
	KRT17high_covid3 <- intersect(KRT17high_spots, covid3spots)
	AGERhigh_covid3 <- intersect(AGERhigh_spots, covid3spots)

	boxplot(scores[KRT17high_spots],
			scores[AGERhigh_spots],
			names=c("KRT17high","AGERhigh"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT17high_spots],scores[AGERhigh_spots],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT17high_covid],
			scores[AGERhigh_covid],
			names=c("KRT17highCovid","AGERhighCovid"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT17high_covid],scores[AGERhigh_covid],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT17high_covid1],
			scores[AGERhigh_covid1],
			names=c("KRT17highCovid1","AGERhighCovid1"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT17high_covid1],scores[AGERhigh_covid1],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT17high_covid2],
			scores[AGERhigh_covid2],
			names=c("KRT17highCovid2","AGERhighCovid2"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT17high_covid2],scores[AGERhigh_covid2],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT17high_covid3],
			scores[AGERhigh_covid3],
			names=c("KRT17highCovid3","AGERhighCovid3"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT17high_covid3],scores[AGERhigh_covid3],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
		#return(c(scores[KRT17high_spots],scores[AGERhigh_spots],scores[KRT17high_covid],scores[AGERhigh_covid]))
}

#KRT17vsAGER_summary <- c()
pdf(file=paste0("scoreBox_KRT17vsAGER.pdf"),height=20*4,width=4*5)
par(mfrow=c(20,5),mar=c(4,4,2,2))
KRT17_AGER_groupBox(Cd8_score,"Cd8 score")
KRT17_AGER_groupBox(Krt_score,"Krt score")
KRT17_AGER_groupBox(Ae_score,"Ae score")
KRT17_AGER_groupBox(Mdm_score,"Mdm score")
KRT17_AGER_groupBox(Am_score,"Am score")
KRT17_AGER_groupBox(CD8A_exp,"CD8A exp")
KRT17_AGER_groupBox(KRT5_exp,"KRT5 exp")
KRT17_AGER_groupBox(KRT8_exp,"KRT8 exp")
KRT17_AGER_groupBox(KRT17_exp,"KRT17 exp")
KRT17_AGER_groupBox(CD14_exp,"CD14 exp")
KRT17_AGER_groupBox(CX3CR1_exp,"CX3CR1 exp")
KRT17_AGER_groupBox(CD68_exp,"CD68 exp")
KRT17_AGER_groupBox(SFTPC_exp,"SFTPC exp")
KRT17_AGER_groupBox(AGER_exp,"AGER exp")
KRT17_AGER_groupBox(TGFb_score,"TGFb_score")
KRT17_AGER_groupBox(IL1R_score,"IL1R_score")
KRT17_AGER_groupBox(Inflammasome_score,"Inflammasome_score")
KRT17_AGER_groupBox(IFNG_score,"IFNG_score")
KRT17_AGER_groupBox(TNF_score,"TNF_score")
KRT17_AGER_groupBox(IFNGTNF_score,"IFNGTNF_score")
dev.off()



KRT17_AGER_groupONpercent <- function(scores, usename, ONcutoff){
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	KRT17high_spots <- colnames(scaleMat)[which(scaleMat["KRT17",] >= 2 & scaleMat["AGER",] < 2 )]
	AGERhigh_spots <- colnames(scaleMat)[which(scaleMat["KRT17",] < 2 & scaleMat["AGER",]  >= 2 )]
	KRT17high_covid <- intersect(KRT17high_spots, covidspots)
	AGERhigh_covid <- intersect(AGERhigh_spots, covidspots)
	KRT17high_covid1 <- intersect(KRT17high_spots, covid1spots)
	AGERhigh_covid1 <- intersect(AGERhigh_spots, covid1spots)
	KRT17high_covid2 <- intersect(KRT17high_spots, covid2spots)
	AGERhigh_covid2 <- intersect(AGERhigh_spots, covid2spots)
	KRT17high_covid3 <- intersect(KRT17high_spots, covid3spots)
	AGERhigh_covid3 <- intersect(AGERhigh_spots, covid3spots)
	barplot(c(percent_ON(scores[KRT17high_spots], ONcutoff),
			percent_ON(scores[AGERhigh_spots], ONcutoff)),
			names=c("KRT17high_spots","AGERhigh_spots"),
			col=c("red","blue"),main=usename,ylab="ONspots%")
	barplot(c(percent_ON(scores[KRT17high_covid], ONcutoff),
			percent_ON(scores[AGERhigh_covid], ONcutoff)),
			names=c("KRT17high_covid","AGERhigh_covid"),col=c("red","blue"),main=usename,ylab="ONspots%")
#	barplot(c(percent_ON(scores[KRT17high_covid1], ONcutoff),
#			percent_ON(scores[AGERhigh_covid1], ONcutoff)),
#			names=c("KRT17high_covid1","AGERhigh_covid1"),col=c("red","blue"),main=usename,ylab="ONspots%")
#	barplot(c(percent_ON(scores[KRT17high_covid2], ONcutoff),
#			percent_ON(scores[AGERhigh_covid2], ONcutoff)),
#			names=c("KRT17high_covid2","AGERhigh_covid2"),col=c("red","blue"),main=usename,ylab="ONspots%")
#	barplot(c(percent_ON(scores[KRT17high_covid3], ONcutoff),
#			percent_ON(scores[AGERhigh_covid3], ONcutoff)),
#			names=c("KRT17high_covid3","AGERhigh_covid3"),col=c("red","blue"),main=usename,ylab="ONspots%")
	return(c(usename,ONcutoff,
		    percent_ON(scores[KRT17high_spots], ONcutoff),
			percent_ON(scores[AGERhigh_spots], ONcutoff),
			percent_ON(scores[KRT17high_covid], ONcutoff),
			percent_ON(scores[AGERhigh_covid], ONcutoff)
			))
}
#	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(CD8A_exp,"CD8A exp",ONcutoff))

KRT17summary <- c()

for(ONcutoff in seq(0,2,0.2)){
	pdf(file=paste0("scoreONpercentBar_KRT17vsAGER_ONcutoff",ONcutoff,".pdf"),height=20*4,width=4*2)
	par(mfrow=c(20,2),mar=c(4,4,2,2))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(Cd8_score,"Cd8 score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(Krt_score,"Krt score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(Ae_score,"Ae score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(Mdm_score,"Mdm score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(Am_score,"Am score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(CD8A_exp,"CD8A exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(KRT5_exp,"KRT5 exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(KRT8_exp,"KRT8 exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(KRT17_exp,"KRT17 exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(CD14_exp,"CD14 exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(CX3CR1_exp,"CX3CR1 exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(CD68_exp,"CD68 exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(SFTPC_exp,"SFTPC exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(AGER_exp,"AGER exp",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(TGFb_score,"TGFb_score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(IL1R_score,"IL1R_score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(Inflammasome_score,"Inflammasome_score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(IFNG_score,"IFNG_score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(TNF_score,"TNF_score",ONcutoff))
	KRT17summary <- rbind(KRT17summary, KRT17_AGER_groupONpercent(IFNGTNF_score,"IFNGTNF_score",ONcutoff))
	dev.off()
}
colnames(KRT17summary) <- c("term","ONcutoff","KRT17high_all","AGERhigh_all","KRT17high_covid","AGERhigh_covid")
write.table(KRT17summary, file="KRT17vsAGER_ONpercent_Summary.txt",row.names=F,col.names=T,sep="\t",quote=F)




#g3 <- rownames(scaleMat)[which(scaleMat[,"KRT8"] >= 2 & scaleMat[,"AGER"]  >= 2 )]

KRT8_AGER_groupBox <- function(scores, usename){	
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	KRT8high_spots <- colnames(scaleMat)[which(scaleMat["KRT8",] >= 2 & scaleMat["AGER",] < 2 )]
	AGERhigh_spots <- colnames(scaleMat)[which(scaleMat["KRT8",] < 2 & scaleMat["AGER",]  >= 2 )]
	KRT8high_covid <- intersect(KRT8high_spots, covidspots)
	AGERhigh_covid <- intersect(AGERhigh_spots, covidspots)
	KRT8high_covid1 <- intersect(KRT8high_spots, covid1spots)
	AGERhigh_covid1 <- intersect(AGERhigh_spots, covid1spots)
	KRT8high_covid2 <- intersect(KRT8high_spots, covid2spots)
	AGERhigh_covid2 <- intersect(AGERhigh_spots, covid2spots)
	KRT8high_covid3 <- intersect(KRT8high_spots, covid3spots)
	AGERhigh_covid3 <- intersect(AGERhigh_spots, covid3spots)

	boxplot(scores[KRT8high_spots],
			scores[AGERhigh_spots],
			names=c("KRT8high","AGERhigh"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT8high_spots],scores[AGERhigh_spots],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT8high_covid],
			scores[AGERhigh_covid],
			names=c("KRT8highCovid","AGERhighCovid"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT8high_covid],scores[AGERhigh_covid],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT8high_covid1],
			scores[AGERhigh_covid1],
			names=c("KRT8highCovid1","AGERhighCovid1"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT8high_covid1],scores[AGERhigh_covid1],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT8high_covid2],
			scores[AGERhigh_covid2],
			names=c("KRT8highCovid2","AGERhighCovid2"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT8high_covid2],scores[AGERhigh_covid2],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
	boxplot(scores[KRT8high_covid3],
			scores[AGERhigh_covid3],
			names=c("KRT8highCovid3","AGERhighCovid3"),outline=F,
			col=c("red","blue"),main=usename,ylab=usename)
		PV <- wilcox.test(scores[KRT8high_covid3],scores[AGERhigh_covid3],,alternative="two.sided")$p.val
		legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
		#return(c(scores[KRT8high_spots],scores[AGERhigh_spots],scores[KRT8high_covid],scores[AGERhigh_covid]))
}

#KRT8vsAGER_summary <- c()
pdf(file=paste0("scoreBox_KRT8vsAGER.pdf"),height=20*4,width=4*5)
par(mfrow=c(20,5),mar=c(4,4,2,2))
KRT8_AGER_groupBox(Cd8_score,"Cd8 score")
KRT8_AGER_groupBox(Krt_score,"Krt score")
KRT8_AGER_groupBox(Ae_score,"Ae score")
KRT8_AGER_groupBox(Mdm_score,"Mdm score")
KRT8_AGER_groupBox(Am_score,"Am score")
KRT8_AGER_groupBox(CD8A_exp,"CD8A exp")
KRT8_AGER_groupBox(KRT5_exp,"KRT5 exp")
KRT8_AGER_groupBox(KRT8_exp,"KRT8 exp")
KRT17_AGER_groupBox(KRT17_exp,"KRT17 exp")
KRT8_AGER_groupBox(CD14_exp,"CD14 exp")
KRT8_AGER_groupBox(CX3CR1_exp,"CX3CR1 exp")
KRT8_AGER_groupBox(CD68_exp,"CD68 exp")
KRT8_AGER_groupBox(SFTPC_exp,"SFTPC exp")
KRT8_AGER_groupBox(AGER_exp,"AGER exp")
KRT8_AGER_groupBox(TGFb_score,"TGFb_score")
KRT8_AGER_groupBox(IL1R_score,"IL1R_score")
KRT8_AGER_groupBox(Inflammasome_score,"Inflammasome_score")
KRT8_AGER_groupBox(IFNG_score,"IFNG_score")
KRT8_AGER_groupBox(TNF_score,"TNF_score")
KRT8_AGER_groupBox(IFNGTNF_score,"IFNGTNF_score")
dev.off()



KRT8_AGER_groupONpercent <- function(scores, usename, ONcutoff){
	covid1spots <- paste0("covid1_",rownames(covid1_highQ_SCT@images$slice1@coordinates))
	covid2spots <- paste0("covid2_",rownames(covid2_highQ_SCT@images$slice1@coordinates))
	covid3spots <- paste0("covid3_",rownames(covid3_highQ_SCT@images$slice1@coordinates))
	ctrl1spots <- paste0("ctrl1_",rownames(ctrl1_highQ_SCT@images$slice1@coordinates))
	ctrl2spots <- paste0("ctrl2_",rownames(ctrl2_highQ_SCT@images$slice1@coordinates))
	covidspots <- c(covid1spots,covid2spots,covid3spots )
	ctrlspots <- c(ctrl1spots,ctrl2spots )
	KRT8high_spots <- colnames(scaleMat)[which(scaleMat["KRT8",] >= 2 & scaleMat["AGER",] < 2 )]
	AGERhigh_spots <- colnames(scaleMat)[which(scaleMat["KRT8",] < 2 & scaleMat["AGER",]  >= 2 )]
	KRT8high_covid <- intersect(KRT8high_spots, covidspots)
	AGERhigh_covid <- intersect(AGERhigh_spots, covidspots)
	KRT8high_covid1 <- intersect(KRT8high_spots, covid1spots)
	AGERhigh_covid1 <- intersect(AGERhigh_spots, covid1spots)
	KRT8high_covid2 <- intersect(KRT8high_spots, covid2spots)
	AGERhigh_covid2 <- intersect(AGERhigh_spots, covid2spots)
	KRT8high_covid3 <- intersect(KRT8high_spots, covid3spots)
	AGERhigh_covid3 <- intersect(AGERhigh_spots, covid3spots)
	barplot(c(percent_ON(scores[KRT8high_spots], ONcutoff),
			percent_ON(scores[AGERhigh_spots], ONcutoff)),
			names=c("KRT8high_spots","AGERhigh_spots"),
			col=c("red","blue"),main=usename,ylab="ONspots%")
	barplot(c(percent_ON(scores[KRT8high_covid], ONcutoff),
			percent_ON(scores[AGERhigh_covid], ONcutoff)),
			names=c("KRT8high_covid","AGERhigh_covid"),col=c("red","blue"),main=usename,ylab="ONspots%")
#	barplot(c(percent_ON(scores[KRT8high_covid1], ONcutoff),
#			percent_ON(scores[AGERhigh_covid1], ONcutoff)),
#			names=c("KRT8high_covid1","AGERhigh_covid1"),col=c("red","blue"),main=usename,ylab="ONspots%")
#	barplot(c(percent_ON(scores[KRT8high_covid2], ONcutoff),
#			percent_ON(scores[AGERhigh_covid2], ONcutoff)),
#			names=c("KRT8high_covid2","AGERhigh_covid2"),col=c("red","blue"),main=usename,ylab="ONspots%")
#	barplot(c(percent_ON(scores[KRT8high_covid3], ONcutoff),
#			percent_ON(scores[AGERhigh_covid3], ONcutoff)),
#			names=c("KRT8high_covid3","AGERhigh_covid3"),col=c("red","blue"),main=usename,ylab="ONspots%")
	return(c(usename,ONcutoff,
		    percent_ON(scores[KRT8high_spots], ONcutoff),
			percent_ON(scores[AGERhigh_spots], ONcutoff),
			percent_ON(scores[KRT8high_covid], ONcutoff),
			percent_ON(scores[AGERhigh_covid], ONcutoff)
			))
}
#	covidCtrlSummary <- rbind(covidCtrlSummary,covid_ctrl_groupONpercent(CD8A_exp,"CD8A exp",ONcutoff))

KRT8summary <- c()

for(ONcutoff in seq(0,2,0.2)){
	pdf(file=paste0("scoreONpercentBar_KRT8vsAGER_ONcutoff",ONcutoff,".pdf"),height=20*4,width=4*2)
	par(mfrow=c(20,2),mar=c(4,4,2,2))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(Cd8_score,"Cd8 score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(Krt_score,"Krt score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(Ae_score,"Ae score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(Mdm_score,"Mdm score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(Am_score,"Am score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(CD8A_exp,"CD8A exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(KRT5_exp,"KRT5 exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(KRT8_exp,"KRT8 exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(KRT17_exp,"KRT17 exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(CD14_exp,"CD14 exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(CX3CR1_exp,"CX3CR1 exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(CD68_exp,"CD68 exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(SFTPC_exp,"SFTPC exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(AGER_exp,"AGER exp",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(TGFb_score,"TGFb_score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(IL1R_score,"IL1R_score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(Inflammasome_score,"Inflammasome_score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(IFNG_score,"IFNG_score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(TNF_score,"TNF_score",ONcutoff))
	KRT8summary <- rbind(KRT8summary, KRT8_AGER_groupONpercent(IFNGTNF_score,"IFNGTNF_score",ONcutoff))
	dev.off()
}
colnames(KRT8summary) <- c("term","ONcutoff","KRT8high_all","AGERhigh_all","KRT8high_covid","AGERhigh_covid")
write.table(KRT8summary, file="KRT8vsAGER_ONpercent_Summary.txt",row.names=F,col.names=T,sep="\t",quote=F)








spotFeature_score <- function(glist){
	gscore <- apply(scaleMat[intersect(toupper(glist),rownames(scaleMat)),],2,mean)
	return(gscore)
}

outdata <- cbind(
spotFeature_score(CD8_refine2),
spotFeature_score(Krt_refine2),
spotFeature_score(MDM_refine2),
spotFeature_score(AE_refine2),
#spotFeature_score(c(AT1_Hglist,AT2_Hglist)),
spotFeature_score(AM_refine2))

colnames(outdata) <- c("CD8","Krt","MDM","AE","AM")
write.table(outdata,file="spotLevel_CTscore_newHg.txt",row.names=T,col.names=T,sep="\t",quote=F)


#cor(t(scaleMat[usegenelist,]))
#outdata <- cbind(
#scaleMat["KRT8",,
#spotFeature_score("KRT17"),
#spotFeature_score("SFTPC"),
#spotFeature_score("AGER"))
#spotFeature_score(c(AT1_Hglist,AT2_Hglist)),
#outdata <- as.data.frame(as.matrix(t(scaleMat[c(usegenelist,TGFb_pathway,IL1R_pathway,Inflammasome_pathway),])))

usegenelist <- c("AGER","KRT8","SFTPC","KRT17","CD8A","CD14","KRT5","CX3CR1","CD68")
TGFb_pathway <- intersect(rownames(scaleMat),c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP"))
IL1R_pathway <- intersect(rownames(scaleMat),c("IL1A","IL1B","IL1RN","TNF","IL6","TGFB1","IL1R1","IRAK1","IRAK1","TGFB2","IRAK3","RELA","NFKB1","IL1RAP","IL1RAP","IL1RAP","IL1RAP","MYD88","MYD88","MYD88","MYD88","IKBKB","IKBKB","RELA","RELA","CHUK","MAPK14","IKBKB","IRAK1","IRAK2","IFNB1","IL1RAP","JUN","MYD88","MAPK8","MAP2K3","MAP2K6","MAP3K7","TGFB2","TGFB3","MAP3K14","NFKB1","TRAF6","MAP3K1","TAB1","IRAK3","NFKBIA","RELA","IFNA1","IL1RAP","MAPK14","MAPK14","MAPK14","MAPK8"))
Inflammasome_pathway <- intersect(rownames(scaleMat),c("NFKB2","P2RX7","NLRC4","NLRP1","HSP90AB1","HMOX1","MEFV","PYCARD","NFKB1","PANX1","TXN","CASP1","PSTPIP1","APP","NLRP3","AIM2","SUGT1","BCL2L1","BCL2","RELA","TXNIP"))
IFNG_pathway <- intersect(rownames(scaleMat),read.table("pathwayInfo/IFNG.txt")[,1])
TNF_pathway <-  intersect(rownames(scaleMat),read.table("pathwayInfo/TNF.txt")[,1])

outdata <- as.data.frame(as.matrix(t(scaleMat[c(usegenelist),])))

covid1_outdata <- outdata[grep("covid1",rownames(outdata)),]
covid2_outdata <- outdata[grep("covid2",rownames(outdata)),]
covid3_outdata <- outdata[grep("covid3",rownames(outdata)),]
ctrl1_outdata <- outdata[grep("ctrl1",rownames(outdata)),]
ctrl2_outdata <- outdata[grep("ctrl2",rownames(outdata)),]
#ctrl3_outdata <- outdata[grep("ctrl3",rownames(outdata)),]
covid_outdata <-  rbind(covid1_outdata,covid2_outdata,covid3_outdata)
ctrl_outdata <-  rbind(ctrl1_outdata,ctrl2_outdata)


#covid1_outdata2 <- covid1_outdata[intersect(rownames(covid1_outdata), Krt_spots),]
#covid2_outdata2<- covid2_outdata[intersect(rownames(covid2_outdata), Krt_spots),]
#covid3_outdata2<- covid3_outdata[intersect(rownames(covid3_outdata), Krt_spots),]
#ctrl1_outdata2<- ctrl1_outdata[intersect(rownames(ctrl1_outdata), Krt_spots),]
#ctrl2_outdata2<- ctrl2_outdata[intersect(rownames(ctrl2_outdata), Krt_spots),]
##ctrl3_outdata
#covid_outdata2 <- covid_outdata[intersect(rownames(covid_outdata), Krt_spots),]
#ctrl_outdata2 <- ctrl_outdata[intersect(rownames(ctrl_outdata), Krt_spots),]

CTscore_cor_all <- cor(outdata,method="spearman")
CTscore_cor_covid1 <- cor(covid1_outdata,method="spearman")
CTscore_cor_covid2 <- cor(covid2_outdata,method="spearman")
CTscore_cor_covid3 <- cor(covid3_outdata,method="spearman")
CTscore_cor_ctrl1 <- cor(ctrl1_outdata,method="spearman")
CTscore_cor_ctrl2 <- cor(ctrl2_outdata,method="spearman")
#CTscore_cor_ctrl3 <- cor(ctrl3_outdata,method="spearman")
CTscore_cor_covid <- cor(covid_outdata,method="spearman")
CTscore_cor_ctrl <- cor(ctrl_outdata,method="spearman")

CTscore_cor_all_sub <- cor(outdata[,c("AGER","CD8A","KRT5","KRT8","KRT17","SFTPC")],method="spearman")
CTscore_cor_covid_sub <- cor(covid_outdata[,c("AGER","CD8A","KRT5","KRT8","KRT17","SFTPC")],method="spearman")




#CTscore2_cor_all <- cor(outdata2,method="spearman")
#CTscore2_cor_covid1 <- cor(covid1_outdata2,method="spearman")
#CTscore2_cor_covid2 <- cor(covid2_outdata2,method="spearman")
#CTscore2_cor_covid3 <- cor(covid3_outdata2,method="spearman")
#CTscore2_cor_ctrl1 <- cor(ctrl1_outdata2,method="spearman")
#CTscore2_cor_ctrl2 <- cor(ctrl2_outdata2,method="spearman")
#CTscore2_cor_covid <- cor(covid_outdata2,method="spearman")
#CTscore2_cor_ctrl <- cor(ctrl_outdata2,method="spearman")


my_palette <- colorRampPalette(c("blue", "white", "red"))

pdf(file="keygenes_corHeat_allspots_cellnote.pdf")
cn <- matrix(as.character(round(CTscore_cor_all,2)),ncol=dim(CTscore_cor_all)[2])
heatmap.2(CTscore_cor_all,col=my_palette(100),trace="none",main="all spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
dev.off()

pdf(file="keygenesSUB_corHeat_allspots_cellnote.pdf")
cn <- matrix(as.character(round(CTscore_cor_all_sub,2)),ncol=dim(CTscore_cor_all_sub)[2])
heatmap.2(CTscore_cor_all_sub,col=my_palette(100),trace="none",main="all spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
dev.off()

pdf(file="keygenes_corHeat_covidspots_cellnote.pdf")
cn <- matrix(as.character(round(CTscore_cor_covid,2)),ncol=dim(CTscore_cor_covid)[2])
heatmap.2(CTscore_cor_covid,col=my_palette(100),trace="none",main="covid spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
dev.off()

pdf(file="keygenesSUB_corHeat_covidspots_cellnote.pdf")
cn <- matrix(as.character(round(CTscore_cor_covid_sub,2)),ncol=dim(CTscore_cor_all_sub)[2])
heatmap.2(CTscore_cor_covid_sub,col=my_palette(100),trace="none",main="covid spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
dev.off()


#### 
percentOV <- function(dataA, dataB, ONcutoff){
	ON_A <- length(which(dataA > ONcutoff))
	ON_B <- length(which(dataB > ONcutoff))
	length(which(dataA > ONcutoff & dataB > ONcutoff))
	fisher_A <- length(which(dataA > ONcutoff & dataB > ONcutoff))
	fisher_B <- length(which(dataA <= ONcutoff & dataB > ONcutoff))
	fisher_C <- length(which(dataA <= ONcutoff & dataB <= ONcutoff))
	fisher_D <- length(which(dataA > ONcutoff & dataB <= ONcutoff))
	fisher_out <- fisher.test(matrix(c(fisher_A,fisher_B,fisher_D,fisher_C),nrow=2),alternative="greater")
	return(-log10(fisher_out$p.value))
}
ONcutoff = 1




n <- 9
my_palette <- colorRampPalette(c("blue", "white", "red"))


for(ONcutoff in seq(0,2,0.2)){

	result <- matrix(0, n, n)
	for (i in 1:n) {
	  for (j in 1:n) {
	    result[i, j] <- percentOV(outdata[,i], outdata[,j],ONcutoff)
	  }
	}
	
	result_covid <- matrix(0, n, n)
	for (i in 1:n) {
	  for (j in 1:n) {
	    result_covid[i, j] <- percentOV(covid_outdata[,i], covid_outdata[,j],ONcutoff)
	  }
	}

	result[result>50] <- 50
	result_covid[result_covid>50] <- 50
	
	rownames(result) <- colnames(outdata)
	colnames(result) <- colnames(outdata)
	rownames(result_covid) <- colnames(outdata)
	colnames(result_covid) <- colnames(outdata)
	
	pdf(file=paste0("keygenes_OVpercentFisherP_allspots_cutoff",ONcutoff,".pdf"))
	cn <- matrix(as.character(round(result,1)),ncol=dim(result)[2])
	heatmap.2(result,col=my_palette(100),trace="none",main="all spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
	dev.off()
	
	pdf(file=paste0("keygenes_OVpercentFisherP_covidspots_cutoff",ONcutoff,".pdf"))
	cn <- matrix(as.character(round(result_covid,1)),ncol=dim(result_covid)[2])
	heatmap.2(result_covid,col=my_palette(100),trace="none",main="covid spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
	dev.off()
}


covid_outdata

n <- ncol(outdata)
result <- matrix(0, n, n)







#barplot(c(CTscore_cor_covid1["AE","Krt"],
#CTscore_cor_covid2["AE","Krt"],
#CTscore_cor_covid3["AE","Krt"],
#CTscore_cor_ctrl1["AE","Krt"],
#CTscore_cor_ctrl2["AE","Krt"]),ylab="R AEvsKrt",names=c("covid1","covid2","covid3","ctrl1","ctrl2"),col=c("red","red","red","blue","blue"))
feature2color <- function(Feature){
	norm01_feature <- norm01(Feature)
	ColorRamp <- colorRampPalette(c("white","red"), bias=1)(101)   #color list
	color_feature <- ColorRamp[norm01_feature]
	data_range <- paste0(round(min(Feature),1)," - ",round(max(Feature),1))
	return(color_feature)
}
#Feature <- outdata[,"CD8A"]
#plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab=data_range, ylab="")

pdf(file="AGER_vs_KRT17_expScatter_colorFeature.pdf",width=12,height=12)
par(bg="grey",mfrow=c(2,2),mar=c(4,4,2,2))
plot(outdata[,"AGER"],outdata[,"KRT17"],
	 col=feature2color(outdata[,"CD8A"]),pch=".",main="allsamples",xlab="AGER",ylab="KRT17")
legend("topright",legend="CD8A",text.col="red",bty="n",cex=2)
plot(outdata[,"AGER"],outdata[,"KRT17"],
	 col=feature2color(apply(outdata[,TGFb_pathway],1,mean)),pch=".",main="allsamples",xlab="AGER",ylab="KRT17")
legend("topright",legend="TGFb",text.col="red",bty="n",cex=2)
plot(outdata[,"AGER"],outdata[,"KRT17"],
	 col=feature2color(apply(outdata[,IL1R_pathway],1,mean)),pch=".",main="allsamples",xlab="AGER",ylab="KRT17")
legend("topright",legend="IL1R",text.col="red",bty="n",cex=2)
plot(outdata[,"AGER"],outdata[,"KRT17"],
	 col=feature2color(apply(outdata[,Inflammasome_pathway],1,mean)),pch=".",main="allsamples",xlab="AGER",ylab="KRT17")
legend("topright",legend="Inflammasome",text.col="red",bty="n",cex=2)
dev.off()



######
outdata$TGFb <- as.numeric(apply(outdata[,TGFb_pathway],1,mean))
outdata$IL1R <- apply(outdata[,IL1R_pathway],1,mean)
outdata$Inflammasome <- apply(outdata[,Inflammasome_pathway],1,mean)

g1 <- rownames(outdata)[grep("covid",rownames(outdata))]
g2 <- rownames(outdata)[grep("ctrl",rownames(outdata))]


pdf(file="covid_vs_ctrl_geneExpBoxplot.pdf",width=16,height=12)
par(bg="grey",mfrow=c(3,4),mar=c(4,4,2,2))
for(Feature in c("CD8A","KRT5","KRT8","KRT17","SFTPC","AGER","CD14","CX3CR1","CD68","TGFb","IL1R","Inflammasome")){
	boxplot(outdata[g1,Feature],outdata[g2,Feature],
			names=c("covid","ctrl"),ylab=paste0(Feature, " normExp"),main=Feature,outline=F,col=c("red","blue"))
	PV <- wilcox.test(outdata[g1,Feature],outdata[g2,Feature],alternative="less")$p.val
	legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
}
dev.off()


g1 <- rownames(outdata)[which(outdata[,"KRT17"] >= 2 & outdata[,"AGER"] < 2 )]
g2 <- rownames(outdata)[which(outdata[,"KRT17"] < 2 & outdata[,"AGER"]  >= 2 )]
g3 <- rownames(outdata)[which(outdata[,"KRT17"] >= 2 & outdata[,"AGER"]  >= 2 )]

pdf(file="AGER_vs_KRT17_expScatter_cutoff.pdf",width=18,height=6)
par(mfrow=c(1,3),mar=c(4,4,2,2))
plot(outdata[,"AGER"],outdata[,"KRT17"],pch=".",main="allsamples",xlab="AGER",ylab="KRT17")
abline(h=2,v=2,col="red")
smoothScatter(outdata[,"AGER"],outdata[,"KRT17"],main="allsamples",xlab="AGER",ylab="KRT17")
abline(h=2,v=2,col="red")
plot(outdata[,"AGER"],outdata[,"KRT17"],pch=".",main="allsamples",xlab="AGER",ylab="KRT17")
points(outdata[g1,"AGER"],outdata[g1,"KRT17"],pch=".",col="red")
points(outdata[g2,"AGER"],outdata[g2,"KRT17"],pch=".",col="blue")
points(outdata[g3,"AGER"],outdata[g3,"KRT17"],pch=".",col="#AAAA00")
dev.off()


pdf(file="KRT17_vs_AGER_uniqGroup_geneExpBoxplot.pdf",width=16,height=12)
par(bg="grey",mfrow=c(3,4),mar=c(4,4,2,2))
for(Feature in c("CD8A","KRT5","KRT8","KRT17","SFTPC","AGER","CD14","CX3CR1","CD68","TGFb","IL1R","Inflammasome")){
	boxplot(outdata[g1,Feature],outdata[g2,Feature],outdata[g3,Feature],
			names=c("KRT17uniq","AGERuniq","both"),ylab=paste0(Feature, " normExp"),main=Feature,outline=F,col=c("red","blue","#AAAA00"))
	PV <- wilcox.test(outdata[g1,Feature],outdata[g2,Feature],alternative="less")$p.val
	legend("topleft",legend=paste0("P=",PV),bty="n",cex=1.2)
}
dev.off()


outdata <- t(scaleMat)
g1 <- rownames(outdata)[which(outdata[,"KRT8"] >= 4 & outdata[,"AGER"] < 4 )]
g2 <- rownames(outdata)[which(outdata[,"KRT8"] < 4 & outdata[,"AGER"]  >= 4 )]
g3 <- rownames(outdata)[which(outdata[,"KRT8"] >= 4 & outdata[,"AGER"]  >= 4 )]

pdf(file="AGER_vs_KRT8_expScatter_cutoff.pdf",width=18,height=6)
par(mfrow=c(1,3),mar=c(4,4,2,2))
plot(outdata[,"AGER"],outdata[,"KRT8"],pch=".",main="allsamples",xlab="AGER",ylab="KRT8")
abline(h=4,v=4,col="red")
smoothScatter(outdata[,"AGER"],outdata[,"KRT8"],main="allsamples",xlab="AGER",ylab="KRT8")
abline(h=4,v=4,col="red")
plot(outdata[,"AGER"],outdata[,"KRT8"],pch=".",main="allsamples",xlab="AGER",ylab="KRT8")
points(outdata[g1,"AGER"],outdata[g1,"KRT8"],pch=".",col="red")
points(outdata[g2,"AGER"],outdata[g2,"KRT8"],pch=".",col="blue")
points(outdata[g3,"AGER"],outdata[g3,"KRT8"],pch=".",col="#AAAA00")
dev.off()

############### DEG
#covid_ctrl <- rep("NA",length(cbdata$sample))
#covid_ctrl[which(cbdata$sample %in% c("covid1","covid2","covid3"))] <- "covid"
cbdata$disease[which(cbdata$sample == "ctrl3")] <- "NA"
cbdata$KRT17AGER <- rep("NA",length(cbdata$sample))
cbdata$KRT17AGER[which(names(cbdata$sample) %in% g1)] <- "KRT17high"
cbdata$KRT17AGER[which(names(cbdata$sample) %in% g2)] <- "AGERhigh"
cbdata$KRT17AGER[which(names(cbdata$sample) %in% g3)] <- "bothhigh"

cbdata$KRT17AGERcovid <- rep("NA",length(cbdata$sample))
cbdata$KRT17AGERcovid[which(names(cbdata$sample) %in% g1 & cbdata$disease == "covid")] <- "KRT17high"
cbdata$KRT17AGERcovid[which(names(cbdata$sample) %in% g2 & cbdata$disease == "covid")] <- "AGERhigh"
cbdata$KRT17AGERcovid[which(names(cbdata$sample) %in% g3 & cbdata$disease == "covid")] <- "bothhigh"


DEG_covid_ctrl <- FindMarkers(cbdata,group.by="disease", ident.1 = "covid", ident.2 = "ctrl")
DEG_KRT17_AGER <- FindMarkers(cbdata,group.by="KRT17AGER", ident.1 = "KRT17high", ident.2 = "AGERhigh")
DEG_KRT17_AGER_covid <- FindMarkers(cbdata,group.by="KRT17AGERcovid", ident.1 = "KRT17high", ident.2 = "AGERhigh")

DEG_covid_g_ctrl <- DEG_covid_ctrl[which(DEG_covid_ctrl[,"p_val_adj"] < 0.001 & (DEG_covid_ctrl[,"avg_log2FC"]) >= log2(1.5)),]
DEG_covid_l_ctrl <- DEG_covid_ctrl[which(DEG_covid_ctrl[,"p_val_adj"] < 0.001 & (DEG_covid_ctrl[,"avg_log2FC"]) <= -log2(1.5)),]
DEG_KRT17_g_AGER <- DEG_KRT17_AGER[which(DEG_KRT17_AGER[,"p_val_adj"] < 0.001 & (DEG_KRT17_AGER[,"avg_log2FC"]) >= log2(1.5)),]
DEG_KRT17_l_AGER <- DEG_KRT17_AGER[which(DEG_KRT17_AGER[,"p_val_adj"] < 0.001 & (DEG_KRT17_AGER[,"avg_log2FC"]) <= -log2(1.5)),]
DEG_KRT17_g_AGER_covid <- DEG_KRT17_AGER_covid[which(DEG_KRT17_AGER_covid[,"p_val_adj"] < 0.001 & (DEG_KRT17_AGER_covid[,"avg_log2FC"]) >= log2(1.5)),]
DEG_KRT17_l_AGER_covid <- DEG_KRT17_AGER_covid[which(DEG_KRT17_AGER_covid[,"p_val_adj"] < 0.001 & (DEG_KRT17_AGER_covid[,"avg_log2FC"]) <= -log2(1.5)),]


write.table(row.names(DEG_covid_g_ctrl),file="DEG_covid_g_ctrl_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_covid_l_ctrl),file="DEG_covid_l_ctrl_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT17_g_AGER),file="DEG_KRT17_g_AGER_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT17_l_AGER),file="DEG_KRT17_l_AGER_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT17_g_AGER_covid),file="DEG_KRT17_g_AGER_covid_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT17_l_AGER_covid),file="DEG_KRT17_l_AGER_covid_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )

write.table((DEG_covid_g_ctrl),file="DEG_covid_g_ctrl_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_covid_l_ctrl),file="DEG_covid_l_ctrl_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT17_g_AGER),file="DEG_KRT17_g_AGER_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT17_l_AGER),file="DEG_KRT17_l_AGER_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT17_g_AGER_covid),file="DEG_KRT17_g_AGER_covid_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT17_l_AGER_covid),file="DEG_KRT17_l_AGER_covid_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )


DEG_covid_ctrl_ranklist <- c(DEG_covid_g_ctrl[,2],DEG_covid_l_ctrl[,2])
names(DEG_covid_ctrl_ranklist) <- c(row.names(DEG_covid_g_ctrl),row.names(DEG_covid_l_ctrl))
DEG_KRT17_AGER_ranklist <- c(DEG_KRT17_g_AGER[,2],DEG_KRT17_l_AGER[,2])
names(DEG_KRT17_AGER_ranklist) <- c(row.names(DEG_KRT17_g_AGER),row.names(DEG_KRT17_l_AGER))
DEG_KRT17_AGER_covid_ranklist <- c(DEG_KRT17_g_AGER_covid[,2],DEG_KRT17_l_AGER_covid[,2])
names(DEG_KRT17_AGER_covid_ranklist) <- c(row.names(DEG_KRT17_g_AGER_covid),row.names(DEG_KRT17_l_AGER_covid))

write.table(DEG_covid_ctrl_ranklist,file="DEG_covid_ctrl_DEG.rnk",row.names=T,col.names=F,sep="\t",quote=F )
write.table(DEG_KRT17_AGER_ranklist,file="DEG_KRT17_AGER_DEG.rnk",row.names=T,col.names=F,sep="\t",quote=F )
write.table(DEG_KRT17_AGER_covid_ranklist,file="DEG_KRT17_AGER_covid_DEG.rnk",row.names=T,col.names=F,sep="\t",quote=F )


### KRT8 ver
cbdata$disease[which(cbdata$sample == "ctrl3")] <- "NA"
cbdata$KRT8AGER <- rep("NA",length(cbdata$sample))
cbdata$KRT8AGER[which(names(cbdata$sample) %in% g1)] <- "KRT8high"
cbdata$KRT8AGER[which(names(cbdata$sample) %in% g2)] <- "AGERhigh"
cbdata$KRT8AGER[which(names(cbdata$sample) %in% g3)] <- "bothhigh"

cbdata$KRT8AGERcovid <- rep("NA",length(cbdata$sample))
cbdata$KRT8AGERcovid[which(names(cbdata$sample) %in% g1 & cbdata$disease == "covid")] <- "KRT8high"
cbdata$KRT8AGERcovid[which(names(cbdata$sample) %in% g2 & cbdata$disease == "covid")] <- "AGERhigh"
cbdata$KRT8AGERcovid[which(names(cbdata$sample) %in% g3 & cbdata$disease == "covid")] <- "bothhigh"

DEG_KRT8_AGER <- FindMarkers(cbdata,group.by="KRT8AGER", ident.1 = "KRT8high", ident.2 = "AGERhigh")
DEG_KRT8_AGER_covid <- FindMarkers(cbdata,group.by="KRT8AGERcovid", ident.1 = "KRT8high", ident.2 = "AGERhigh")

DEG_KRT8_g_AGER <- DEG_KRT8_AGER[which( (DEG_KRT8_AGER[,"avg_log2FC"]) > 0),]
DEG_KRT8_l_AGER <- DEG_KRT8_AGER[which( (DEG_KRT8_AGER[,"avg_log2FC"]) < 0),]
DEG_KRT8_g_AGER_covid <- DEG_KRT8_AGER_covid[which((DEG_KRT8_AGER_covid[,"avg_log2FC"]) > 0),]
DEG_KRT8_l_AGER_covid <- DEG_KRT8_AGER_covid[which((DEG_KRT8_AGER_covid[,"avg_log2FC"]) <0 ),]

write.table(row.names(DEG_KRT8_g_AGER),file="DEG_KRT8_g_AGER_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT8_l_AGER),file="DEG_KRT8_l_AGER_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT8_g_AGER_covid),file="DEG_KRT8_g_AGER_covid_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT8_l_AGER_covid),file="DEG_KRT8_l_AGER_covid_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )


DEG_KRT8_g_AGER <- DEG_KRT8_AGER[which(DEG_KRT8_AGER[,"p_val_adj"] < 0.001 & (DEG_KRT8_AGER[,"avg_log2FC"]) >= log2(1.5)),]
DEG_KRT8_l_AGER <- DEG_KRT8_AGER[which(DEG_KRT8_AGER[,"p_val_adj"] < 0.001 & (DEG_KRT8_AGER[,"avg_log2FC"]) <= -log2(1.5)),]
DEG_KRT8_g_AGER_covid <- DEG_KRT8_AGER_covid[which(DEG_KRT8_AGER_covid[,"p_val_adj"] < 0.001 & (DEG_KRT8_AGER_covid[,"avg_log2FC"]) >= log2(1.5)),]
DEG_KRT8_l_AGER_covid <- DEG_KRT8_AGER_covid[which(DEG_KRT8_AGER_covid[,"p_val_adj"] < 0.001 & (DEG_KRT8_AGER_covid[,"avg_log2FC"]) <= -log2(1.5)),]

write.table(row.names(DEG_KRT8_g_AGER),file="DEG_KRT8_g_AGER_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT8_l_AGER),file="DEG_KRT8_l_AGER_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT8_g_AGER_covid),file="DEG_KRT8_g_AGER_covid_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(DEG_KRT8_l_AGER_covid),file="DEG_KRT8_l_AGER_covid_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )

write.table((DEG_KRT8_g_AGER),file="DEG_KRT8_g_AGER_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT8_l_AGER),file="DEG_KRT8_l_AGER_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT8_g_AGER_covid),file="DEG_KRT8_g_AGER_covid_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )
write.table((DEG_KRT8_l_AGER_covid),file="DEG_KRT8_l_AGER_covid_highQ_summary.txt",row.names=T,col.names=T,sep="\t",quote=F )

DEG_KRT8_AGER_ranklist <- c(DEG_KRT8_g_AGER[,2],DEG_KRT8_l_AGER[,2])
names(DEG_KRT8_AGER_ranklist) <- c(row.names(DEG_KRT8_g_AGER),row.names(DEG_KRT8_l_AGER))
DEG_KRT8_AGER_covid_ranklist <- c(DEG_KRT8_g_AGER_covid[,2],DEG_KRT8_l_AGER_covid[,2])
names(DEG_KRT8_AGER_covid_ranklist) <- c(row.names(DEG_KRT8_g_AGER_covid),row.names(DEG_KRT8_l_AGER_covid))

write.table(DEG_KRT8_AGER_ranklist,file="DEG_KRT8_AGER_DEG.rnk",row.names=T,col.names=F,sep="\t",quote=F )
write.table(DEG_KRT8_AGER_covid_ranklist,file="DEG_KRT8_AGER_covid_DEG.rnk",row.names=T,col.names=F,sep="\t",quote=F )









##### cell type score cor
#useCT <- c("CD8","MDM","Krt","AE","AM")

idx <- names(Krt_score)
outdata <- cbind(Krt_score[idx],
                 Ae_score[idx],
                 Mdm_score[idx],
                 Am_score[idx]
                  )
rownames(outdata) <- idx
colnames(outdata) <- c("Krt","Ae","Mdm","Am")

covid1_outdata <- outdata[grep("covid1",rownames(outdata)),]
covid2_outdata <- outdata[grep("covid2",rownames(outdata)),]
covid3_outdata <- outdata[grep("covid3",rownames(outdata)),]
ctrl1_outdata <- outdata[grep("ctrl1",rownames(outdata)),]
ctrl2_outdata <- outdata[grep("ctrl2",rownames(outdata)),]
#ctrl3_outdata <- outdata[grep("ctrl3",rownames(outdata)),]
covid_outdata <-  rbind(covid1_outdata,covid2_outdata,covid3_outdata)
ctrl_outdata <-  rbind(ctrl1_outdata,ctrl2_outdata)


CTscore_cor_all <- cor(outdata,method="spearman")
CTscore_cor_covid <- cor(covid_outdata,method="spearman")

my_palette <- colorRampPalette(c( "white", "red"))
pdf(file="CTscore_corHeat_allspots_cellnote.pdf")
cn <- matrix(as.character(round(CTscore_cor_all,2)),ncol=dim(CTscore_cor_all)[2])
heatmap.2(CTscore_cor_all,col=my_palette(100),trace="none",main="all spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
dev.off()

pdf(file="CTscore_corHeat_covidspots_cellnote.pdf")
cn <- matrix(as.character(round(CTscore_cor_covid,2)),ncol=dim(CTscore_cor_covid)[2])
heatmap.2(CTscore_cor_covid,col=my_palette(100),trace="none",main="covid spots",margins=c(8,8),cellnote=cn, notecol="black",cexRow=0.8,cexCol=0.8)
dev.off()

#pdf(file="CTscore_5CTcorHeatmap_newHg_all.pdf")
#cn <- matrix(as.character(round(CTscore_cor_all,2)),ncol=dim(CTscore_cor_all)[2])
#heatmap.2(CTscore_cor_all,col=my_palette(100),trace="none",main="all",margins=c(8,8))
#dev.off()
#pdf(file="CTscore_5CTcorHeatmap_newHg_covid.pdf")
#heatmap.2(CTscore_cor_covid,col=my_palette(100),trace="none",main="covid",margins=c(8,8))
#dev.off()

#pdf(file="CTscore_5CTcorHeatmap_newHg_ctrl.pdf")
#heatmap.2(CTscore_cor_ctrl,col=my_palette(100),trace="none",main="ctrl",margins=c(8,8))
#dev.off()
#


#pdf(file="CTscore_5CTcorHeatmap_newHg_ctrl1.pdf")
#heatmap.2(CTscore_cor_ctrl1,col=my_palette(100),trace="none",main="ctrl1",margins=c(8,8))
#dev.off()
#
#pdf(file="CTscore_5CTcorHeatmap_newHg_ctrl2.pdf")
#heatmap.2(CTscore_cor_ctrl2,col=my_palette(100),trace="none",main="ctrl2",margins=c(8,8))
#dev.off()
#
#pdf(file="CTscore_5CTcorHeatmap_newHg_covid1.pdf")
#heatmap.2(CTscore_cor_covid1,col=my_palette(100),trace="none",main="covid1",margins=c(8,8))
#dev.off()
#
#pdf(file="CTscore_5CTcorHeatmap_newHg_covid2.pdf")
#heatmap.2(CTscore_cor_covid2,col=my_palette(100),trace="none",main="covid2",margins=c(8,8))
#dev.off()
#
#pdf(file="CTscore_5CTcorHeatmap_newHg_covid3.pdf")
#heatmap.2(CTscore_cor_covid3,col=my_palette(100),trace="none",main="covid3",margins=c(8,8))
#dev.off()





pdf(file="CTscore_5CTcorHeatmap_newHg_ctrl1_v2.pdf")
heatmap.2(CTscore_cor_ctrl1,col=my_palette(100),trace="none",main="ctrl1",margins=c(8,8))
dev.off()

pdf(file="CTscore_5CTcorHeatmap_newHg_ctrl2.pdf")
heatmap.2(CTscore_cor_ctrl2,col=my_palette(100),trace="none",main="ctrl2",margins=c(8,8))
dev.off()

pdf(file="CTscore_5CTcorHeatmap_newHg_covid1.pdf")
heatmap.2(CTscore_cor_covid1,col=my_palette(100),trace="none",main="covid1",margins=c(8,8))
dev.off()

pdf(file="CTscore_5CTcorHeatmap_newHg_covid2.pdf")
heatmap.2(CTscore_cor_covid2,col=my_palette(100),trace="none",main="covid2",margins=c(8,8))
dev.off()

pdf(file="CTscore_5CTcorHeatmap_newHg_covid3.pdf")
heatmap.2(CTscore_cor_covid3,col=my_palette(100),trace="none",main="covid3",margins=c(8,8))
dev.off()

















Krt_spots <- colnames(scaleMat)[which(Krt_score > 0)]
Ae_spots <- colnames(scaleMat)[which(Ae_score > 0)]
Cd8_spots <- colnames(scaleMat)[which(Cd8_score > 0)]
Mdm_spots <- colnames(scaleMat)[which(Mdm_score > 0)]
Am_spots <- colnames(scaleMat)[which(Am_score > 0)]


#Krt_spots_use <- setdiff(Krt_spots, Ae_spots)#colnames(scaleMat)[which(KrtExp > -0.4)]

### AE vs Krt
covid3spot <- colnames(scaleMat)[grep("covid3",colnames(scaleMat))]
covid3_Ae <- intersect(covid3spot, Ae_spots)
covid3_Krt <- intersect(covid3spot, Krt_spots)
Cd8_score <- scaleMat["CD8A",]
names(Cd8_score) <- colnames(scaleMat)
pdf(file="AEvsKrt_diffscoreEnrich.pdf",width=15,height=4)
par(mfrow=c(1,5),mar=c(6,4,2,2))
boxplot(Cd8_score[covid3_Ae],
	    Cd8_score[covid3_Krt],
	    Cd8_score[setdiff(covid3_Ae,covid3_Krt)],
	    Cd8_score[setdiff(covid3_Krt,covid3_Ae)],las=2,
	    names=c("Ae","Krt","Ae uniq","Krt uniq"),main="covid3",ylab="cd8 score")
boxplot(Mdm_score[covid3_Ae],
	    Mdm_score[covid3_Krt],
	    Mdm_score[setdiff(covid3_Ae,covid3_Krt)],
	    Mdm_score[setdiff(covid3_Krt,covid3_Ae)],las=2,
	    names=c("Ae","Krt","Ae uniq","Krt uniq"),main="covid3",ylab="Mdm score")
boxplot(Am_score[covid3_Ae],
	    Am_score[covid3_Krt],
	    Am_score[setdiff(covid3_Ae,covid3_Krt)],
	    Am_score[setdiff(covid3_Krt,covid3_Ae)],las=2,
	    names=c("Ae","Krt","Ae uniq","Krt uniq"),main="covid3",ylab="Am score")
boxplot(Ae_score[covid3_Ae],
	    Ae_score[covid3_Krt],
	    Ae_score[setdiff(covid3_Ae,covid3_Krt)],
	    Ae_score[setdiff(covid3_Krt,covid3_Ae)],las=2,
	    names=c("Ae","Krt","Ae uniq","Krt uniq"),main="covid3",ylab="Ae score")
boxplot(Krt_score[covid3_Ae],
	    Krt_score[covid3_Krt],
	    Krt_score[setdiff(covid3_Ae,covid3_Krt)],
	    Krt_score[setdiff(covid3_Krt,covid3_Ae)],las=2,
	    names=c("Ae","Krt","Ae uniq","Krt uniq"),main="covid3",ylab="Krt score")
dev.off()



Krt_spots_covid1 <- fetch_sample_cell(Krt_spots, "covid1")
Cd8_spots_covid1 <- fetch_sample_cell(Cd8_spots, "covid1")
Ae_spots_covid1 <- fetch_sample_cell(Ae_spots, "covid1")
Mdm_spots_covid1 <- fetch_sample_cell(Mdm_spots, "covid1")
Am_spots_covid1 <- fetch_sample_cell(Am_spots, "covid1")

Krt_spots_covid2 <- fetch_sample_cell(Krt_spots, "covid2")
Cd8_spots_covid2 <- fetch_sample_cell(Cd8_spots, "covid2")
Ae_spots_covid2 <- fetch_sample_cell(Ae_spots, "covid2")
Mdm_spots_covid2 <- fetch_sample_cell(Mdm_spots, "covid2")
Am_spots_covid2 <- fetch_sample_cell(Am_spots, "covid2")

Krt_spots_covid3 <- fetch_sample_cell(Krt_spots, "covid3")
Cd8_spots_covid3 <- fetch_sample_cell(Cd8_spots, "covid3")
Ae_spots_covid3 <- fetch_sample_cell(Ae_spots, "covid3")
Mdm_spots_covid3 <- fetch_sample_cell(Mdm_spots, "covid3")
Am_spots_covid3 <- fetch_sample_cell(Am_spots, "covid3")

Krt_spots_ctrl1 <- fetch_sample_cell(Krt_spots, "ctrl1")
Cd8_spots_ctrl1 <- fetch_sample_cell(Cd8_spots, "ctrl1")
Ae_spots_ctrl1 <- fetch_sample_cell(Ae_spots, "ctrl1")
Mdm_spots_ctrl1 <- fetch_sample_cell(Mdm_spots, "ctrl1")
Am_spots_ctrl1 <- fetch_sample_cell(Am_spots, "ctrl1")

Krt_spots_ctrl2 <- fetch_sample_cell(Krt_spots, "ctrl2")
Cd8_spots_ctrl2 <- fetch_sample_cell(Cd8_spots, "ctrl2")
Ae_spots_ctrl2 <- fetch_sample_cell(Ae_spots, "ctrl2")
Mdm_spots_ctrl2 <- fetch_sample_cell(Mdm_spots, "ctrl2")
Am_spots_ctrl2 <- fetch_sample_cell(Am_spots, "ctrl2")

#Krt_spots_ctrl3 <- fetch_sample_cell(Krt_spots, "ctrl3")
#Cd8_spots_ctrl3 <- fetch_sample_cell(Cd8_spots, "ctrl3")
#Ae_spots_ctrl3 <- fetch_sample_cell(Ae_spots, "ctrl3")
#Mdm_spots_ctrl3 <- fetch_sample_cell(Mdm_spots, "ctrl3")
#At_spots_ctrl3 <- fetch_sample_cell(At_spots, "ctrl3")



pdf(file=paste0("spotFeature_CTscoreSeparate.pdf"),height=20,width=20)
par(bg="grey",mfrow=c(5,5),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots_covid1, M="cd8 covid1", usecol="red")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots_covid2, M="cd8 covid2", usecol="red")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots_covid3, M="cd8 covid3", usecol="red")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Cd8_spots_ctrl1,  M="cd8 ctrl1" , usecol="red")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Cd8_spots_ctrl2,  M="cd8 ctrl2" , usecol="red")
#Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  Cd8_spots_ctrl3,  M="cd8 ctrl3" , usecol="red")

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Krt_spots_covid1, M="Krt covid1", usecol="blue")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Krt_spots_covid2, M="Krt covid2", usecol="blue")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Krt_spots_covid3, M="Krt covid3", usecol="blue")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Krt_spots_ctrl1,  M="Krt ctrl1" , usecol="blue")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Krt_spots_ctrl2,  M="Krt ctrl2" , usecol="blue")
#Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  Krt_spots_ctrl3,  M="Krt ctrl3" , usecol="blue")

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Ae_spots_covid1, M="Ae covid1", usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Ae_spots_covid2, M="Ae covid2", usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Ae_spots_covid3, M="Ae covid3", usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Ae_spots_ctrl1,  M="Ae ctrl1" , usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Ae_spots_ctrl2,  M="Ae ctrl2" , usecol="#AAAA00")
#Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  Ae_spots_ctrl3,  M="Ae ctrl3" , usecol="#AAAA00")

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Mdm_spots_covid1, M="Mdm covid1", usecol="green")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Mdm_spots_covid2, M="Mdm covid2", usecol="green")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Mdm_spots_covid3, M="Mdm covid3", usecol="green")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Mdm_spots_ctrl1,  M="Mdm ctrl1" , usecol="green")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Mdm_spots_ctrl2,  M="Mdm ctrl2" , usecol="green")
#Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  Mdm_spots_ctrl3,  M="Mdm ctrl3" , usecol="green")

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Am_spots_covid1, M="Am covid1", usecol="brown")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Am_spots_covid2, M="Am covid2", usecol="brown")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Am_spots_covid3, M="Am covid3", usecol="brown")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Am_spots_ctrl1,  M="Am ctrl1" , usecol="brown")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Am_spots_ctrl2,  M="Am ctrl2" , usecol="brown")
#Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  At_spots_ctrl3,  M="At ctrl3" , usecol="brown")

dev.off()





pdf(file=paste0("spotFeature_CTscore_Cd8KrtMix.pdf"),height=16,width=24)
par(bg="grey",mfrow=c(4,6),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature2col(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots_covid1, Krt_spots_covid1, M="cd8Krt covid1")
Cd8Krt8Mix_spotFeature2col(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots_covid2, Krt_spots_covid2, M="cd8Krt covid2")
Cd8Krt8Mix_spotFeature2col(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots_covid3, Krt_spots_covid3, M="cd8Krt covid3")
Cd8Krt8Mix_spotFeature2col(ctrl1_highQ_SCT@images$slice1@coordinates,  Cd8_spots_ctrl1,  Krt_spots_ctrl1,  M="cd8Krt ctrl1" )
Cd8Krt8Mix_spotFeature2col(ctrl2_highQ_SCT@images$slice1@coordinates,  Cd8_spots_ctrl2,  Krt_spots_ctrl2,  M="cd8Krt ctrl2" )
Cd8Krt8Mix_spotFeature2col(ctrl3_highQ_SCT@images$slice1@coordinates,  Cd8_spots_ctrl3,  Krt_spots_ctrl3,  M="cd8Krt ctrl3" )

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Ae_spots_covid1, M="Ae covid1", usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Ae_spots_covid2, M="Ae covid2", usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Ae_spots_covid3, M="Ae covid3", usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Ae_spots_ctrl1,  M="Ae ctrl1" , usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Ae_spots_ctrl2,  M="Ae ctrl2" , usecol="#AAAA00")
Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  Ae_spots_ctrl3,  M="Ae ctrl3" , usecol="#AAAA00")

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, Mdm_spots_covid1, M="Mdm covid1", usecol="green")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, Mdm_spots_covid2, M="Mdm covid2", usecol="green")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, Mdm_spots_covid3, M="Mdm covid3", usecol="green")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  Mdm_spots_ctrl1,  M="Mdm ctrl1" , usecol="green")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  Mdm_spots_ctrl2,  M="Mdm ctrl2" , usecol="green")
Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  Mdm_spots_ctrl3,  M="Mdm ctrl3" , usecol="green")

Cd8Krt8Mix_spotFeature_single(covid1_highQ_SCT@images$slice1@coordinates, At_spots_covid1, M="At covid1", usecol="brown")
Cd8Krt8Mix_spotFeature_single(covid2_highQ_SCT@images$slice1@coordinates, At_spots_covid2, M="At covid2", usecol="brown")
Cd8Krt8Mix_spotFeature_single(covid3_highQ_SCT@images$slice1@coordinates, At_spots_covid3, M="At covid3", usecol="brown")
Cd8Krt8Mix_spotFeature_single(ctrl1_highQ_SCT@images$slice1@coordinates,  At_spots_ctrl1,  M="At ctrl1" , usecol="brown")
Cd8Krt8Mix_spotFeature_single(ctrl2_highQ_SCT@images$slice1@coordinates,  At_spots_ctrl2,  M="At ctrl2" , usecol="brown")
Cd8Krt8Mix_spotFeature_single(ctrl3_highQ_SCT@images$slice1@coordinates,  At_spots_ctrl3,  M="At ctrl3" , usecol="brown")

dev.off()












Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_covid1, Krt_spots_use_covid1, Ae_spots_covid1, M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_covid2, Krt_spots_use_covid2, Ae_spots_covid2, M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_covid3, Krt_spots_use_covid3, Ae_spots_covid3, M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_ctrl1, Krt_spots_use_ctrl1, Ae_spots_ctrl1, M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_ctrl2, Krt_spots_use_ctrl2, Ae_spots_ctrl2, M="ctrl2")









Krt_spots_use_covid1 <- fetch_sample_cell(Krt_spots_use, "covid1")
Ae_spots_covid1 <- fetch_sample_cell(Ae_spots, "covid1")
Cd8_spots1_use_covid1 <- fetch_sample_cell(Cd8_spots1_use, "covid1")
Cd8_spots2_use_covid1 <- fetch_sample_cell(Cd8_spots2_use, "covid1")
Cd8_spots3_use_covid1 <- fetch_sample_cell(Cd8_spots3_use, "covid1")
Cd8_spots4_use_covid1 <- fetch_sample_cell(Cd8_spots4_use, "covid1")
Cd8_spots5_use_covid1 <- fetch_sample_cell(Cd8_spots5_use, "covid1")

Krt_spots_use_covid2 <- fetch_sample_cell(Krt_spots_use, "covid2")
Ae_spots_covid2 <- fetch_sample_cell(Ae_spots, "covid2")
Cd8_spots1_use_covid2 <- fetch_sample_cell(Cd8_spots1_use, "covid2")
Cd8_spots2_use_covid2 <- fetch_sample_cell(Cd8_spots2_use, "covid2")
Cd8_spots3_use_covid2 <- fetch_sample_cell(Cd8_spots3_use, "covid2")
Cd8_spots4_use_covid2 <- fetch_sample_cell(Cd8_spots4_use, "covid2")
Cd8_spots5_use_covid2 <- fetch_sample_cell(Cd8_spots5_use, "covid2")

Krt_spots_use_covid3 <- fetch_sample_cell(Krt_spots_use, "covid3")
Ae_spots_covid3 <- fetch_sample_cell(Ae_spots, "covid3")
Cd8_spots1_use_covid3 <- fetch_sample_cell(Cd8_spots1_use, "covid3")
Cd8_spots2_use_covid3 <- fetch_sample_cell(Cd8_spots2_use, "covid3")
Cd8_spots3_use_covid3 <- fetch_sample_cell(Cd8_spots3_use, "covid3")
Cd8_spots4_use_covid3 <- fetch_sample_cell(Cd8_spots4_use, "covid3")
Cd8_spots5_use_covid3 <- fetch_sample_cell(Cd8_spots5_use, "covid3")

Krt_spots_use_ctrl1 <- fetch_sample_cell(Krt_spots_use, "ctrl1")
Ae_spots_ctrl1 <- fetch_sample_cell(Ae_spots, "ctrl1")
Cd8_spots1_use_ctrl1 <- fetch_sample_cell(Cd8_spots1_use, "ctrl1")
Cd8_spots2_use_ctrl1 <- fetch_sample_cell(Cd8_spots2_use, "ctrl1")
Cd8_spots3_use_ctrl1 <- fetch_sample_cell(Cd8_spots3_use, "ctrl1")
Cd8_spots4_use_ctrl1 <- fetch_sample_cell(Cd8_spots4_use, "ctrl1")
Cd8_spots5_use_ctrl1 <- fetch_sample_cell(Cd8_spots5_use, "ctrl1")

Krt_spots_use_ctrl2 <- fetch_sample_cell(Krt_spots_use, "ctrl2")
Ae_spots_ctrl2 <- fetch_sample_cell(Ae_spots, "ctrl2")
Cd8_spots1_use_ctrl2 <- fetch_sample_cell(Cd8_spots1_use, "ctrl2")
Cd8_spots2_use_ctrl2 <- fetch_sample_cell(Cd8_spots2_use, "ctrl2")
Cd8_spots3_use_ctrl2 <- fetch_sample_cell(Cd8_spots3_use, "ctrl2")
Cd8_spots4_use_ctrl2 <- fetch_sample_cell(Cd8_spots4_use, "ctrl2")
Cd8_spots5_use_ctrl2 <- fetch_sample_cell(Cd8_spots5_use, "ctrl2")

Krt_spots_use_ctrl3 <- fetch_sample_cell(Krt_spots_use, "ctrl3")
Ae_spots_ctrl3 <- fetch_sample_cell(Ae_spots, "ctrl3")
Cd8_spots1_use_ctrl3 <- fetch_sample_cell(Cd8_spots1_use, "ctrl3")
Cd8_spots2_use_ctrl3 <- fetch_sample_cell(Cd8_spots2_use, "ctrl3")
Cd8_spots3_use_ctrl3 <- fetch_sample_cell(Cd8_spots3_use, "ctrl3")
Cd8_spots4_use_ctrl3 <- fetch_sample_cell(Cd8_spots4_use, "ctrl3")
Cd8_spots5_use_ctrl3 <- fetch_sample_cell(Cd8_spots5_use, "ctrl3")


Cd8Krt8Mix_spotFeature <- function(Coordinate, Cd8cell, Krt8cell, AEcell, M){
	XLYL <- c(0,130)
	#g1_cell_names <- g1cell#unlist(lapply(g1cell, trimcellname))
	#g2_cell_names <- g2cell#unlist(lapply(g2cell, trimcellname))	
	#g3_cell_names <- g3cell#unlist(lapply(g3cell, trimcellname))	
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[Cd8cell] <- "red"
	color_feature[Krt8cell] <- "blue"
	color_feature[intersect(Cd8cell, Krt8cell)] <- "purple"
	color_feature[AEcell] <- "#AAAA00"
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}


pdf(file=paste0("spotFeature_CTscore_cd8v1.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_covid1, Krt_spots_use_covid1, Ae_spots_covid1, M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_covid2, Krt_spots_use_covid2, Ae_spots_covid2, M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_covid3, Krt_spots_use_covid3, Ae_spots_covid3, M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_ctrl1, Krt_spots_use_ctrl1, Ae_spots_ctrl1, M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_ctrl2, Krt_spots_use_ctrl2, Ae_spots_ctrl2, M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, Cd8_spots1_use_ctrl3, Krt_spots_use_ctrl3, Ae_spots_ctrl3, M="ctrl3")
dev.off()


pdf(file=paste0("spotFeature_CTscore_cd8v2.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots2_use_covid1, Krt_spots_use_covid1, Ae_spots_covid1, M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots2_use_covid2, Krt_spots_use_covid2, Ae_spots_covid2, M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots2_use_covid3, Krt_spots_use_covid3, Ae_spots_covid3, M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots2_use_ctrl1, Krt_spots_use_ctrl1, Ae_spots_ctrl1, M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots2_use_ctrl2, Krt_spots_use_ctrl2, Ae_spots_ctrl2, M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, Cd8_spots2_use_ctrl3, Krt_spots_use_ctrl3, Ae_spots_ctrl3, M="ctrl3")
dev.off()


pdf(file=paste0("spotFeature_CTscore_cd8v3.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots3_use_covid1, Krt_spots_use_covid1, Ae_spots_covid1, M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots3_use_covid2, Krt_spots_use_covid2, Ae_spots_covid2, M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots3_use_covid3, Krt_spots_use_covid3, Ae_spots_covid3, M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots3_use_ctrl1, Krt_spots_use_ctrl1, Ae_spots_ctrl1, M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots3_use_ctrl2, Krt_spots_use_ctrl2, Ae_spots_ctrl2, M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, Cd8_spots3_use_ctrl3, Krt_spots_use_ctrl3, Ae_spots_ctrl3, M="ctrl3")
dev.off()


pdf(file=paste0("spotFeature_CTscore_cd8v4.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots4_use_covid1, Krt_spots_use_covid1, Ae_spots_covid1, M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots4_use_covid2, Krt_spots_use_covid2, Ae_spots_covid2, M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots4_use_covid3, Krt_spots_use_covid3, Ae_spots_covid3, M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots4_use_ctrl1, Krt_spots_use_ctrl1, Ae_spots_ctrl1, M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots4_use_ctrl2, Krt_spots_use_ctrl2, Ae_spots_ctrl2, M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, Cd8_spots4_use_ctrl3, Krt_spots_use_ctrl3, Ae_spots_ctrl3, M="ctrl3")
dev.off()


pdf(file=paste0("spotFeature_CTscore_cd8v5.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots5_use_covid1, Krt_spots_use_covid1, Ae_spots_covid1, M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots5_use_covid2, Krt_spots_use_covid2, Ae_spots_covid2, M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots5_use_covid3, Krt_spots_use_covid3, Ae_spots_covid3, M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots5_use_ctrl1, Krt_spots_use_ctrl1, Ae_spots_ctrl1, M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots5_use_ctrl2, Krt_spots_use_ctrl2, Ae_spots_ctrl2, M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, Cd8_spots5_use_ctrl3, Krt_spots_use_ctrl3, Ae_spots_ctrl3, M="ctrl3")
dev.off()



pdf(file=paste0("spotFeature_CTscore_cd8v4_cd8only.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, Cd8_spots4_covid1, c(), c(), M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, Cd8_spots4_covid2, c(), c(), M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, Cd8_spots4_covid3, c(), c(), M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, Cd8_spots4_ctrl1, c(), c(), M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, Cd8_spots4_ctrl2, c(), c(), M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, Cd8_spots4_ctrl3, c(), c(), M="ctrl3")
dev.off()


pdf(file=paste0("spotFeature_CTscore_cd8v4_Krtonly.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_covid1, c(), M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_covid2, c(), M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_covid3, c(), M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_ctrl1, c(), M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_ctrl2, c(), M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_ctrl3, c(), M="ctrl3")
dev.off()


pdf(file=paste0("spotFeature_CTscore_cd8v4_Krtonly.pdf"),height=9,width=14)
par(bg="grey",mfrow=c(2,3),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(covid1_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_covid1, c(), M="covid1")
Cd8Krt8Mix_spotFeature(covid2_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_covid2, c(), M="covid2")
Cd8Krt8Mix_spotFeature(covid3_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_covid3, c(), M="covid3")
Cd8Krt8Mix_spotFeature(ctrl1_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_ctrl1, c(), M="ctrl1")
Cd8Krt8Mix_spotFeature(ctrl2_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_ctrl2, c(), M="ctrl2")
Cd8Krt8Mix_spotFeature(ctrl3_highQ_SCT@images$slice1@coordinates, c(), Krt_spots_ctrl3, c(), M="ctrl3")
dev.off()







pdf(file=paste0("spotFeature_ImmcolorOnly_exGlist.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(S1A_highQ_SCT@images$slice1@coordinates, Cd8on_S1A, c(), S1A_AE_cells_new, M="S1A")
Cd8Krt8Mix_spotFeature(S2A_highQ_SCT@images$slice1@coordinates, Cd8on_S2A, c(), S2A_AE_cells_new, M="S2A")
Cd8Krt8Mix_spotFeature(S1D_highQ_SCT@images$slice1@coordinates, Cd8on_S1D, c(), S1D_AE_cells_new, M="S1D")
Cd8Krt8Mix_spotFeature(S2D_highQ_SCT@images$slice1@coordinates, Cd8on_S2D, c(), S2D_AE_cells_new, M="S2D")
dev.off()


pdf(file=paste0("spotFeature_KrtcolorOnly_exGlist.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(S1A_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S1A,S1A_AE_cells_new, M="S1A")
Cd8Krt8Mix_spotFeature(S2A_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S2A,S2A_AE_cells_new, M="S2A")
Cd8Krt8Mix_spotFeature(S1D_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S1D,S1D_AE_cells_new, M="S1D")
Cd8Krt8Mix_spotFeature(S2D_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S2D,S2D_AE_cells_new, M="S2D")
dev.off()


# long list
#MacrophageAM_glist <- c("Rnase1","Vcan","Fcn1","S100a8","Coro1a","Tmem176b","Ms4a6a","Timp1","C15orf48","Zfp36l1","Ifitm3","Lgmn","Cd14","Tmem176a","Fgl2","Tymp","Fcgr2b","Emp3","Fpr3","Ccl18","Fyb1","Mafb","Ap1s2","Cebpd","Emp1","Trappc5","Ifitm2","Neat1","S100a9","Sod2","Ncf1","Apoe","Sgk1","Lilrb4","Ctsb","Csf1r","Tnfrsf1b","Cd99","Lmna","Plekho1","Klf6","Ahnak","Anpep","Blvrb","Hm13","Tgfbi","Lgals1","Pmp22","Lilrb2","Itgam","Rnase6","Vsir","Gapdh","Hif1a","Calhm6","Dusp6","Eno1","Csf3r","Lsp1","Cd36","Vamp5","Fos","S100a10","Rhoc","Clec5a","Tnfsf10","Cpm","Lims1","Cd84","Junb","Lhfpl2","Npc2","Ctsz","Hmga1","Rgs10","Il3ra","Ccr1","Mtrnr2l12","Tsc22d3","Pea15","Kctd12","Fxyd5","Plp2","Tpm4","Tnfsf13b","Il7r","Klf2","Plec","Zeb2","Eif4a1","Ptafr","Naip","Ptma","Zfp36l2","Pak1","Cybb","Sorl1","Mpeg1","Cd44","Ezr","Actb","Nfkbia","Gabarap","Zfp36","Sh3kbp1","Asgr1","Ppib","Ninj1","Trem2","App","Fgr","Ptpre","Capn2","Mcl1","Mgat1","Rap2b","Itgb2","Lilrb3","Ubc","Itga4","Lap3","Rbpj","Adgre5","Rala","Alcam","Aqp9","Rpl36a","Rps17","Slamf8","Slc43a3","Ier2","Samhd1","Ap2s1","Il10ra","Ahr","H3f3b","Il17ra","Vmp1","Fn1","Mthfd2","Prcp","C3ar1","Mettl9","Clec7a","Myadm","Fus","Sh3bgrl3","Dazap2","Cd300a","Tcirg1","Zyx","Myo1f","Adap2","Eef1g","Fpr1","Slc8a1","Icam1","Klf10","Trmt1","Jak1","Ifi16","Tpm3","Bcat1","Hmox1","Rab32","Irf1","Ifi30","Arhgap4","Dnph1","Psmb9","Krtcap2","Eif5a","Idh2","Actg1","Stat1","Rnf130","Glipr1","Smap2","Rpl22l1","Gnb4","Ifngr2","Cotl1","Rps7","Dnmt1","Aif1","Epsti1","Mpp1","Crip1","Ctsa","Hla-e","Dse","Tap1","Atf5","Apbb1ip","Ciita","Gbp1","Ywhah","Sec61g","Hnrnpu","Rb1","Gadd45b","Rnf213","Id3","Tpi1","Abi3","Metrnl","Cdkn1a","Lyn","Tln1","Tubb","Rps24","Cmtm3","Adam15","Eef1d","Eif4ebp1","Rplp1","Kiaa0930","H1fx","Agtrap","Rplp0","Pcbp1","Rpl39","Arf6","Rpl18","Capzb","Ralgds","Pdxk","Btg1","Srsf3","Rsrp1","Me2","Igflr1","Odf3b","Gnai2","Cnpy3","Prex1","Hnrnpa1","Iqgap1")
#AM_glist <- c("Spocd1","Dtx4","Sec11a","Cd46","Gga2","Acvrl1","Cldn7","Smim14","Coro2a","Galnt12","Myl12a","Ubash3b","Ac020656.1","Mgll","Cdc42ep3","Fabp5","Eef1a1","Lyz","Septin11","Itm2b","Tmsb4x","Tuba1b","Arhgef28","Msmo1","Ccnd3","Edem2","Lamtor1","Usp30-as1","Olr1","Acp5","Amigo2","Litaf","Lta4h","Akr1c3","Cd59","Lrpap1","Tbc1d4","Ethe1","Mpc1","Cd37","Cd101","Trhde","Mtpn","Tgoln2","Aqp3","Axl","Aoc3","Pnpla6","Stac","Fam3b","Ly86","Cyb5a","Folr3","Tpt1","Camp","Sparc","Slc47a1","Ms4a7","Acot4","Qsox1","Macc1","Tyrobp","Mt-cyb","11-sep","Lsamp","Acot2","Sh3bgrl","Apip","Rac2","Paral1","Kptn","Spn","Plbd1","Cfd","Ffar4","Smim25","Scp2","Scpep1","Glrx","Card16","Ptpmt1","Hcar2","Sod1","Myo6","Hddc2","Rasal2","Csf1","Cat","Folr1","Adtrp","Selenom","Ctsc","Decr1","Thbs1","Ptms","Alox5","Ca2","Apoc1","Dst","Irs2","Abcg1","Tcf7l2","Slc19a3","Ctsw","Tnni2","Fth1","Ubb","Stmn1","Lpl","Lima1","Bhlhe41","H2afz","Rhbdd2","Jpt1","Oasl","Cd81","Prss21","Ftl","C2","Pla2g16","Acot7","Abhd5","Mgst3","Plin2","Lgals3","C8b","Hbegf","Msr1","Stxbp2","Cyp27a1","Nmb","Tfrc","Pebp1","Ly6e","Itih5","Pparg","Fhl1","Rps20","Cxcl16","Fdft1","Evl","C9orf16","Cpe","Bsg","Apoc2","Fam89a","Ac026369.3","Hpgd","Svil","Cd9","Akr1b1","Aldh2","Plaat3","Vmo1","Alas1","Serpina1","Mlph","Aco1","Phlda3","Fbp1","Tspan3","Pdlim1","Atp1b1","C1qa","Trem1","Ppic","Gypc","C1qb","Alox5ap","Mcemp1","Gldn","S100a13","Scd","Linc02154","Tgm2","Lgals3bp","Ccl24","Mme","Nupr1","Ldhb","Hp","Igfbp2","Cd52","Ces1","Pcolce2","Gpd1","Inhba","Rbp4","Serping1","Fabp4","Ifi27")
#ProfibroticAM_glist <- c("A2m","Abca1","Abhd2","Ac145110.1","Acp2","Actb","Actg1","Adgre5","Alk","Anxa2","Ap000695.4","Apoc1","Apoe","Arhgdib","Arid5b","Arpc1b","Atox1","Atp13a3","Atp6v0d2","Atp6v1b2","Basp1","Bcat1","C15orf48","C6orf62","Calm3","Calr","Cap1","Ccl18","Ccl2","Ccl4","Cd14","Cd36","Cd4","Cd44","Cd48","Cd84","Cd9","Cdc42","Cfl1","Chi3l1","Chit1","Cltc","Cmklr1","Cmtm3","Cmtm6","Col4a2?as2","Cotl1","Cpm","Crip1","Crybb1","Csf1r","Ctsb","Ctsk","Ctsl","Ctsz","Cybb","Cyfip1","Dab2","Dfna5","Dopey2","Dpysl2","Emp1","Evl","Fam20c","Fcgr2a","Fcgr2b","Fcgr3a","Fcho2","Fgl2","Fkbp1a","Flna","Fmn1","Fn1","Fnip2","Folr2","Fpr3","Fuca1","Fxyd5","Fyb","Gbp1","Glipr1","Gm2a","Gpnmb","Gsdma","Gsn","Hck","Hcst","Hif1a","Hla-a","Hla-dqa2","Hmgn1","Hs3st2","Htra4","Igf2r","Ighg1","Ighg3","Ighg4","Igkc","Iglc2","Igsf6","Il1rn","Il7r","Itgam","Itgax","Itgb1","Itgb2","Kcnj5","Kiaa0930","L3mbtl4?as1","Lap3","Lasp1","Lcp1","Ldha","Lep","Lgmn","Lhfpl2","Lilrb4","Lipa","Litaf","Lsp1","Lst1","Mafb","Marcks","Mertk","Mmp7","Mmp9","Mnda","Mreg","Ms4a4a","Ms4a6a","Ms4a6e","Mt-atp6","Mt-atp8","Mt-co3","Mt-nd3","Mt-nd6","Nhsl2","Npc2","Nr1h3","Nrp2","Olr1","P4hb","Pdia4","Pea15","Pgk1","Pla2g7","Plac8","Plek","Ppt1","Psap","Qsox1","Rab31","Rala","Rgs1","Rnase1","Rp11-1143g9.4","Rp11-184m15.1","Rp11-20g13.3","Rp4?644l1.2","S100a8","S100a9","Samhd1","Scarb1","Scgb3a1","Sdc2","Sdc3","Sds","Sepp1","Sgk1","Sh3bgrl3","Siglec15","Sla","Slc16a10","Slc1a3","Slc26a11","Slc29a1","Slc29a3","Slc36a1","Smim3","Smpdl3a","Sod2","Sparc","Spp1","Srgn","Srsf2","Stat1","Tgfbi","Tgm2","Timp1","Tln1","Tm4sf19","Tm4sf19?tctex1","Tmem176b","Tmem37","Tmsb10","Tmsb4x","Tnfsf13b","Tpm3","Tpm4","Tpp1","Trem2","Ttyh3","Tuba1b","Tymp","Ucp2","Vamp5","Vil1","Vmp1","Wars","Wasf2","Zfp36l1","Zmiz1-as1","Znf385a")

#### curated_list 
#convertToMouseGeneSymbol <- function(gene_symbol) {
#  return(paste0(toupper(substr(gene_symbol, 1, 1)), 
#                tolower(substr(gene_symbol, 2, nchar(gene_symbol)))))
#}
#convertToHumanGeneSymbol <- function(gene_symbol) {
#  return(toupper(gene_symbol))
#}
#
# Example usage:
#mouse_gene = "Actb"
#human_gene = convertToHumanGeneSymbol(mouse_gene)
#print(human_gene)  # Expected output: "ACTB"
#
## Example
#gene <- "VCAN"
#convertToMouseGeneSymbol(gene) # Outputs: Vcan
#




######### combine samples
#cbdata <- merge(covid1_highQ, y=c(covid2_highQ,covid3_highQ,ctrl1_highQ,ctrl2_highQ,ctrl3_highQ),add.cell.ids=c("covid1","covid2","covid3","ctrl1","ctrl2","ctrl3"),project="Spatial")

#trimsamplename <- function(inname){
#  return(strsplit(inname, "_")[[1]][1])
#}
#usesample <- unlist(lapply(names(cbdata$orig.ident), trimsamplename))
#
#cbdata$sample <- usesample
#cbdata <- RunPCA(cbdata, assay = "Spatial", verbose = FALSE)
#cbdata <- FindNeighbors(cbdata, reduction = "pca", dims = 1:30)
#cbdata <- FindClusters(cbdata, verbose = FALSE)
#cbdata <- RunUMAP(cbdata, reduction = "pca", dims = 1:30)
#
#
#pdf(file="cbdata_UMAP_features.pdf",width=8,height=4)
#p1 <- DimPlot(cbdata, label = TRUE,group.by="sample") + ggtitle("sample")
#p2 <- DimPlot(cbdata, reduction = "umap", group.by = "seurat_clusters")
##p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
#p1+p2
#dev.off()




#### harmony
#cbdata <- merge(S1A_highQ, y=c(S1D_highQ,S2A_highQ,S2D_highQ),add.cell.ids=c("S1A","S1D","S2A","S2D"),project="spatial")
#
#trimsamplename <- function(inname){
#  return(strsplit(inname, "_")[[1]][1])
#}
#usesample <- unlist(lapply(names(cbdata$orig.ident), trimsamplename))
#cbdata$sample <- usesample
#cbdata$replicate <- rep("S1",length(cbdata$sample))
#cbdata$replicate[which(cbdata$sample %in% c("S2A","S2D"))] <- "S2"
#
##brain.combined <- NormalizeData(brain.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
##brain.combined <- RunUMAP(brain.combined, dims = 1:30)
#
#cbdata.combined <- SCTransform(cbdata, assay = "Spatial", verbose = FALSE)
#cbdata.combined <- RunPCA(cbdata.combined, npcs = 30, verbose = FALSE, assay="SCT")
#cbdataHAR <- RunHarmony(cbdata.combined, group.by.vars = "replicate")
#cbdataHAR <- RunUMAP(cbdataHAR, reduction = "harmony", dims = 1:30)
#cbdataHAR <- FindNeighbors(cbdataHAR, reduction = "harmony", dims = 1:30)
##cbdataHAR <- FindClusters(cbdataHAR, resolution = 0.5)
#
#
#pdf(file="cbdata_SCT_HarmonybatchRemove_UMAP_features.pdf",width=8,height=4)
#p1 <- DimPlot(cbdataHAR, label = TRUE,group.by="sample") + ggtitle("sample")
#p2 <- DimPlot(cbdataHAR, reduction = "umap", group.by = "seurat_clusters")
##p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
#p1+p2
#dev.off()



#cbdataHAR2 <- RunHarmony(cbdata, group.by.vars = "sample")
#cbdataHAR2 <- RunUMAP(cbdataHAR2, reduction = "harmony", dims = 1:30)
#cbdataHAR2 <- FindClusters(cbdataHAR2, verbose = FALSE)
#
##cbdataHAR2dis <- RunHarmony(cbdata, group.by.vars = "disease")
##cbdataHAR2dis <- RunUMAP(cbdataHAR2dis, reduction = "harmony", dims = 1:30)
#
#
#pdf(file="cbdata_SCTharmony_UMAP_features.pdf",width=12,height=4)
#p1 <- DimPlot(cbdataHAR2, label = TRUE,group.by="sample" ,reduction = "umap") + ggtitle("Harmony sample")
#p2 <- DimPlot(cbdataHAR2, label = TRUE,group.by="disease" ,reduction = "umap") + ggtitle("Harmony disease")
#p3 <- DimPlot(cbdataHAR2, label = TRUE,group.by="seurat_clusters" ,reduction = "umap") + ggtitle("Harmony cluster")
##p3 <- DimPlot(cbdataHAR2dis, label = TRUE,group.by="sample" ,reduction = "umap") + ggtitle("HarmonyDisease sample")
##p4 <- DimPlot(cbdataHAR2dis, label = TRUE,group.by="disease" ,reduction = "umap") + ggtitle("HarmonyDisease disease")
##p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
#p1+p2+p3
#dev.off()


#cbdataHAR2_clusterRes0.5 <- FindClusters(cbdataHAR2, resolution = 0.5)
#cbdataHAR2_clusterRes0.6 <- FindClusters(cbdataHAR2, resolution = 0.6)
#cbdataHAR2_clusterRes0.7 <- FindClusters(cbdataHAR2, resolution = 0.7)
#cbdataHAR2_clusterRes0.8 <- FindClusters(cbdataHAR2, resolution = 0.8)
#cbdataHAR2_clusterRes0.9 <- FindClusters(cbdataHAR2, resolution = 0.9)
#cbdataHAR2_clusterRes1.0 <- FindClusters(cbdataHAR2, resolution = 1.0)
#cbdataHAR2_clusterRes1.1 <- FindClusters(cbdataHAR2, resolution = 1.1)
#cbdataHAR2_clusterRes1.2 <- FindClusters(cbdataHAR2, resolution = 1.2)
#cbdataHAR2_clusterRes1.3 <- FindClusters(cbdataHAR2, resolution = 1.3)
#cbdataHAR2_clusterRes1.4 <- FindClusters(cbdataHAR2, resolution = 1.4)
#cbdataHAR2_clusterRes1.5 <- FindClusters(cbdataHAR2, resolution = 1.5)
#cbdataHAR2_clusterRes1.6 <- FindClusters(cbdataHAR2, resolution = 1.6)
#cbdataHAR2_clusterRes1.7 <- FindClusters(cbdataHAR2, resolution = 1.7)
#cbdataHAR2_clusterRes1.8 <- FindClusters(cbdataHAR2, resolution = 1.8)
#cbdataHAR2_clusterRes1.9 <- FindClusters(cbdataHAR2, resolution = 1.9)
#cbdataHAR2_clusterRes2.0 <- FindClusters(cbdataHAR2, resolution = 2.0)
#
##pdf(file="cbdata_SCT_HarmonySamplebatchRemove_UMAP_features.pdf",width=4,height=4)
#p1 <- DimPlot(cbdataHAR2, label = TRUE,group.by="sample") + ggtitle("sample")
##p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
#p1
#dev.off()



#p0.5 <- DimPlot(cbdataHAR2_clusterRes0.5, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.5")
#p0.6 <- DimPlot(cbdataHAR2_clusterRes0.6, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.6")
#p0.7 <- DimPlot(cbdataHAR2_clusterRes0.7, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.7")
#p0.8 <- DimPlot(cbdataHAR2_clusterRes0.8, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.8")
#p0.9 <- DimPlot(cbdataHAR2_clusterRes0.9, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.9")
#p1.0 <- DimPlot(cbdataHAR2_clusterRes1.0, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.0")
#p1.1 <- DimPlot(cbdataHAR2_clusterRes1.1, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.1")
#p1.2 <- DimPlot(cbdataHAR2_clusterRes1.2, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.2")
#p1.3 <- DimPlot(cbdataHAR2_clusterRes1.3, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.3")
#p1.4 <- DimPlot(cbdataHAR2_clusterRes1.4, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.4")
#p1.5 <- DimPlot(cbdataHAR2_clusterRes1.5, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.5")
#p1.6 <- DimPlot(cbdataHAR2_clusterRes1.6, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.6")
#p1.7 <- DimPlot(cbdataHAR2_clusterRes1.7, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.7")
#p1.8 <- DimPlot(cbdataHAR2_clusterRes1.8, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.8")
#p1.9 <- DimPlot(cbdataHAR2_clusterRes1.9, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.9")
#p2.0 <- DimPlot(cbdataHAR2_clusterRes2.0, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res2.0")
#
#pdf(file="cbdataSCTharmony_clusterRes_UMAP_features.pdf",width=20,height=20)
#plot_grid(p0.5,
#p0.6,
#p0.7,
#p0.8,
#p0.9,
#p1.0,
#p1.1,
#p1.2,
#p1.3,
#p1.4,
#p1.5,
#p1.6,
#p1.7,
#p1.8,
#p1.9,
#p2.0,nrow=4,ncol=4)
#dev.off()
#

#pdf(file="cbdataSCTharmony_clusterRes2.0_UMAP_features.pdf",width=5,height=4)
#p2.0
#dev.off()
#
#cbdataHAR2$cluster <- cbdataHAR2_clusterRes2.0$seurat_clusters



### cluster v.s. sample
#confmat <- function(inmat){
#  feature1 <- unique(inmat[,1])
#  feature2 <- unique(inmat[,2])
#  outmat <- matrix(rep(0, length(feature1)*length(feature2)),nrow=length(feature1),ncol=length(feature2))
#  rownames(outmat) <- feature1
#  colnames(outmat) <- feature2
#  for(i in 1:nrow(inmat)){
#    outmat[inmat[i,1], inmat[i,2]] <- outmat[inmat[i,1], inmat[i,2]] + 1
#  }
#  return(outmat)
#}
#
#confmat_sample_cluster <- confmat(cbind(as.character(cbdataHAR2$sample),cbdataHAR2$seurat_clusters))#[,as.character(1:13)]
#colnames(confmat_sample_cluster) <- paste0("C", as.numeric(colnames(confmat_sample_cluster))-1)
#confmat_sample_cluster_percent <- confmat_sample_cluster / apply(confmat_sample_cluster, 1, sum)
#
#usecolor <- c("white","red")
#ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
#pdf(file="confmat_sampleVScluster_heatmap.pdf")
#heatmap.2(confmat_sample_cluster_percent, trace="none", col=ColorRamp,main="%cell in sample")
#dev.off()
#
#write.table(confmat_sample_cluster, file="confmat_sample_cluster.txt",row.names=T,col.names=T,sep="\t",quote=F)
#epi_cells: C4,C7
#Imm_cells: C2,  C12, C13


#### curated_list 
#convertToMouseGeneSymbol <- function(gene_symbol) {
#  return(paste0(toupper(substr(gene_symbol, 1, 1)), 
#                tolower(substr(gene_symbol, 2, nchar(gene_symbol)))))
#}

#cbdataHAR2_scale_allgenes <-  SCTransform(cbdataHAR2, assay = "Spatial", verbose = FALSE, ,variable.features.n = 36156)
#
#
##saveRDS(cbdataHAR2_scale_allgenes, file="cbdataHAR2_scale_allgenes.rds")
##cbdataHAR2 <- readRDS(file="cbdataHAR2.rds")
#cbdataHAR2_scale_allgenes <- readRDS(file="cbdataHAR2_scale_allgenes.rds")
#
#dim(cbdataHAR2@assays$SCT@data)
#dim(cbdataHAR2_scale_allgenes@assays$SCT@scale.data)
#

#expmat_scale <- cbdataHAR2_scale_allgenes@assays$SCT@scale.data
#keygeneALL_scale <- intersect( toupper(unique(c(
#CD8_glist ,
#Krt_glist ,
#AT1_glist   ,
#AT2_glist ,
#Mono_glist,
#MDM_glist  ,
#IM_glist  ,
#AM_glist ))),rownames(expmat_scale))
#
#keygeneALL_scale_major <- intersect( toupper(unique(c(
#CD8_glist ,
#Krt_glist ,
#AT1_glist   ,
#AT2_glist 
#))),rownames(expmat_scale))
#
#
#aveExpMat_scale <-  matrix(rep(0, length(keygeneALL_scale)*length(unique(cbdataHAR2$seurat_clusters)) ),nrow=length(keygeneALL_scale))
#rownames(aveExpMat_scale) <- keygeneALL_scale
#colnames(aveExpMat_scale) <- paste0("C",0:(length(unique(cbdataHAR2$seurat_clusters)) -1))
#medianExpMat_scale <- aveExpMat_scale
#for(i in 0:(length(unique(cbdataHAR2$seurat_clusters)) -1)){
#	usecell <- names(cbdataHAR2$seurat_clusters)[which(cbdataHAR2$seurat_clusters==i)]
#	tempMat <- expmat_scale[keygeneALL_scale,usecell]
#	aveExpMat_scale[, paste0("C",i)] <- apply(tempMat, 1,mean)
#	medianExpMat_scale[, paste0("C",i)] <- apply(tempMat, 1,median)
#}
#
#

#usecolor <- c("blue","white","red")
#ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
#pdf(file="keygeneExp_cluster_heatmap.pdf",height=12)
#heatmap.2(aveExpMat_scale, trace="none",Rowv=F, col=ColorRamp,main="aveExp scale",cexRow=0.6)
#dev.off()
#
#pdf(file="keygeneExpMajor_cluster_heatmap.pdf")
#heatmap.2(aveExpMat_scale[keygeneALL_scale_major,], trace="none",Rowv=F, col=ColorRamp,main="aveExp scale",cexRow=1)
#dev.off()
#



#CD8MKG<-read.table("CD8Tmarker_FCge2.txt",row.names=1,header=T)
#CD8MKG_top10 <- rownames(CD8MKG)[1:10]
#CD8MKG_top20 <- rownames(CD8MKG)[1:20]
#CD8MKG_FCge2 <- rownames(CD8MKG)[which(CD8MKG[,"avg_log2FC"]>=1)]#[1:10]
#