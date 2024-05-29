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
S1A <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/visium/mouseLungCD8/10x_analysis_6433-HN/Sample_6433-HN-S1-A-GEX_TCAGAACT-GTGGACTG")
S1D <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/visium/mouseLungCD8/10x_analysis_6433-HN/Sample_6433-HN-S1-D-GEX_CGGCCATA-GAACGAAG")
S2A <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/visium/mouseLungCD8/10x_analysis_6433-HN/Sample_6433-HN-S2-A-GEX_CGCGCCAA-CCTTGGAT")
S2D <- Load10X_Spatial("/project/zanglab_project/sh8tv/SunLab/Data/visium/mouseLungCD8/10x_analysis_6433-HN/Sample_6433-HN-S2-D-GEX_ACTAAATG-ACAGATTT")


sampleName <- c("S1A","S1D","S2A","S2D")
sampleCol <- c("red","blue","darkred","darkblue")

pdf(file="sample_legend.pdf")
plot(1,1,type="n",axes=F,xlab="",ylab="")
legend("center",legend=sampleName,col=sampleCol,bty="n",cex=2,lwd=5)
dev.off()


### QC

pdf(file="sampleQC_outlier.pdf",width=12,height=4)
par(mfrow=c(1,3),mar=c(4,4,2,2))
barplot(c(length(S1A$nCount_Spatial), 
	    length(S1D$nCount_Spatial),
	    length(S2A$nCount_Spatial),
	    length(S2D$nCount_Spatial)),
	    names=sampleName,ylab="#spot",main="#spot per spot",col=sampleCol)
boxplot(log10(as.numeric(S1A$nCount_Spatial)), 
	    log10(as.numeric(S1D$nCount_Spatial)),
	    log10(as.numeric(S2A$nCount_Spatial)),
	    log10(as.numeric(S2D$nCount_Spatial)),
	    names=sampleName,ylab="log10 reads count",main="#reads per spot",col=sampleCol, outline=T)

boxplot(log10(as.numeric(S1A$nFeature_Spatial)), 
	    log10(as.numeric(S1D$nFeature_Spatial)),
	    log10(as.numeric(S2A$nFeature_Spatial)),
	    log10(as.numeric(S2D$nFeature_Spatial)),
	    names=sampleName,ylab="log10 covered gene",main="#genes per spot",col=sampleCol, outline=T)
abline(h=log10(2000),lwd=2)
dev.off()


pdf(file="sampleQC.pdf",width=12,height=4)
par(mfrow=c(1,3),mar=c(4,4,2,2))
barplot(c(length(S1A$nCount_Spatial), 
	    length(S1D$nCount_Spatial),
	    length(S2A$nCount_Spatial),
	    length(S2D$nCount_Spatial)),
	    names=sampleName,ylab="#spot",main="#spot per spot",col=sampleCol)
boxplot(log10(as.numeric(S1A$nCount_Spatial)), 
	    log10(as.numeric(S1D$nCount_Spatial)),
	    log10(as.numeric(S2A$nCount_Spatial)),
	    log10(as.numeric(S2D$nCount_Spatial)),
	    names=sampleName,ylab="log10 reads count",main="#reads per spot",col=sampleCol, outline=F)

boxplot(log10(as.numeric(S1A$nFeature_Spatial)), 
	    log10(as.numeric(S1D$nFeature_Spatial)),
	    log10(as.numeric(S2A$nFeature_Spatial)),
	    log10(as.numeric(S2D$nFeature_Spatial)),
	    names=sampleName,ylab="log10 covered gene",main="#genes per spot",col=sampleCol, outline=F)
dev.off()

### QC filtering
# nFeature >= 2000
S1A_highQ <- S1A[,which(S1A$nFeature_Spatial>=2000)]
S1D_highQ <- S1D[,which(S1D$nFeature_Spatial>=2000)]
S2A_highQ <- S2A[,which(S2A$nFeature_Spatial>=2000)]
S2D_highQ <- S2D[,which(S2D$nFeature_Spatial>=2000)]

#saveRDS(S1A_highQ, file="S1A_highQ.rds")
#saveRDS(S1D_highQ, file="S1D_highQ.rds")
#saveRDS(S2A_highQ, file="S2A_highQ.rds")
#saveRDS(S2D_highQ, file="S2D_highQ.rds")
#
#S1A_highQ <- readRDS(file="S1A_highQ.rds")
#S1D_highQ <- readRDS(file="S1D_highQ.rds")
#S2A_highQ <- readRDS(file="S2A_highQ.rds")
#S2D_highQ <- readRDS(file="S2D_highQ.rds")

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

pdf(file="spotFeature_nCount.pdf",height=9,width=9)
par(bg="grey",mfrow=c(2,2),mar=c(4,4,2,2))
spotFeature_plot(S1A_highQ@images$slice1@coordinates, S1A_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S1A log10 nCount")
spotFeature_plot(S1D_highQ@images$slice1@coordinates, S1D_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S1D log10 nCount")
spotFeature_plot(S2A_highQ@images$slice1@coordinates, S2A_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S2A log10 nCount")
spotFeature_plot(S2D_highQ@images$slice1@coordinates, S2D_highQ$nCount_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S2D log10 nCount")
dev.off()


pdf(file="spotFeature_nGene.pdf",height=9,width=9)
par(bg="grey",mfrow=c(2,2),mar=c(4,4,2,2))
spotFeature_plot(S1A_highQ@images$slice1@coordinates, S1A_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S1A log10 nGene")
spotFeature_plot(S1D_highQ@images$slice1@coordinates, S1D_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S1D log10 nGene")
spotFeature_plot(S2A_highQ@images$slice1@coordinates, S2A_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S2A log10 nGene")
spotFeature_plot(S2D_highQ@images$slice1@coordinates, S2D_highQ$nFeature_Spatial, usecolor = c("blue","white","red"), log_scale = 1, M="S2D log10 nGene")
dev.off()


######### combine samples
cbdata <- merge(S1A_highQ, y=c(S1D_highQ,S2A_highQ,S2D_highQ),add.cell.ids=c("S1A","S1D","S2A","S2D"),project="spatial")

trimsamplename <- function(inname){
  return(strsplit(inname, "_")[[1]][1])
}
usesample <- unlist(lapply(names(cbdata$orig.ident), trimsamplename))

cbdata$sample <- usesample

condition <- rep("Cd8a", length(cbdata$sample))
condition[which(cbdata$sample %in% c("S1A","S2A"))] <- "IgG"
cbdata$condition <- condition


cbdata <- RunPCA(cbdata, assay = "SCT", verbose = FALSE)
cbdata <- FindNeighbors(cbdata, reduction = "pca", dims = 1:30)
cbdata <- FindClusters(cbdata, verbose = FALSE)
cbdata <- RunUMAP(cbdata, reduction = "pca", dims = 1:30)


pdf(file="mouse_cbdata_SCTraw_UMAP_features.pdf",width=12,height=4)
p1 <- DimPlot(cbdata, label = TRUE,group.by="sample") + ggtitle("sample")
p2 <- DimPlot(cbdata, label = TRUE,group.by="condition", reduction = "umap") + ggtitle("condition")
p3 <- DimPlot(cbdata, reduction = "umap", group.by = "seurat_clusters")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1+p2+p3
dev.off()

#pdf(file="cbdata_UMAP_features.pdf",width=8,height=4)
#p1 <- DimPlot(cbdata, label = TRUE,group.by="sample") + ggtitle("sample")
#p2 <- DimPlot(cbdata, reduction = "umap", group.by = "seurat_clusters")
##p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
#p1+p2
#dev.off()




#cbdataHAR2_scale_allgenes <- readRDS(file="cbdataHAR2_scale_allgenes.rds")
#DefaultAssay(cbdata) <- "SCT"
#condition <- rep("Cd8a", length(cbdataHAR2_scale_allgenes$sample))
#condition[which(cbdataHAR2_scale_allgenes$sample %in% c("S1A","S2A"))] <- "IgG"
#cbdataHAR2_scale_allgenes$condition <- condition
#
#
#pdf(file="mouse_cbdataHAR_SCTraw_UMAP_features.pdf",width=12,height=4)
#p1 <- DimPlot(cbdataHAR2_scale_allgenes, label = TRUE,group.by="sample") + ggtitle("sample")
#p2 <- DimPlot(cbdataHAR2_scale_allgenes, label = TRUE,group.by="condition", reduction = "umap") + ggtitle("condition")
#p3 <- DimPlot(cbdataHAR2_scale_allgenes, reduction = "umap", group.by = "cluster")
##p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
#p1+p2+p3
#dev.off()


###### condition wise DEG
DefaultAssay(cbdata) <- "SCT"
condition_DEG <- FindMarkers(cbdata,group.by="condition", ident.1 = "Cd8a", ident.2 = "IgG")
#condition_DEG_filter <- condition_DEG[which(abs(condition_DEG[,"avg_log2FC"]) >= 1 & condition_DEG[,"p_val_adj"]< 0.01),]
top10_condition_DEG <- condition_DEG[order(condition_DEG[,"avg_log2FC"]),][c(1:10, (nrow(condition_DEG)-9):nrow(condition_DEG)),]

pdf(file="mouse_UMAP_conditionDEGtop10.pdf",width=12,height=15)
FeaturePlot(cbdata, features=rownames(top10_condition_DEG), reduction="umap")
dev.off()

write.table(condition_DEG, file="mouse_condition_DEG.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(rownames(condition_DEG)[which(condition_DEG[,"avg_log2FC"]>0)], file="mouse_condition_Cd8aHigh_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(rownames(condition_DEG)[which(condition_DEG[,"avg_log2FC"]<0)], file="mouse_condition_IgGHigh_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)


###### mouse AE vs Krt DEG
AE_cells_new <- colnames(cbdataHAR2_scale_allgenes)[which(cbdataHAR2_scale_allgenes$cluster %in% c("0","5","14","18","3","10","19","22","23","1","8","20","26","17","6","9","24"))]

scaleMat <- cbdataHAR2_scale_allgenes@assays$SCT@scale.data
ImmExp <- apply(scaleMat[c("Cd8a","Cx3cr1","Mertk"),],2,mean)
KrtExp <- apply(scaleMat[c("Krt8","Krt5","Krt17"),],2,mean)
Cd8on_spots <- colnames(scaleMat)[which(ImmExp > -0.4)]
Krt8on_spots <- colnames(scaleMat)[which(KrtExp > -0.4)]

celltype3 <- rep("AE",ncol(cbdataHAR2_scale_allgenes))
names(celltype3) <- colnames(cbdataHAR2_scale_allgenes)
celltype3[Cd8on_spots] <- "Cd8"
celltype3[Krt8on_spots] <- "Krt"
celltype3[AE_cells_new] <- "AE"
cbdataHAR2_scale_allgenes$celltype3 <- celltype3

celltype3 <- rep("AE",ncol(cbdata))
names(celltype3) <- colnames(cbdata)
celltype3[Cd8on_spots] <- "Cd8"
celltype3[Krt8on_spots] <- "Krt"
celltype3[AE_cells_new] <- "AE"
cbdata$celltype3 <- celltype3


pdf(file="mouse_cbdata_SCTraw_UMAP_features_celltyp3.pdf",width=4,height=4)
p1 <- DimPlot(cbdata, label = TRUE,group.by="celltype3",cols=c("#AAAA00","red","blue")) + ggtitle("CT3")
#p2 <- DimPlot(cbdata, label = TRUE,group.by="condition", reduction = "umap") + ggtitle("condition")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "seurat_clusters")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1
dev.off()

pdf(file="mouse_cbdataHAR_SCTraw_UMAP_features_celltyp3.pdf",width=4,height=4)
p1 <- DimPlot(cbdataHAR2_scale_allgenes, label = TRUE,group.by="celltype3",cols=c("#AAAA00","red","blue")) + ggtitle("CT3")
#p2 <- DimPlot(cbdata, label = TRUE,group.by="condition", reduction = "umap") + ggtitle("condition")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "seurat_clusters")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1
dev.off()


AEKrt_DEG <- FindMarkers(cbdata,group.by="celltype3", ident.1 = "AE", ident.2 = "Krt")
#condition_DEG_filter <- condition_DEG[which(abs(condition_DEG[,"avg_log2FC"]) >= 1 & condition_DEG[,"p_val_adj"]< 0.01),]
top10_AEKrt_DEG <- AEKrt_DEG[order(AEKrt_DEG[,"avg_log2FC"]),][c(1:10, (nrow(AEKrt_DEG)-9):nrow(AEKrt_DEG)),]

pdf(file="mouse_UMAP_AEKrtDEGtop10.pdf",width=12,height=15)
FeaturePlot(cbdata, features=rownames(top10_AEKrt_DEG), reduction="umap")
dev.off()

write.table(AEKrt_DEG, file="mouse_AEKrt_DEG.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(rownames(AEKrt_DEG)[which(AEKrt_DEG[,"avg_log2FC"]>0)], file="mouse_celltype3_AEHigh_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(rownames(AEKrt_DEG)[which(AEKrt_DEG[,"avg_log2FC"]<0)], file="mouse_celltype3_KrtHigh_DEGname.txt",row.names=F,col.names=F,sep="\t",quote=F)

AEKrt_DEG <- read.table("mouse_AEKrt_DEG.txt",row.names=1,header=T)



### harmony
cbdata <- merge(S1A_highQ, y=c(S1D_highQ,S2A_highQ,S2D_highQ),add.cell.ids=c("S1A","S1D","S2A","S2D"),project="spatial")

trimsamplename <- function(inname){
  return(strsplit(inname, "_")[[1]][1])
}
usesample <- unlist(lapply(names(cbdata$orig.ident), trimsamplename))
cbdata$sample <- usesample
cbdata$replicate <- rep("S1",length(cbdata$sample))
cbdata$replicate[which(cbdata$sample %in% c("S2A","S2D"))] <- "S2"

#brain.combined <- NormalizeData(brain.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
#brain.combined <- RunUMAP(brain.combined, dims = 1:30)

cbdata.combined <- SCTransform(cbdata, assay = "Spatial", verbose = FALSE)
cbdata.combined <- RunPCA(cbdata.combined, npcs = 30, verbose = FALSE, assay="SCT")
cbdataHAR <- RunHarmony(cbdata.combined, group.by.vars = "replicate")
cbdataHAR <- RunUMAP(cbdataHAR, reduction = "harmony", dims = 1:30)
cbdataHAR <- FindNeighbors(cbdataHAR, reduction = "harmony", dims = 1:30)
#cbdataHAR <- FindClusters(cbdataHAR, resolution = 0.5)


pdf(file="cbdata_SCT_HarmonybatchRemove_UMAP_features.pdf",width=8,height=4)
p1 <- DimPlot(cbdataHAR, label = TRUE,group.by="sample") + ggtitle("sample")
p2 <- DimPlot(cbdataHAR, reduction = "umap", group.by = "seurat_clusters")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1+p2
dev.off()



cbdata.combined <- SCTransform(cbdata, assay = "Spatial", verbose = FALSE)
cbdata.combined <- RunPCA(cbdata.combined, npcs = 30, verbose = FALSE, assay="SCT")
cbdataHAR2 <- RunHarmony(cbdata.combined, group.by.vars = "sample")
cbdataHAR2 <- RunUMAP(cbdataHAR2, reduction = "harmony", dims = 1:30)
cbdataHAR2_clusterRes0.5 <- FindClusters(cbdataHAR2, resolution = 0.5)
cbdataHAR2_clusterRes0.6 <- FindClusters(cbdataHAR2, resolution = 0.6)
cbdataHAR2_clusterRes0.7 <- FindClusters(cbdataHAR2, resolution = 0.7)
cbdataHAR2_clusterRes0.8 <- FindClusters(cbdataHAR2, resolution = 0.8)
cbdataHAR2_clusterRes0.9 <- FindClusters(cbdataHAR2, resolution = 0.9)
cbdataHAR2_clusterRes1.0 <- FindClusters(cbdataHAR2, resolution = 1.0)
cbdataHAR2_clusterRes1.1 <- FindClusters(cbdataHAR2, resolution = 1.1)
cbdataHAR2_clusterRes1.2 <- FindClusters(cbdataHAR2, resolution = 1.2)
cbdataHAR2_clusterRes1.3 <- FindClusters(cbdataHAR2, resolution = 1.3)
cbdataHAR2_clusterRes1.4 <- FindClusters(cbdataHAR2, resolution = 1.4)
cbdataHAR2_clusterRes1.5 <- FindClusters(cbdataHAR2, resolution = 1.5)
cbdataHAR2_clusterRes1.6 <- FindClusters(cbdataHAR2, resolution = 1.6)
cbdataHAR2_clusterRes1.7 <- FindClusters(cbdataHAR2, resolution = 1.7)
cbdataHAR2_clusterRes1.8 <- FindClusters(cbdataHAR2, resolution = 1.8)
cbdataHAR2_clusterRes1.9 <- FindClusters(cbdataHAR2, resolution = 1.9)
cbdataHAR2_clusterRes2.0 <- FindClusters(cbdataHAR2, resolution = 2.0)

pdf(file="cbdata_SCT_HarmonySamplebatchRemove_UMAP_features.pdf",width=4,height=4)
p1 <- DimPlot(cbdataHAR2, label = TRUE,group.by="sample") + ggtitle("sample")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1
dev.off()



p0.5 <- DimPlot(cbdataHAR2_clusterRes0.5, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.5")
p0.6 <- DimPlot(cbdataHAR2_clusterRes0.6, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.6")
p0.7 <- DimPlot(cbdataHAR2_clusterRes0.7, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.7")
p0.8 <- DimPlot(cbdataHAR2_clusterRes0.8, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.8")
p0.9 <- DimPlot(cbdataHAR2_clusterRes0.9, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res0.9")
p1.0 <- DimPlot(cbdataHAR2_clusterRes1.0, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.0")
p1.1 <- DimPlot(cbdataHAR2_clusterRes1.1, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.1")
p1.2 <- DimPlot(cbdataHAR2_clusterRes1.2, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.2")
p1.3 <- DimPlot(cbdataHAR2_clusterRes1.3, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.3")
p1.4 <- DimPlot(cbdataHAR2_clusterRes1.4, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.4")
p1.5 <- DimPlot(cbdataHAR2_clusterRes1.5, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.5")
p1.6 <- DimPlot(cbdataHAR2_clusterRes1.6, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.6")
p1.7 <- DimPlot(cbdataHAR2_clusterRes1.7, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.7")
p1.8 <- DimPlot(cbdataHAR2_clusterRes1.8, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.8")
p1.9 <- DimPlot(cbdataHAR2_clusterRes1.9, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res1.9")
p2.0 <- DimPlot(cbdataHAR2_clusterRes2.0, reduction = "umap", group.by = "seurat_clusters") + ggtitle("res2.0")

pdf(file="cbdataSCTharmony_clusterRes_UMAP_features.pdf",width=20,height=20)
plot_grid(p0.5,
p0.6,
p0.7,
p0.8,
p0.9,
p1.0,
p1.1,
p1.2,
p1.3,
p1.4,
p1.5,
p1.6,
p1.7,
p1.8,
p1.9,
p2.0,nrow=4,ncol=4)
dev.off()


pdf(file="cbdataSCTharmony_clusterRes2.0_UMAP_features.pdf",width=5,height=4)
p2.0
dev.off()

cbdataHAR2$cluster <- cbdataHAR2_clusterRes2.0$seurat_clusters



### cluster v.s. sample
confmat <- function(inmat){
  feature1 <- unique(inmat[,1])
  feature2 <- unique(inmat[,2])
  outmat <- matrix(rep(0, length(feature1)*length(feature2)),nrow=length(feature1),ncol=length(feature2))
  rownames(outmat) <- feature1
  colnames(outmat) <- feature2
  for(i in 1:nrow(inmat)){
    outmat[inmat[i,1], inmat[i,2]] <- outmat[inmat[i,1], inmat[i,2]] + 1
  }
  return(outmat)
}

confmat_sample_cluster <- confmat(cbind(as.character(cbdataHAR2$sample),cbdataHAR2$cluster))#[,as.character(1:13)]
colnames(confmat_sample_cluster) <- paste0("C", as.numeric(colnames(confmat_sample_cluster))-1)
confmat_sample_cluster_percent <- confmat_sample_cluster / apply(confmat_sample_cluster, 1, sum)

usecolor <- c("white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
pdf(file="confmat_sampleVScluster_heatmap.pdf")
heatmap.2(confmat_sample_cluster_percent, trace="none", col=ColorRamp,main="%cell in sample")
dev.off()

write.table(confmat_sample_cluster, file="confmat_sample_cluster.txt",row.names=T,col.names=T,sep="\t",quote=F)
#epi_cells: C4,C7
#Imm_cells: C2,  C12, C13


### RNA projection

anchors <- FindTransferAnchors(reference = scRNA_process, query = cbdataHAR2, normalization.method = "SCT")

predictions_VagueCluster <- TransferData(anchorset = anchors, refdata = scRNA_process$Vague_ident, prediction.assay = TRUE,
    weight.reduction = cbdataHAR2[["pca"]], dims = 1:30)
predictions_Cluster <- TransferData(anchorset = anchors, refdata = scRNA_process$ident, prediction.assay = TRUE,
    weight.reduction = cbdataHAR2[["pca"]], dims = 1:30)

score_matVague <- predictions_VagueCluster@data
score_mat <- predictions_Cluster@data

pdf(file="scRNA_max_projectionScore_hist.pdf",width=8,height=4)
par(mar=c(4,4,2,2),mfrow=c(1,2))
hist(as.numeric(score_matVague["max",]),n=200, xlab="max projection score per spot", main="scRNA vague cluster")
abline(v=0.75,lwd=2,col="red")
hist(as.numeric(score_mat["max",]),n=200, xlab="max projection score per spot", main="scRNA cluster")
abline(v=0.75,lwd=2,col="red")
dev.off()

projection_cutoff = 0.4
assign_projection <- function(SC_mat, pcutoff){
	projection_CT <- rep("NA",ncol(SC_mat))
	names(projection_CT) <- colnames(SC_mat)
	for(i in colnames(SC_mat)){
		tmpdata <- SC_mat[, i]
		if(tmpdata["max"] >= pcutoff){
			maxCTname <- names(which(tmpdata[-length(tmpdata)] == tmpdata[length(tmpdata)]))[1]
			projection_CT[i] <- maxCTname
		}
	}
	return(projection_CT)
}

projection_vagueCT_0.4 <- assign_projection(score_matVague,0.4)
projection_vagueCT_0.5 <- assign_projection(score_matVague,0.5)
projection_vagueCT_0.6 <- assign_projection(score_matVague,0.6)
projection_vagueCT_0.7 <- assign_projection(score_matVague,0.7)
projection_vagueCT_0.75 <- assign_projection(score_matVague,0.75)

projection_CT_0.4 <- assign_projection(score_mat,0.4)
projection_CT_0.5 <- assign_projection(score_mat,0.5)
projection_CT_0.6 <- assign_projection(score_mat,0.6)
projection_CT_0.7 <- assign_projection(score_mat,0.7)
projection_CT_0.75 <- assign_projection(score_mat,0.75)

cbdataHAR2_projection <- cbdataHAR2
cbdataHAR2_projection$projection_vagueCT <- projection_vagueCT_0.75
cbdataHAR2_projection$projection_CT <- projection_CT_0.75


pdf(file="cbdata_HAR_projectionCT_UMAP_features.pdf",width=10,height=4)
p1 <- DimPlot(cbdataHAR2_projection, label = TRUE,group.by="projection_vagueCT") + ggtitle("vagueCT")
p2 <- DimPlot(cbdataHAR2_projection, label = TRUE,group.by="projection_CT") + ggtitle("CT")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1+p2
dev.off()

cbdataHAR2$projection <- projection_vagueCT_0.75




### tmp summary
p0 <- DimPlot(cbdataHAR2, reduction="umap",label=TRUE,group.by="sample") + ggtitle("sample")
p1 <- DimPlot(cbdataHAR2, reduction="umap",label=TRUE,group.by="cluster") + ggtitle("cluster")
p2 <- DimPlot(cbdataHAR2, reduction="umap",label=TRUE,group.by="projection") + ggtitle("projection CT")
pdf(file="cbdata_UMAP_sample_summary.pdf",width=5,height=4)
p0
dev.off()

pdf(file="cbdata_UMAP_cluster_summary.pdf",width=5,height=4)
p1
dev.off()

pdf(file="cbdata_UMAP_projection_summary.pdf",width=5,height=4)
p2
dev.off()

### key gene UMAP
# CD8
CD8_glist <- c("Cd3d","Cd8a","Cd8b1","Trbc1","Trbc2","Sell","Cd44","Cd69","Itgae")
CD8v2_glist <- c("Ms4a4b","Cd8b1","Ccl5","Il7r","Ms4a6b","Cd3d","Trbc2","Nkg7")
NK_glist <- c("Klrd1","Klra7","Klrc1","Klrk1")
CD4_glist <- c("Cd3d","Cd4","Trdc","Trbc1","Trbc2","Cd44","Tbx21","Ifng","Gata3","Stat4","Il4","Il5",
	           "Rorc","Rora","Il17a","Foxp3","Bcl6","Pdcd1","Izumo1r","Cxcr5")
B_glist <- c("Ptprc","Cd19","Apex1","Mzb1","Ly6a","Prdm1","Irf4","Ighm","Ighd")
Mono_glist <- c("Spn","Nr4a1","Cx3cr1","Cd36","Cd300e","Cd300a")
AM_glist <- c("Fabp4","Pparg","Siglecf","Ccl6","Atp6v0d2","Mrc1","Fabp1","Ctsd","Abcg1") #alveolar macrophages
IM_glist <- c("C1qb","C1qc","C1qa","Ms4a7","Pf4","Aif1","Ccl12","Cd163","Mgl2") #interstitial
DC_glist <- c("Irf4","Ccr7","Itgae","Xcr1","Zbtb46")
epi_glist <- c("Epcam","Sftpc","Sftpd","Pecam1")
Krt8_glist <- c("Krt8","Cldn4","Fosl1","Areg","Lgals3","Clu","Krt18","Sprr1a","Krt7","Cryab","Krt19")
Krt5_glist <- c("Krt5","Trp63","Krt17","Krt14","Aqp3","Ngfr","Ly6d")
AT1_glist <- c("Ager","Akap5","Clic5","Rtkn2","Vegfa","Cldn18","Emp2","Spock2","Timp3","Lmo7","Col4a3","Lgfbp2","Hopx")
AT2_glist <- c("Sftpa1","Slc34a2","Sftpb","Hc","Chi3l1","Lamp3","Sfta2","S100g","Sftpd")
activatedAT2_glist <- c("Lcn2","Sftpd","Chi3l1","Scd1","Cxcl15","Lyz1","Ager","Napsa","Lrg1","Cystm1","Rgcc")
summary_glist <- c("Cldn18","Clic5","Emp2","Vegfa","Sftpc","Lamp3","Hc","S100g","Slc34a2","Krt14","Aqp3","Krt17","Krt5","Cldn4","Krt8","Krt19","Cd8a","Cd8b1","Cd103","Cd3","C1qa","C1qb","Ms4a7","Cx3cr1","Il1b","Il6","Ndrg1","Ifng","Tnf","Lum","Sma","Col1a1","Scgb3a2","scgb1a1","Cc10","Foxj1","Pecam1","Plvap","Gpihbp1","Vwf","Trp63","Cd104","Aqp5","Krt15","Csf1","Ccl2","Csf1r")
IL1bSig_glist <- c("Ikbkb","Il1a","Il1b","Il1r1","Il1rap","Il6","Irak1","Irak2","Irak3","Irak4","Mapk3","Myd88")
inflammasome_glist <- c("Casp1","Ddx3x","Dhx33","Gsdmd","Aim2","Nlrp3","Pycard","Nlrp1a","Nlrp1b")
posRegTCR_glist <- c("Ada","Bcl10","Card11","Ccr7","Cd81","Cd226","Cyld","Ikbkg","Kcnn4","Lipa","Nectin2","Prkd2","Rab29","Rela","Rps3","Tespa1","Trat1","Usp12","Usp46")
inflamComplex_glist <- c("Aim2","Casp1","Casp4","Ddx3x","Dhx33","Gsdmd","Mefv","Naip1","Naip2","Naip5","Naip6","Naip7","Nlrc4","Nlrp1a","Nlrp1b","Nlrp3","Nlrp6","Nlrp9a","Nlrp9b","Nlrp9c","Pycard")
AE_glist <- c("Cldn18","Clic5","Emp2","Vegfa",
			  "Sftpc","Lamp3","Hc","S100g","Slc34a2")
TGFb_glist <- c("Gzmk","Eomes","Nkg7","Tox","Ccl2","Ccr2","Cxcl13","Lta","Ltb","Hif1a","Ifng","Ifna","Ifnb1","Cxcl16")

Aberrant_basaloid_glist <- c("Epcam","Chd1","Vim","Fn1","Col1a1","Chd2","Tnc","Vcan","Pcp4","Cux2","Spink1","Prss2","Cpa6","Ctse","Mmp7","Mdk","Gdf15","Ptgs2","Slco2a1","Ephb2","Itgb8","Itgav","Itgb6","Tgfbl","Kcnn4","Kcnq5","Kcns3","Cdkn1a","Cdkn2a","Cdkn2b","Ccnd1","Ccnd2","Mdm2","Hmga2","Ptchd4","Ociad2","Tp63","Krt17","Lab3","Lamc2")
Profibrotic_glist <- c("Csf1","Adamts9","Cdkn1a","Timp3","Col4a1","Col4a2","Slco2a1","Plvap","Mmp14","Fbn1","Fn1","Col5a1","Col1a2","Col1a1","Col3a1","Sparc","Sparcl1","Col6a1","Cd74","Mfap5","Spp1","Chil1")
ADI_DATP_PATS_glist <- c("Anxa2","Bcl2l1","Itga3","Jun","Krt8","Krt19","F3","Klf6","Msn","Tapbp","Cdkn1a","Icam1","Timp3","Anxa1","Krt18","Tnfrsf1a","Sdc4","Cd63","Npc2","Myh9","Actg1","B2m","Jund")


cbdataHAR2_scale_allgenes <-  SCTransform(cbdataHAR2, assay = "Spatial", verbose = FALSE, ,variable.features.n = 19465)


#saveRDS(cbdataHAR2_scale_allgenes, file="cbdataHAR2_scale_allgenes.rds")
cbdataHAR2_scale_allgenes <- readRDS(file="cbdataHAR2_scale_allgenes.rds")

dim(cbdataHAR2@assays$SCT@data)
dim(cbdataHAR2_scale_allgenes@assays$SCT@scale.data)

pdf(file="UMAP_keygeneExp_scatter_Krt8test.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=Krt8_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_CD8.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=CD8_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_CD8v2.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=CD8v2_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_CD4.pdf",width=24,height=18)
FeaturePlot(cbdataHAR2, features=CD4_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_NK.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=NK_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_B.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=B_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_Mono.pdf",width=9,height=9)
FeaturePlot(cbdataHAR2, features=Mono_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_DC.pdf",width=9,height=9)
FeaturePlot(cbdataHAR2, features=DC_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_AM.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=AM_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_IM.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=IM_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_epi.pdf",width=9,height=9)
FeaturePlot(cbdataHAR2, features=epi_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_Krt8.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2, features=Krt8_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_Krt5.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=Krt5_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_AT1.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2, features=AT1_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_AT2.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2, features=AT2_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_activatedAT2.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2, features=activatedAT2_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_summaryList.pdf",width=15,height=25)
FeaturePlot(cbdataHAR2_scale_allgenes, features=summary_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_Il1bSigList.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=IL1bSig_glist, reduction="umap")
dev.off()


pdf(file="UMAP_keygeneExp_scatter_inflammasomeList.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=inflammasome_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_posRegTCRlist.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=posRegTCR_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_inflamComplexlist.pdf",width=12,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=inflamComplex_glist, reduction="umap")
dev.off()

pdf(file="UMAP_keygeneExp_scatter_TGFb.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=TGFb_glist, reduction="umap")
dev.off()


expmat <- cbdataHAR2@assays$SCT@data

keygeneALL <- intersect( unique(c(CD8_glist,
CD8v2_glist ,
NK_glist ,
CD4_glist   ,
B_glist ,
Mono_glist,
AM_glist  ,
IM_glist  ,
DC_glist  ,
epi_glist ,
Krt8_glist,
Krt5_glist,
AT1_glist,
AT2_glist,
activatedAT2_glist)),rownames(expmat))


aveExpMat <-  matrix(rep(0, length(keygeneALL)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL))
rownames(aveExpMat) <- keygeneALL
colnames(aveExpMat) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
medianExpMat <- aveExpMat
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat[keygeneALL,usecell]
	aveExpMat[, paste0("C",i)] <- apply(tempMat, 1,mean)
	medianExpMat[, paste0("C",i)] <- apply(tempMat, 1,median)
}


expmat_ownscale <- scale(cbdataHAR2@assays$SCT@data)

keygeneALL <- intersect(unique(c(CD8_glist,
CD8v2_glist ,
NK_glist ,
CD4_glist   ,
B_glist ,
Mono_glist,
AM_glist  ,
IM_glist  ,
DC_glist  ,
epi_glist ,
Krt8_glist,
Krt5_glist,
AT1_glist,
AT2_glist,
activatedAT2_glist)),rownames(expmat))


aveExpMat_ownscale <-  matrix(rep(0, length(keygeneALL)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL))
rownames(aveExpMat_ownscale) <- keygeneALL
colnames(aveExpMat_ownscale) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
medianExpMat_ownscale <- aveExpMat_ownscale
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_ownscale[keygeneALL,usecell]
	aveExpMat_ownscale[, paste0("C",i)] <- apply(tempMat, 1,mean)
	medianExpMat_ownscale[, paste0("C",i)] <- apply(tempMat, 1,median)
}

intersect(c("Ifna1"), rownames(expmat_scale))

expmat_scale <- cbdataHAR2_scale_allgenes@assays$SCT@scale.data
keygeneALL_scale <- intersect( unique(c(CD8_glist,
CD8v2_glist ,
NK_glist ,
CD4_glist   ,
B_glist ,
Mono_glist,
AM_glist  ,
IM_glist  ,
DC_glist  ,
epi_glist ,
Krt8_glist,
Krt5_glist,
AT1_glist,
AT2_glist,
activatedAT2_glist)),rownames(expmat_scale))

keygeneALL_scale_epi <- intersect( unique(c(
epi_glist ,
Krt8_glist,
Krt5_glist,
AT1_glist,
AT2_glist,
activatedAT2_glist)),rownames(expmat_scale))

keygeneALL_scale_summary <- intersect( unique(c(
summary_glist )),rownames(expmat_scale))


keygeneALL_scale_summaryIL1AM <- intersect( unique(c(
summary_glist,
IL1bSig_glist,
AM_glist )),rownames(expmat_scale))


keygeneALL_scale_summaryAM <- intersect( unique(c(
summary_glist,
AM_glist )),rownames(expmat_scale))


keygeneALL_scale_IL1 <- intersect( unique(c(
IL1bSig_glist)),rownames(expmat_scale))

aveExpMat_scale <-  matrix(rep(0, length(keygeneALL_scale)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL_scale))
rownames(aveExpMat_scale) <- keygeneALL_scale
colnames(aveExpMat_scale) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
medianExpMat_scale <- aveExpMat_scale
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_scale[keygeneALL_scale,usecell]
	aveExpMat_scale[, paste0("C",i)] <- apply(tempMat, 1,mean)
	medianExpMat_scale[, paste0("C",i)] <- apply(tempMat, 1,median)
}

aveExpMat_scale_summary <-  matrix(rep(0, length(keygeneALL_scale_summary)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL_scale_summary))
rownames(aveExpMat_scale_summary) <- keygeneALL_scale_summary
colnames(aveExpMat_scale_summary) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_scale[keygeneALL_scale_summary,usecell]
	aveExpMat_scale_summary[, paste0("C",i)] <- apply(tempMat, 1,mean)
}

aveExpMat_scale_summaryIL1AM <-  matrix(rep(0, length(keygeneALL_scale_summaryIL1AM)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL_scale_summaryIL1AM))
rownames(aveExpMat_scale_summaryIL1AM) <- keygeneALL_scale_summaryIL1AM
colnames(aveExpMat_scale_summaryIL1AM) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_scale[keygeneALL_scale_summaryIL1AM,usecell]
	aveExpMat_scale_summaryIL1AM[, paste0("C",i)] <- apply(tempMat, 1,mean)
}

aveExpMat_scale_summaryAM <-  matrix(rep(0, length(keygeneALL_scale_summaryAM)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL_scale_summaryAM))
rownames(aveExpMat_scale_summaryAM) <- keygeneALL_scale_summaryAM
colnames(aveExpMat_scale_summaryAM) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_scale[keygeneALL_scale_summaryAM,usecell]
	aveExpMat_scale_summaryAM[, paste0("C",i)] <- apply(tempMat, 1,mean)
}






usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
pdf(file="keygeneExp_cluster_heatmap.pdf",height=12)
#heatmap.2(aveExpMat, trace="none",Rowv=F, col=ColorRamp,main="aveExp")
#heatmap.2(medianExpMat, trace="none",Rowv=F, col=ColorRamp,main="medianExp")
#heatmap.2(aveExpMat_ownscale, trace="none",Rowv=F, col=ColorRamp,main="aveExp ownscale")
#heatmap.2(medianExpMat_ownscale, trace="none",Rowv=F, col=ColorRamp,main="medianExp ownscale")
heatmap.2(aveExpMat_scale, trace="none",Rowv=F, col=ColorRamp,main="aveExp scale",cexRow=0.4)
heatmap.2(aveExpMat_scale[keygeneALL_scale_epi,], trace="none",Rowv=F, col=ColorRamp,main="aveExp scale epiGene")
heatmap.2(aveExpMat_scale_summary, trace="none",Rowv=F, col=ColorRamp,main="aveExp scale summaryGene")
heatmap.2(aveExpMat_scale_summaryIL1AM, trace="none",Rowv=F, col=ColorRamp,main="aveExp scale summary+IL1+AM")
heatmap.2(aveExpMat_scale_summaryAM, trace="none",Rowv=F, col=ColorRamp,main="aveExp scale summary+AM")
dev.off()





aveExpMat_scale_IL1 <-  matrix(rep(0, length(keygeneALL_scale_IL1)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneALL_scale_IL1))
rownames(aveExpMat_scale_IL1) <- keygeneALL_scale_IL1
colnames(aveExpMat_scale_IL1) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_scale[keygeneALL_scale_IL1,usecell]
	aveExpMat_scale_IL1[, paste0("C",i)] <- apply(tempMat, 1,mean)
}

 ### 3 color ver
Imm_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("2","12","13"))]
Krt_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("16","4","7","11","21"))]
AE_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("0","5","14","18","3","10","19","22","25","23","1","8","20","26","17","6","9","15","24"))]

Imm_cluster <- c("2","12","13")
Krt_cluster <- c("16","4","7","11","21")
AE_cluster <- c("0","5","14","18","3","10","19","22","25","23","1","8","20","26","17","6","9","15","24")

usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
pdf(file="keygeneExp_IL1_celltype_heatmap.pdf",height=12)
heatmap.2(aveExpMat_scale_IL1, trace="none",Rowv=T, Colv=T, col=ColorRamp,main="clustering")
heatmap.2(aveExpMat_scale_IL1[,paste0("C",c(Imm_cluster, Krt_cluster, AE_cluster))], trace="none",Rowv=T, Colv=F, col=ColorRamp,main="celltype3 order")
dev.off()

pdf(file="keygeneExp_IL1_celltype_boxplot.pdf",height=12,width=16)
par(mfrow=c(3,4),mar=c(4,4,2,2))
for( thisGene in keygeneALL_scale_IL1){
	boxplot(expmat_scale[thisGene,Imm_cells ],
		    expmat_scale[thisGene,Krt_cells ],
		    expmat_scale[thisGene,AE_cells ], 
		    main=thisGene,colnames="normalized gene exp",
		    names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
}
dev.off()



## image analysis
S1A_highQ_SCT <- SCTransform(S1A_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
S2A_highQ_SCT <- SCTransform(S2A_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
S1D_highQ_SCT <- SCTransform(S1D_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)
S2D_highQ_SCT <- SCTransform(S2D_highQ, assay = "Spatial", verbose = FALSE,variable.features.n = 19465)

#saveRDS(S1A_highQ_SCT, file="S1A_highQ_SCT.rds")
#saveRDS(S1D_highQ_SCT, file="S1D_highQ_SCT.rds")
#saveRDS(S2A_highQ_SCT, file="S2A_highQ_SCT.rds")
#saveRDS(S2D_highQ_SCT, file="S2D_highQ_SCT.rds")
#
#S1A_highQ_SCT <- readRDS(file="S1A_highQ_SCT.rds")
#S1D_highQ_SCT <- readRDS(file="S1D_highQ_SCT.rds")
#S2A_highQ_SCT <- readRDS(file="S2A_highQ_SCT.rds")
#S2D_highQ_SCT <- readRDS(file="S2D_highQ_SCT.rds")

allgenes <- intersect(intersect(intersect(rownames(S1A_highQ_SCT@assays$SCT@scale.data), rownames(S2A_highQ_SCT@assays$SCT@scale.data)),rownames(S1D_highQ_SCT@assays$SCT@scale.data)), rownames(S2D_highQ_SCT@assays$SCT@scale.data))
markerGeneList_spotFeatures <- function(Glist, usename){
	dirname <- paste0("spotFeatures_keygenes/",usename)
	if (!file.exists(dirname)){
	    dir.create(dirname)
	}

	for(this_keygene in intersect(allgenes, unique(Glist)) ){
		print(this_keygene)
		pdf(file=paste0(dirname,"/spotFeature_scaleExp_",this_keygene,".pdf"),height=9,width=9)
		par(bg="grey",mfrow=c(2,2),mar=c(4,4,2,2))
		spotFeature_plot(S1A_highQ_SCT@images$slice1@coordinates, S1A_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S1A scaleExp ",this_keygene))
		spotFeature_plot(S1D_highQ_SCT@images$slice1@coordinates, S1D_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S1D scaleExp ",this_keygene))
		spotFeature_plot(S2A_highQ_SCT@images$slice1@coordinates, S2A_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S2A scaleExp ",this_keygene))
		spotFeature_plot(S2D_highQ_SCT@images$slice1@coordinates, S2D_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S2D scaleExp ",this_keygene))
		dev.off()
	}
}

markerGeneList_spotFeatures(c(CD8_glist,CD8v2_glist),"CD8")
markerGeneList_spotFeatures(c(NK_glist),"NK")
markerGeneList_spotFeatures(c(CD4_glist),"CD4")
markerGeneList_spotFeatures(c(B_glist),"B")
markerGeneList_spotFeatures(c(Mono_glist),"Mono")
markerGeneList_spotFeatures(c(AM_glist),"AM")
markerGeneList_spotFeatures(c(IM_glist),"IM")
markerGeneList_spotFeatures(c(DC_glist),"DC")
markerGeneList_spotFeatures(c(epi_glist),"epi")
markerGeneList_spotFeatures(c(Krt8_glist),"Krt8")
markerGeneList_spotFeatures(c(Krt5_glist),"Krt5")
markerGeneList_spotFeatures(c(AT1_glist),"AT1")
markerGeneList_spotFeatures(c(AT2_glist),"AT2")
markerGeneList_spotFeatures(c(activatedAT2_glist),"activatedAT2")
markerGeneList_spotFeatures(c(summary_glist),"summary")
markerGeneList_spotFeatures(c(inflammasome_glist),"inflammasome")
markerGeneList_spotFeatures(c(IL1bSig_glist),"IL1bSig")
markerGeneList_spotFeatures(c(posRegTCR_glist),"posRegTCR")
markerGeneList_spotFeatures(c(inflamComplex_glist),"inflamComplex")
markerGeneList_spotFeatures(c(TGFb_glist),"TGFb")

markerGeneList_spotFeatures(c("Mertk"),"summary")
markerGeneList_spotFeatures(c("Krt17"),"summary")
markerGeneList_spotFeatures(c("Krt8"),"Krt8")

dirname <- "spotFeatures_keygenes/Krt5"
this_keygene <- "Krt17"
		pdf(file=paste0(dirname,"/spotFeature_scaleExp_",this_keygene,".pdf"),height=9,width=9)
		par(bg="grey",mfrow=c(2,2),mar=c(4,4,2,2))
		spotFeature_plot(S1A_highQ_SCT@images$slice1@coordinates, S1A_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S1A scaleExp ",this_keygene))
		plot(1,1,type="n",axes=F,xlab="",ylab="")
#		spotFeature_plot(S1D_highQ_SCT@images$slice1@coordinates, S1D_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S1D scaleExp ",this_keygene))
		spotFeature_plot(S2A_highQ_SCT@images$slice1@coordinates, S2A_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S2A scaleExp ",this_keygene))
		plot(1,1,type="n",axes=F,xlab="",ylab="")
#		spotFeature_plot(S2D_highQ_SCT@images$slice1@coordinates, S2D_highQ_SCT@assays$SCT@scale.data[this_keygene,], usecolor = c("blue","white","red"), log_scale = 0, M=paste0("S2D scaleExp ",this_keygene))
		dev.off()


#intersect(rownames(cbdataHAR2),c("Cd103", "Cd3", "Sma", "scgb1a1", "Cc10", "Cd104"))
#epi_cells: C4,C7,C11
#Imm_cells: C2,  C12, C13

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

#### 6 color ver
Imm_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("2","12","13"))]
Krt8_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("16"))]
Krt58_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("4","7","11","21"))]
AT1_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("0","5","14","18"))]
AT2_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("3","10","19","22","25"))]
AM_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("23"))]

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

S1A_Imm_cells <- fetch_sample_cell(Imm_cells, "S1A")
S2A_Imm_cells <- fetch_sample_cell(Imm_cells, "S2A")
S1D_Imm_cells <- fetch_sample_cell(Imm_cells, "S1D")
S2D_Imm_cells <- fetch_sample_cell(Imm_cells, "S2D")

S1A_Krt8_cells <- fetch_sample_cell(Krt8_cells, "S1A")
S2A_Krt8_cells <- fetch_sample_cell(Krt8_cells, "S2A")
S1D_Krt8_cells <- fetch_sample_cell(Krt8_cells, "S1D")
S2D_Krt8_cells <- fetch_sample_cell(Krt8_cells, "S2D")

S1A_Krt58_cells <- fetch_sample_cell(Krt58_cells, "S1A")
S2A_Krt58_cells <- fetch_sample_cell(Krt58_cells, "S2A")
S1D_Krt58_cells <- fetch_sample_cell(Krt58_cells, "S1D")
S2D_Krt58_cells <- fetch_sample_cell(Krt58_cells, "S2D")

S1A_AT1_cells <- fetch_sample_cell(AT1_cells, "S1A")
S2A_AT1_cells <- fetch_sample_cell(AT1_cells, "S2A")
S1D_AT1_cells <- fetch_sample_cell(AT1_cells, "S1D")
S2D_AT1_cells <- fetch_sample_cell(AT1_cells, "S2D")

S1A_AT2_cells <- fetch_sample_cell(AT2_cells, "S1A")
S2A_AT2_cells <- fetch_sample_cell(AT2_cells, "S2A")
S1D_AT2_cells <- fetch_sample_cell(AT2_cells, "S1D")
S2D_AT2_cells <- fetch_sample_cell(AT2_cells, "S2D")

S1A_AM_cells <- fetch_sample_cell(AM_cells, "S1A")
S2A_AM_cells <- fetch_sample_cell(AM_cells, "S2A")
S1D_AM_cells <- fetch_sample_cell(AM_cells, "S1D")
S2D_AM_cells <- fetch_sample_cell(AM_cells, "S2D")


Color_spotFeature <- function(Coordinate, g1cell, g2cell,g3cell,g4cell, g5cell,g6cell, M){
	XLYL <- c(0,130)
	g1_cell_names <- g1cell#unlist(lapply(g1cell, trimcellname))
	g2_cell_names <- g2cell#unlist(lapply(g2cell, trimcellname))	
	g3_cell_names <- g3cell#unlist(lapply(g3cell, trimcellname))	
	g4_cell_names <- g4cell#unlist(lapply(g1cell, trimcellname))
	g5_cell_names <- g5cell#unlist(lapply(g2cell, trimcellname))	
	g6_cell_names <- g6cell#unlist(lapply(g3cell, trimcellname))	
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[g1_cell_names] <- "red"
	color_feature[g2_cell_names] <- "blue"
	color_feature[g3_cell_names] <- "#AAAA00"
	color_feature[g4_cell_names] <- "darkgreen"
	color_feature[g5_cell_names] <- "purple"
	color_feature[g6_cell_names] <- "black"
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}


pdf(file=paste0("spotFeature_ImmKrtAT_cluster.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Color_spotFeature(S1A_highQ_SCT@images$slice1@coordinates, S1A_Imm_cells, S1A_Krt8_cells, S1A_Krt58_cells, S1A_AT1_cells, S1A_AT2_cells, S1A_AM_cells, M="S1A")
Color_spotFeature(S2A_highQ_SCT@images$slice1@coordinates, S2A_Imm_cells, S2A_Krt8_cells, S2A_Krt58_cells, S2A_AT1_cells, S2A_AT2_cells, S2A_AM_cells, M="S2A")
Color_spotFeature(S1D_highQ_SCT@images$slice1@coordinates, S1D_Imm_cells, S1D_Krt8_cells, S1D_Krt58_cells, S1D_AT1_cells, S1D_AT2_cells, S1D_AM_cells, M="S1D")
Color_spotFeature(S2D_highQ_SCT@images$slice1@coordinates, S2D_Imm_cells, S2D_Krt8_cells, S2D_Krt58_cells, S2D_AT1_cells, S2D_AT2_cells, S2D_AM_cells, M="S2D")
dev.off()


Imm_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("2","12","13"))]
Krt8_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("16"))]
Krt58_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("4","7","11","21"))]
AT1_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("0","5","14","18"))]
AT2_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("3","10","19","22","25"))]
AM_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("23"))]

celltype6 <- rep("unknown",ncol(cbdataHAR2))
names(celltype6) <- colnames(cbdataHAR2)
celltype6[Imm_cells] <- "Imm"
celltype6[Krt8_cells] <- "Krt8"
celltype6[Krt58_cells] <- "Krt5&8"
celltype6[AT1_cells] <- "AT1"
celltype6[AT2_cells] <- "AT2"
celltype6[AM_cells] <- "AM"

cbdataHAR2$celltype6 <- celltype6



 ### 3 color ver
Imm_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("2","12","13"))]
Krt_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("16","4","7","11","21"))]
AE_cells <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("0","5","14","18","3","10","19","22","25","23","1","8","20","26","17","6","9","15","24"))]


S1A_Imm_cells <- fetch_sample_cell(Imm_cells, "S1A")
S2A_Imm_cells <- fetch_sample_cell(Imm_cells, "S2A")
S1D_Imm_cells <- fetch_sample_cell(Imm_cells, "S1D")
S2D_Imm_cells <- fetch_sample_cell(Imm_cells, "S2D")

S1A_Krt_cells <- fetch_sample_cell(Krt_cells, "S1A")
S2A_Krt_cells <- fetch_sample_cell(Krt_cells, "S2A")
S1D_Krt_cells <- fetch_sample_cell(Krt_cells, "S1D")
S2D_Krt_cells <- fetch_sample_cell(Krt_cells, "S2D")

S1A_AE_cells <- fetch_sample_cell(AE_cells, "S1A")
S2A_AE_cells <- fetch_sample_cell(AE_cells, "S2A")
S1D_AE_cells <- fetch_sample_cell(AE_cells, "S1D")
S2D_AE_cells <- fetch_sample_cell(AE_cells, "S2D")


Color3_spotFeature <- function(Coordinate, g1cell, g2cell,g3cell, M){
	XLYL <- c(0,130)
	g1_cell_names <- g1cell#unlist(lapply(g1cell, trimcellname))
	g2_cell_names <- g2cell#unlist(lapply(g2cell, trimcellname))	
	g3_cell_names <- g3cell#unlist(lapply(g3cell, trimcellname))	
	spot_names <- rownames(Coordinate)
	x_pos <- Coordinate[spot_names,"row"]
	y_pos <- Coordinate[spot_names,"col"]
	color_feature <- rep("white", length(spot_names))
	names(color_feature) <- spot_names
	color_feature[g1_cell_names] <- "red"
	color_feature[g2_cell_names] <- "blue"
	color_feature[g3_cell_names] <- "#AAAA00"
	plot(x_pos,y_pos,col=color_feature,pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL, main=M, xlab="", ylab="")
	#return(list(x_pos, y_pos, color_feature))
}



pdf(file=paste0("spotFeature_ImmKrtAE3color_cluster.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Color3_spotFeature(S1A_highQ_SCT@images$slice1@coordinates, S1A_Imm_cells, S1A_Krt_cells, S1A_AE_cells, M="S1A")
Color3_spotFeature(S2A_highQ_SCT@images$slice1@coordinates, S2A_Imm_cells, S2A_Krt_cells, S2A_AE_cells, M="S2A")
Color3_spotFeature(S1D_highQ_SCT@images$slice1@coordinates, S1D_Imm_cells, S1D_Krt_cells, S1D_AE_cells, M="S1D")
Color3_spotFeature(S2D_highQ_SCT@images$slice1@coordinates, S2D_Imm_cells, S2D_Krt_cells, S2D_AE_cells, M="S2D")
dev.off()


### UMAP
celltype3 <- rep("Imm",ncol(cbdataHAR2))
names(celltype3) <- colnames(cbdataHAR2)
celltype3[Krt_cells] <- "Krt"
celltype3[AE_cells] <- "Epi"
cbdataHAR2$celltype3 <- celltype3


pdf(file="celltype3_UMAP_features.pdf",width=10,height=4)
p1 <- DimPlot(cbdataHAR2, label = TRUE,group.by="celltype3") + ggtitle("cell type3")
p2 <- DimPlot(cbdataHAR2, label = TRUE,group.by="celltype6") + ggtitle("cell type6")
p1 + p2
dev.off()


pdf(file="celltype3_UMAP_features_v2.pdf",width=5,height=4)
p1 <- DimPlot(cbdataHAR2, label = TRUE,group.by="celltype3") + ggtitle("cell type3")
#p2 <- DimPlot(cbdataHAR2, label = TRUE,group.by="celltype6") + ggtitle("cell type6")
p1 
dev.off()



#### slide & celltype gene boxplot
pdf(file="keygeneExp_CD8_SlideCelltype_boxplot.pdf",height=6,width=30)
par(mfcol=c(2,10),mar=c(4,4,2,2))
#thisGene <- "Cd8a"
for(thisGene in c(CD8_glist,"Ccl5")){
	boxplot(expmat_scale[thisGene,paste0("S1A_",S1A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_AE_cells) ], 
			main=paste0(thisGene," S1A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S2A_",S2A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_AE_cells )], 
			main=paste0(thisGene," S2A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )	
}
dev.off()


pdf(file="keygeneExp_IL1bSig_SlideCelltype_boxplot.pdf",height=9,width=36)
par(mfcol=c(3,12),mar=c(4,4,2,2))
#thisGene <- "Cd8a"
for(thisGene in c(IL1bSig_glist)){
	boxplot(expmat_scale[thisGene,Imm_cells ],
			expmat_scale[thisGene,Krt_cells ],
			expmat_scale[thisGene,AE_cells  ], 
			main=paste0(thisGene," combine"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S1A_",S1A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_AE_cells) ], 
			main=paste0(thisGene," S1A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S2A_",S2A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_AE_cells )], 
			main=paste0(thisGene," S2A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )	
}
dev.off()


pdf(file="keygeneExp_inflammasome_SlideCelltype_boxplot.pdf",height=9,width=24)
par(mfcol=c(3,8),mar=c(4,4,2,2))
#thisGene <- "Cd8a"
for(thisGene in intersect(rownames(expmat_scale),inflammasome_glist)){
	boxplot(expmat_scale[thisGene,Imm_cells ],
			expmat_scale[thisGene,Krt_cells ],
			expmat_scale[thisGene,AE_cells  ], 
			main=paste0(thisGene," combine"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S1A_",S1A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_AE_cells) ], 
			main=paste0(thisGene," S1A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S2A_",S2A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_AE_cells )], 
			main=paste0(thisGene," S2A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )	
}
dev.off()


pdf(file="keygeneExp_TCRSig_SlideCelltype_boxplot.pdf",height=9,width=18*3)
par(mfcol=c(3,18),mar=c(4,4,2,2))
#thisGene <- "Cd8a"
for(thisGene in intersect(rownames(expmat_scale),posRegTCR_glist)){
	boxplot(expmat_scale[thisGene,Imm_cells ],
			expmat_scale[thisGene,Krt_cells ],
			expmat_scale[thisGene,AE_cells  ], 
			main=paste0(thisGene," combine"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S1A_",S1A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_AE_cells) ], 
			main=paste0(thisGene," S1A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S2A_",S2A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_AE_cells )], 
			main=paste0(thisGene," S2A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )	
}
dev.off()



pdf(file="keygeneExp_TCRSig_SlideCelltype_boxplot_v2.pdf",height=9,width=18*3)
par(mfcol=c(3,18),mar=c(4,4,2,2))
#thisGene <- "Cd8a"
for(thisGene in intersect(rownames(expmat_scale),posRegTCR_glist)){
	boxplot(expmat_scale[thisGene,Imm_cells ],
			expmat_scale[thisGene,Krt_cells ],
			expmat_scale[thisGene,AE_cells  ], 
			main=paste0(thisGene," combine"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S1A_",S1A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S1A_",S1A_AE_cells) ], 
			main=paste0(thisGene," S1A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	boxplot(expmat_scale[thisGene,paste0("S2A_",S2A_Imm_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_Krt_cells) ],
			expmat_scale[thisGene,paste0("S2A_",S2A_AE_cells )], 
			main=paste0(thisGene," S2A"),colnames="normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )	
}
dev.off()


#c(CD8_glist,"Ccl5")
#inflammasome_glist
#IL1bSig_glist
#posRegTCR_glist
#Aberrant_basaloid_glist
#Profibrotic_glist
#ADI_DATP_PATS_glist 
gene_lists <- list(CD8new_glist, 
					inflammasome_glist,
					IL1bSig_glist,
					posRegTCR_glist,
					Aberrant_basaloid_glist,
					Profibrotic_glist,
					ADI_DATP_PATS_glist
					)

sumdata <- c()

for(in_glist in gene_lists){
	thisGeneList <- intersect(rownames(expmat_scale),in_glist)
	AE_geneAveExp_this <- apply(expmat_scale[thisGeneList,AE_cells  ],1,mean)
	Krt_geneAveExp_this <- apply(expmat_scale[thisGeneList,Krt_cells ],1,mean)
	ave_AE <- mean(AE_geneAveExp_this)
	ave_Krt <- mean(Krt_geneAveExp_this)
	diffexp <- ave_Krt - ave_AE
	if(diffexp > 0){
		PV <- -log10(t.test(Krt_geneAveExp_this,AE_geneAveExp_this,alternative="greater")$p.val)
	}else{
		PV <- -log10(t.test(Krt_geneAveExp_this,AE_geneAveExp_this,alternative="less")$p.val)
	}
	sumdata <- rbind(sumdata, c(diffexp, PV))
}
colnames(sumdata) <- c("diffexp","PV")
rownames(sumdata) <- c("CD8","inflammasome","Il1bSig","posRegTCR","ABbasaloid","Profibrotic","ADI_DATP_PATS")



CT3_expBox_glist<- function(in_glist, inname){
	thisGeneList <- intersect(rownames(expmat_scale),in_glist)
	boxplot(apply(expmat_scale[thisGeneList,Imm_cells ],2,mean),
			apply(expmat_scale[thisGeneList,Krt_cells ],2,mean),
			apply(expmat_scale[thisGeneList,AE_cells  ],2,mean), 
			main=paste0(inname," genes"),ylab="mean normalized gene exp",
			names=c("Imm","Krt","AE"),col=c("red","blue","#AAAA00") )
	PV12 <- wilcox.test(apply(expmat_scale[thisGeneList,Imm_cells ],2,mean),
						apply(expmat_scale[thisGeneList,Krt_cells ],2,mean))$p.val
	PV13 <- wilcox.test(apply(expmat_scale[thisGeneList,Imm_cells ],2,mean),
						apply(expmat_scale[thisGeneList,AE_cells ],2,mean))$p.val
	PV23 <- wilcox.test(apply(expmat_scale[thisGeneList,Krt_cells ],2,mean),
						apply(expmat_scale[thisGeneList,AE_cells ],2,mean))$p.val
	legend("topleft",paste0(c("Imm_vs_Krt","Imm_vs_AE","Krt_vs_AE"),": P=", c(PV12,PV13,PV23)),
					 bty="n",pch=16,col=c("red","blue","#AAAA00") )
}



pdf(file="keygeneExp_genelistAVE_SlideCelltype_boxplot.pdf",height=8,width=16)
par(mfrow=c(2,4),mar=c(4,4,2,2))
CT3_expBox_glist(c(CD8_glist,"Ccl5"), "CD8")
CT3_expBox_glist(inflammasome_glist, "inflammasome")
CT3_expBox_glist(IL1bSig_glist, "IL1bSig")
CT3_expBox_glist(posRegTCR_glist, "posRegTCR")
CT3_expBox_glist(Aberrant_basaloid_glist, "Aberrant_basaloid")
CT3_expBox_glist(Profibrotic_glist, "Profibrotic")
CT3_expBox_glist(ADI_DATP_PATS_glist, "ADI_DATP_PATS")
dev.off()

CD8new_glist <- c(CD8_glist,"Ccl5")

# Assuming you have a matrix named gene_expression_matrix.
# rows should correspond to genes and columns to cells.
# Your gene lists and cell lists should be vectors with the names of the genes and cells.

# Assuming you have a matrix named gene_expression_matrix.
# rows should correspond to genes and columns to cells.
# Your gene lists and cell lists should be lists of vectors with the names of the genes and cells.
gene_lists <- list(CD8new_glist, 
					inflammasome_glist,
					IL1bSig_glist,
					posRegTCR_glist,
					Aberrant_basaloid_glist,
					Profibrotic_glist,
					ADI_DATP_PATS_glist
					)
set.seed(1)
cell_lists <- list(Imm_cells, Krt_cells,AE_cells)
#usegenes <- intersect(unlist(gene_lists),rownames(expmat_scale))
usegenes <-c(intersect(CD8new_glist,rownames(expmat_scale)),
			 intersect(inflammasome_glist,rownames(expmat_scale)),
			 intersect(IL1bSig_glist,rownames(expmat_scale)),
			 intersect(posRegTCR_glist,rownames(expmat_scale)),
			 intersect(Aberrant_basaloid_glist,rownames(expmat_scale)),
			 intersect(Profibrotic_glist,rownames(expmat_scale)),
			 intersect(ADI_DATP_PATS_glist,rownames(expmat_scale)))
usegeneNum <- c(length(intersect(CD8new_glist,rownames(expmat_scale))),
			 length(intersect(inflammasome_glist,rownames(expmat_scale))),
			 length(intersect(IL1bSig_glist,rownames(expmat_scale))),
			 length(intersect(posRegTCR_glist,rownames(expmat_scale))),
			 length(intersect(Aberrant_basaloid_glist,rownames(expmat_scale))),
			 length(intersect(Profibrotic_glist,rownames(expmat_scale))),#,
			 length(intersect(ADI_DATP_PATS_glist,rownames(expmat_scale))))
gene_expression_matrix<-expmat_scale[usegenes,unlist(cell_lists)]

row_annotation <- c(rep("CD8",length(intersect(CD8new_glist,rownames(expmat_scale)))),
					rep("inflammasome",length(intersect(inflammasome_glist,rownames(expmat_scale)))),
					rep("IL1bSig",length(intersect(IL1bSig_glist,rownames(expmat_scale)))),
					rep("posRegTCR",length(intersect(posRegTCR_glist,rownames(expmat_scale)))),
					rep("Aberrant_basaloid",length(intersect(Aberrant_basaloid_glist,rownames(expmat_scale)))),
					rep("Profibrotic",length(intersect(Profibrotic_glist,rownames(expmat_scale)))),
					rep("ADI_DATP_PATS",length(intersect(ADI_DATP_PATS_glist,rownames(expmat_scale)))))



col_annotation <- c(rep("Imm",length(Imm_cells)),
					rep("Krt",length(Krt_cells)),
					rep("AE",length(AE_cells)))


# Create annotation colors list for the rows and columns# Define colors for the cell groups
cell_group_colors <- c("Imm" = "red", "Krt" = "blue", "AE" = "#AAAA00")
rain <- rainbow(10)
gene_group_colors <- c("CD8"=rain[1],
					   "inflammasome"=rain[2],
					   "IL1bSig"=rain[3],
					   "posRegTCR"=rain[4],
					   "Aberrant_basaloid"=rain[5],
					   "Profibrotic"=rain[6],
					   "ADI_DATP_PATS"=rain[7]
							)
annotation_colors <- list(Gene_Group = gene_group_colors, Cell_Group = cell_group_colors)

gene_expression_matrix_trim<-gene_expression_matrix
gene_expression_matrix_trim[gene_expression_matrix_trim>5] <- 5
gene_expression_matrix_trim[gene_expression_matrix_trim< -5] <- -5
# Create the heatmap
pdf(file="keygeneExp_genelist_Celltyp3_heatmap.pdf",width=16,height=16)
par(mar=c(8,8,8,8))
pheatmap(gene_expression_matrix, breaks = seq(-5, 5, length.out = 101),
		 gaps_col=cumsum(c(length(Imm_cells),length(Krt_cells),length(AE_cells))),gaps_row=cumsum(usegeneNum),
         annotation_row = data.frame(Gene_Group=row_annotation),
         annotation_col = data.frame(Cell_Group=col_annotation),
         annotation_colors = annotation_colors,
         show_colnames = FALSE,
         show_rownames = FALSE,
         cluster_rows = FALSE, 
         cluster_cols = FALSE)
dev.off()







#### Celltype3 marker gene
#cbdataHAR2
#CT3markers <- FindAllMarkers(cbdataHAR2,group.by="cluster", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use="t")

#Clustermarkers <- FindAllMarkers(cbdataHAR2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DefaultAssay(cbdataHAR2) <- "SCT"
CT3markers_AT_vs_Krt <- FindMarkers(cbdataHAR2,group.by="celltype3", ident.1 = "Epi", ident.2 = "Krt")

CT3markers_AT_vs_Krt_defaultRNA <- CT3markers_AT_vs_Krt

pdf(file="UMAP_keygeneExp_scatter_AExKrtDEG.pdf",width=15,height=9)
FeaturePlot(cbdataHAR2_scale_allgenes, features=Krt8_glist, reduction="umap")
dev.off()


CT3markers_AT_g_Krt_highQ <- CT3markers_AT_vs_Krt[which(CT3markers_AT_vs_Krt[,"p_val_adj"] < 0.001 & (CT3markers_AT_vs_Krt[,"avg_log2FC"]) >= log2(1.5)),]
CT3markers_AT_l_Krt_highQ <- CT3markers_AT_vs_Krt[which(CT3markers_AT_vs_Krt[,"p_val_adj"] < 0.001 & (CT3markers_AT_vs_Krt[,"avg_log2FC"]) <= -log2(1.5)),]

write.table(row.names(CT3markers_AT_g_Krt_highQ),file="CT3markers_AT_g_Krt_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )
write.table(row.names(CT3markers_AT_l_Krt_highQ),file="CT3markers_AT_l_Krt_highQ_glist.txt",row.names=F,col.names=F,sep="\t",quote=F )

CT3markers_AT_g_Krt_highQ[order(CT3markers_AT_g_Krt_highQ[,2],decreasing=T),c(2,5)]
CT3markers_AT_l_Krt_highQ[order(CT3markers_AT_l_Krt_highQ[,2],decreasing=F),c(2,5)]

CT3marker_ATvsKrt_DEGhighQ <- rbind(CT3markers_AT_g_Krt_highQ, CT3markers_AT_l_Krt_highQ)
outdat<- CT3marker_ATvsKrt_DEGhighQ[order(CT3marker_ATvsKrt_DEGhighQ[,2]),]
write.table(outdat,file="CT3marker_ATvsKrt_DEGhighQ.txt",row.names=T,col.names=T,sep="\t",quote=F)


### add annotation, CT3
CT3 <- rep("AE",length(cbdataHAR2$orig.ident))
names(CT3) <- names(cbdataHAR2$orig.ident)
CT3[Imm_cells] <- "Imm"
CT3[Krt_cells] <- "Krt"
cbdataHAR2$CT3anno <- CT3

CT3 <- rep("AE",length(cbdataHAR2$orig.ident))
names(CT3) <- names(cbdataHAR2$orig.ident)
CT3[Imm_cells] <- "Imm"
CT3[Krt_cells] <- "Krt"
cbdataHAR2$CT3anno <- CT3



############### use CD8 and Krt8 expression as celltype score, CD8 for immune cell, krt8 for krt cell

rawMat <- cbdataHAR2@assays$SCT@data
scaleMat <- cbdataHAR2_scale_allgenes@assays$SCT@scale.data

### old ver , CD8a; Krt8
ImmExp <- scaleMat[c("Cd8a",),]
KrtExp <- scaleMat[c("Krt8",),]

pdf(file="ImmKrtSubgroup1gene_keyGeneExp_Hist.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2))
#hist(rawMat["Cd8a",],n=200, main="Cd8a",xlab="normalized exp")
#hist(rawMat["Krt8",],n=200, main="Krt8",xlab="normalized exp")
hist(ImmExp,n=200, main="Cd8a",xlab="SCT normalized exp",xlim=c(-2,5))
abline(v=0,col="red",lwd=2,lty=2)
hist(KrtExp,n=200, main="Krt8",xlab="SCT normalized exp",xlim=c(-2,5))
abline(v=-0.6,col="red",lwd=2,lty=2)
dev.off()












######### new ver, CD8a+; Krt8+Krt5+Krt17
ImmExp <- apply(scaleMat[c("Cd8a","Cx3cr1","Mertk"),],2,mean)
KrtExp <- apply(scaleMat[c("Krt8","Krt5","Krt17"),],2,mean)
pdf(file="ImmKrtSubgroup_keyGeneExp_Hist.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2))
#hist(rawMat["Cd8a",],n=200, main="Cd8a",xlab="normalized exp")
#hist(rawMat["Krt8",],n=200, main="Krt8",xlab="normalized exp")
hist(ImmExp,n=200, main="mean(Cd8a,Cx3cr1,Mertk)",xlab="SCT normalized exp",xlim=c(-2,5))
abline(v=-0.4,col="red",lwd=2,lty=2)
hist(KrtExp,n=200, main="mean(Krt8,Krt5,Krt17)",xlab="SCT normalized exp",xlim=c(-2,5))
abline(v=-0.4,col="red",lwd=2,lty=2)
dev.off()


Cd8on_spots <- colnames(scaleMat)[which(ImmExp > -0.4)]
Krt8on_spots <- colnames(scaleMat)[which(KrtExp > -0.4)]


### AE gene ave heatmap
AE_glist <- c("Cldn18","Clic5","Emp2","Vegfa",
			  "Sftpc","Lamp3","Hc","S100g","Slc34a2")

expmat_scale <- cbdataHAR2_scale_allgenes@assays$SCT@scale.data
keygeneAE_scale <- intersect( unique(c(AE_glist)),rownames(expmat_scale))

aveExpMat_scale_AE <-  matrix(rep(0, length(keygeneAE_scale)*length(unique(cbdataHAR2$cluster)) ),nrow=length(keygeneAE_scale))
rownames(aveExpMat_scale_AE) <- keygeneAE_scale
colnames(aveExpMat_scale_AE) <- paste0("C",0:(length(unique(cbdataHAR2$cluster)) -1))
medianExpMat_scale_AE <- aveExpMat_scale_AE
for(i in 0:(length(unique(cbdataHAR2$cluster)) -1)){
	usecell <- names(cbdataHAR2$cluster)[which(cbdataHAR2$cluster==i)]
	tempMat <- expmat_scale[keygeneAE_scale,usecell]
	aveExpMat_scale_AE[, paste0("C",i)] <- apply(tempMat, 1,mean)
	medianExpMat_scale_AE[, paste0("C",i)] <- apply(tempMat, 1,median)
}

usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
pdf(file="AEkeygeneExp_cluster_heatmap.pdf",width=12)
heatmap.2(t(aveExpMat_scale_AE), trace="none",Colv=F, col=ColorRamp,main="aveExp scale AEgene")
dev.off()


clusterAVE_AEgeneExp <- apply(aveExpMat_scale_AE,2,mean)
pdf(file="AEkeygeneExp_cluster_aveBarplot.pdf")
barplot(sort(clusterAVE_AEgeneExp),las=2,ylab="mean normalized expression of AE genes")
dev.off()


### here remove C15, C24 from AE cells 
AE_cells_new <- colnames(cbdataHAR2)[which(cbdataHAR2$cluster %in% c("0","5","14","18","3","10","19","22","23","1","8","20","26","17","6","9","24"))]

Cd8on_S1A <- fetch_sample_cell(Cd8on_spots, "S1A")
Cd8on_S1D <- fetch_sample_cell(Cd8on_spots, "S1D")
Cd8on_S2A <- fetch_sample_cell(Cd8on_spots, "S2A")
Cd8on_S2D <- fetch_sample_cell(Cd8on_spots, "S2D")
Krt8on_S1A <- fetch_sample_cell(Krt8on_spots, "S1A")
Krt8on_S1D <- fetch_sample_cell(Krt8on_spots, "S1D")
Krt8on_S2A <- fetch_sample_cell(Krt8on_spots, "S2A")
Krt8on_S2D <- fetch_sample_cell(Krt8on_spots, "S2D")


S1A_AE_cells_new <- fetch_sample_cell(AE_cells_new, "S1A")
S2A_AE_cells_new <- fetch_sample_cell(AE_cells_new, "S2A")
S1D_AE_cells_new <- fetch_sample_cell(AE_cells_new, "S1D")
S2D_AE_cells_new <- fetch_sample_cell(AE_cells_new, "S2D")

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




#pdf(file=paste0("spotFeature_ImmKrtcolor_new_cmp.pdf"),height=9,width=9)
pdf(file=paste0("spotFeature_ImmKrtcolor_new.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(S1A_highQ_SCT@images$slice1@coordinates, Cd8on_S1A, Krt8on_S1A, S1A_AE_cells_new, M="S1A")
Cd8Krt8Mix_spotFeature(S2A_highQ_SCT@images$slice1@coordinates, Cd8on_S2A, Krt8on_S2A, S2A_AE_cells_new, M="S2A")
Cd8Krt8Mix_spotFeature(S1D_highQ_SCT@images$slice1@coordinates, Cd8on_S1D, Krt8on_S1D, S1D_AE_cells_new, M="S1D")
Cd8Krt8Mix_spotFeature(S2D_highQ_SCT@images$slice1@coordinates, Cd8on_S2D, Krt8on_S2D, S2D_AE_cells_new, M="S2D")
dev.off()


pdf(file=paste0("spotFeature_ImmcolorOnly_new.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(S1A_highQ_SCT@images$slice1@coordinates, Cd8on_S1A, c(), S1A_AE_cells_new, M="S1A")
Cd8Krt8Mix_spotFeature(S2A_highQ_SCT@images$slice1@coordinates, Cd8on_S2A, c(), S2A_AE_cells_new, M="S2A")
Cd8Krt8Mix_spotFeature(S1D_highQ_SCT@images$slice1@coordinates, Cd8on_S1D, c(), S1D_AE_cells_new, M="S1D")
Cd8Krt8Mix_spotFeature(S2D_highQ_SCT@images$slice1@coordinates, Cd8on_S2D, c(), S2D_AE_cells_new, M="S2D")
dev.off()


pdf(file=paste0("spotFeature_KrtcolorOnly_new.pdf"),height=9,width=9)
par(bg="grey",mfcol=c(2,2),mar=c(4,4,2,2))
Cd8Krt8Mix_spotFeature(S1A_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S1A,S1A_AE_cells_new, M="S1A")
Cd8Krt8Mix_spotFeature(S2A_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S2A,S2A_AE_cells_new, M="S2A")
Cd8Krt8Mix_spotFeature(S1D_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S1D,S1D_AE_cells_new, M="S1D")
Cd8Krt8Mix_spotFeature(S2D_highQ_SCT@images$slice1@coordinates,  c(), Krt8on_S2D,S2D_AE_cells_new, M="S2D")
dev.off()





### new annotation, CD8 ov Krt8

CTov <- rep("Cd8on",length(cbdataHAR2$orig.ident))
names(CTov) <- names(cbdataHAR2$orig.ident)
CTov[Krt8on_spots] <- "Krt8on"
CTov[intersect(Cd8on_spots, Krt8on_spots)] <- "CD8Krt8on"
CTov[AE_cells_new] <- "AE"
cbdataHAR2$CTOVanno <- CTov

### save data
#saveRDS(cbdataHAR2, file="cbdataHAR2.rds")
cbdataHAR2 <- readRDS(file="cbdataHAR2.rds")
DefaultAssay(cbdataHAR2) <- "SCT"


#cbdataHAR2$CTOVanno 
#CTOVmarkers <- FindAllMarkers(cbdataHAR2,group.by="CTOVanno", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#,test.use="t"
Cd8on_markers <- FindMarkers(object = cbdataHAR2,group.by="CTOVanno", ident.1 = "Cd8on", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
Krt8on_markers <- FindMarkers(object = cbdataHAR2,group.by="CTOVanno", ident.1 = "Krt8on", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
Cd8Krt8on_markers <- FindMarkers(object = cbdataHAR2,group.by="CTOVanno", ident.1 = "CD8Krt8on", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
AE_markers <- FindMarkers(object = cbdataHAR2,group.by="CTOVanno", ident.1 = "AE", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")

write.table(Cd8on_markers, file="Cd8on_markers_newMarkerGene.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(Krt8on_markers, file="Krt8on_markers_newMarkerGene.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(Cd8Krt8on_markers, file="Cd8Krt8on_markers_newMarkerGene.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(AE_markers, file="AE_markers_newMarkerGene.txt",row.names=T,col.names=T,sep="\t",quote=F)

Cd8on_markers_significant <- rownames(Cd8on_markers)[which(Cd8on_markers[,"avg_log2FC"]>=log2(1.5) & Cd8on_markers[,"p_val_adj"] < 0.05)]
Krt8on_markers_significant <- rownames(Krt8on_markers)[which(Krt8on_markers[,"avg_log2FC"]>=log2(1.5) & Krt8on_markers[,"p_val_adj"] < 0.05)]
Cd8Krt8on_markers_significant <- rownames(Cd8Krt8on_markers)[which(Cd8Krt8on_markers[,"avg_log2FC"]>=log2(1.5) & Cd8Krt8on_markers[,"p_val_adj"] < 0.05)]

write.table(Cd8on_markers_significant,file="Cd8on_markers_significant.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(Krt8on_markers_significant,file="Krt8on_markers_significant.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(Cd8Krt8on_markers_significant,file="Cd8Krt8on_markers_significant.txt",row.names=F,col.names=F,sep="\t",quote=F)

Cd8on_markers_FCorder <- Cd8on_markers[order(Cd8on_markers[,"avg_log2FC"],decreasing=T),]
Krt8on_markers_FCorder <- Krt8on_markers[order(Krt8on_markers[,"avg_log2FC"],decreasing=T),]
Cd8Krt8on_markers_FCorder <- Cd8Krt8on_markers[order(Cd8Krt8on_markers[,"avg_log2FC"],decreasing=T),]

write.table(cbind(rownames(Cd8on_markers_FCorder), Cd8on_markers_FCorder[,2]),file="Cd8on_markers_FCorder.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(cbind(rownames(Krt8on_markers_FCorder), Krt8on_markers_FCorder[,2]),file="Krt8on_markers_FCorder.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(cbind(rownames(Cd8Krt8on_markers_FCorder), Cd8Krt8on_markers_FCorder[,2]),file="Cd8Krt8on_markers_FCorder.txt",row.names=F,col.names=F,sep="\t",quote=F)


topgenes <- c(rownames(Cd8on_markers_FCorder)[1:20],
			  rownames(Krt8on_markers_FCorder)[1:20],
			  rownames(Cd8Krt8on_markers_FCorder)[1:20])



expmat_scale <- cbdataHAR2_scale_allgenes@assays$SCT@scale.data
keygeneAE_scale <- intersect( unique(c(AE_glist)),rownames(expmat_scale))

aveExpMat_topgenes <-  matrix(rep(0, length(topgenes)*length(unique(cbdataHAR2$CTOVanno)) ),nrow=length(topgenes))
rownames(aveExpMat_topgenes) <- topgenes
colnames(aveExpMat_topgenes) <- unique(cbdataHAR2$CTOVanno)
#medianExpMat_scale_AE <- aveExpMat_scale_AE
for(i in unique(cbdataHAR2$CTOVanno)) {
	usecell <- names(cbdataHAR2$CTOVanno)[which(cbdataHAR2$CTOVanno==i)]
	tempMat <- expmat_scale[topgenes,usecell]
	aveExpMat_topgenes[, i] <- apply(tempMat, 1,mean)
#	medianExpMat_scale_AE[, paste0("C",i)] <- apply(tempMat, 1,median)
}

usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
pdf(file="topgeneExp_CTov_heatmap.pdf",width=15)
heatmap.2(t(aveExpMat_topgenes), trace="none",Colv=F,Rowv=F, col=ColorRamp,main="aveExp scale topgenes")
dev.off()
















# seurat v3 batch removal 

cbdata <- merge(S1A_highQ, y=c(S1D_highQ,S2A_highQ,S2D_highQ),add.cell.ids=c("S1A","S1D","S2A","S2D"),project="spatial")

trimsamplename <- function(inname){
  return(strsplit(inname, "_")[[1]][1])
}
usesample <- unlist(lapply(names(cbdata$orig.ident), trimsamplename))
cbdata$sample <- usesample

#cbdata <- SCTransform(cbdata, assay = "Spatial", verbose = FALSE)

ifnb.list <- SplitObject(cbdata, split.by = "sample")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- SCTransform(x, assay="Spatial", verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)
cbdata.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
cbdata.combined <- IntegrateData(anchorset = cbdata.anchors)
DefaultAssay(cbdata.combined) <- "integrated"

cbdata.combined <- SCTransform(cbdata.combined, assay = "Spatial", verbose = FALSE)
#cbdata.combined <- ScaleData(cbdata.combined, verbose = FALSE)
cbdata.combined <- RunPCA(cbdata.combined, npcs = 30, verbose = FALSE, assay="SCT")
cbdata.combined <- RunUMAP(cbdata.combined, reduction = "pca", dims = 1:30)
cbdata.combined <- FindNeighbors(cbdata.combined, reduction = "pca", dims = 1:30)
cbdata.combined <- FindClusters(cbdata.combined, resolution = 0.5)

pdf(file="cbdata_SCT_sv3batchRemove_UMAP_features.pdf",width=8,height=4)
p1 <- DimPlot(cbdata.combined, label = TRUE,group.by="sample") + ggtitle("sample")
p2 <- DimPlot(cbdata.combined, reduction = "umap", group.by = "seurat_clusters")
#p3 <- DimPlot(cbdata, reduction = "umap", group.by = "predicted.id")
p1+p2
dev.off()


pdf(file="cbdata_spatial_features.pdf",width=16,height=4)
p1 <- SpatialDimPlot(cbdata, label = TRUE, label.size = 3)
p1
dev.off()






#################### human covid vs ctrl gene , exp here
allgenes <- rownames(cbdataHAR2_scale_allgenes@assays$SCT@data)
expTable <- cbdataHAR2_scale_allgenes@assays$SCT@data

mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
# separate human and mouse 
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
# remove some columns
mouse <- mouse[,c(1,4)]
human <- human[,c(1,4)]
# merge the 2 dataset  (note that the human list is longer than the mouse one)
mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 

covidHigh_filter <- read.table("/sfs/weka/scratch/sh8tv/lungST/humanST/humanLung_10x_spaceranger/human_condition_covidHigh_filter_DEGname.txt")[,1]
ctrlHigh_filter <- read.table("/sfs/weka/scratch/sh8tv/lungST/humanST/humanLung_10x_spaceranger/human_condition_ctrlHigh_filter_DEGname.txt")[,1]


covidHigh <- read.table("/sfs/weka/scratch/sh8tv/lungST/humanST/humanLung_10x_spaceranger/human_condition_covidHigh_DEGname.txt")[,1]
ctrlHigh <- read.table("/sfs/weka/scratch/sh8tv/lungST/humanST/humanLung_10x_spaceranger/human_condition_ctrlHigh_DEGname.txt")[,1]

covidHigh_mouseVer_raw <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% covidHigh),"Symbol.x"],allgenes)
ctrlHigh_mouseVer_raw <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% ctrlHigh),"Symbol.x"],allgenes)
covidHigh_mouseVer_filter <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% covidHigh_filter),"Symbol.x"],allgenes)
ctrlHigh_mouseVer_filter <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% ctrlHigh_filter),"Symbol.x"],allgenes)







Cd8cells <- names(cbdataHAR2_scale_allgenes$sample)[which(cbdataHAR2_scale_allgenes$sample %in% c("S1D","S2D"))]
IgGcells <- names(cbdataHAR2_scale_allgenes$sample)[which(cbdataHAR2_scale_allgenes$sample %in% c("S1A","S2A"))]

pdf(file="humanDEG_in_mouseCd8aIgG_box_raw.pdf")
par(mar=c(10,4,4,2),mfrow=c(1,2))
boxplot(apply(expTable[covidHigh_mouseVer_raw, IgGcells],1,mean),
				apply(expTable[covidHigh_mouseVer_raw, Cd8cells],1,mean),
				outline=F,ylab="aveExp",las=2,main='covidHigh defaultCutoff',
				names=c("IgG","Cd8a"),col=c("red","blue"))
boxplot(apply(expTable[ctrlHigh_mouseVer_raw, IgGcells],1,mean),
				apply(expTable[ctrlHigh_mouseVer_raw, Cd8cells],1,mean),
				outline=F,ylab="aveExp",las=2,main='ctrlHigh defaultCutoff',
				names=c("IgG","Cd8a"),col=c("red","blue"))
dev.off()


pdf(file="humanDEG_in_mouseCd8aIgG_box_filter.pdf")
par(mar=c(10,4,4,2),mfrow=c(1,2))
boxplot(apply(expTable[covidHigh_mouseVer_filter, IgGcells],1,mean),
				apply(expTable[covidHigh_mouseVer_filter, Cd8cells],1,mean),
				outline=F,ylab="aveExp",las=2,main='covidHigh filter',
				names=c("IgG","Cd8a"),col=c("red","blue"))
boxplot(apply(expTable[ctrlHigh_mouseVer_filter, IgGcells],1,mean),
				apply(expTable[ctrlHigh_mouseVer_filter, Cd8cells],1,mean),
				outline=F,ylab="aveExp",las=2,main='ctrlHigh filter',
				names=c("IgG","Cd8a"),col=c("red","blue"))
dev.off()

wilcox.test(apply(expTable[covidHigh_mouseVer_raw, IgGcells],1,mean),
						apply(expTable[covidHigh_mouseVer_raw, Cd8cells],1,mean),alternative="greater")
wilcox.test(apply(expTable[ctrlHigh_mouseVer_raw, IgGcells],1,mean),
						apply(expTable[ctrlHigh_mouseVer_raw, Cd8cells],1,mean),alternative="greater")


wilcox.test(apply(expTable[covidHigh_mouseVer_filter, IgGcells],1,mean),
						apply(expTable[covidHigh_mouseVer_filter, Cd8cells],1,mean),alternative="greater")
wilcox.test(apply(expTable[ctrlHigh_mouseVer_filter, IgGcells],1,mean),
						apply(expTable[ctrlHigh_mouseVer_filter, Cd8cells],1,mean),alternative="greater")




#### human covid > ctrl PATHWAY
get_clean_vector <- function(indata){
	my_vector <-as.vector(indata)
	clean_vector <- my_vector[nzchar(my_vector) & !is.na(my_vector)]
	return(clean_vector)
}

human_covidHigh_pathway_raw <- read.table("/sfs/weka/scratch/sh8tv/lungST/humanST/humanLung_10x_spaceranger/human_covid_g_ctrl_GSEA_geneLists.txt",row.names=1,sep="\t")


	cb_glist_all <- c()
pdf(file="covidHigh_pathwayGene_in_mouseCd8aIgG_box.pdf",width=3*12,height=5*2)
par(mar=c(10,4,6,2),mfcol=c(2,12))
for( i in 1:11){
	pathway_name <- rownames(human_covidHigh_pathway_raw)[i]
	glist_all <- get_clean_vector(human_covidHigh_pathway_raw[i,])
	cb_glist_all <- c(cb_glist_all, glist_all)
	glist_all_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% glist_all),"Symbol.x"],allgenes)
	glist_all_mouseVer_covidHigh <- intersect(glist_all_mouseVer, covidHigh_mouseVer_raw)
	print(c(pathway_name, length(glist_all_mouseVer),length(glist_all_mouseVer_covidHigh)))
	PV1 <- wilcox.test(apply(expTable[glist_all_mouseVer, IgGcells],1,mean),
										 apply(expTable[glist_all_mouseVer, Cd8cells],1,mean),alternative="greater")$p.val
	PV2 <- wilcox.test(apply(expTable[glist_all_mouseVer_covidHigh, IgGcells],1,mean),
										 apply(expTable[glist_all_mouseVer_covidHigh, Cd8cells],1,mean),alternative="greater")$p.val
	boxplot(apply(expTable[glist_all_mouseVer, IgGcells],1,mean),
					apply(expTable[glist_all_mouseVer, Cd8cells],1,mean),
					outline=F,ylab="aveExp all pathway genes",las=2,main=paste0(strsplit(pathway_name, "_")[[1]],collapse="\n"),
					names=c("IgG","Cd8a"),col=c("red","blue"),xlab=PV1)
	boxplot(apply(expTable[glist_all_mouseVer_covidHigh, IgGcells],1,mean),
					apply(expTable[glist_all_mouseVer_covidHigh, Cd8cells],1,mean),
					outline=F,ylab="aveExp pathway DEG",las=2,main=paste0(strsplit(pathway_name, "_")[[1]],collapse="\n"),
					names=c("IgG","Cd8a"),col=c("red","blue"),xlab=PV2)
}
	pathway_name <- "combine"
	glist_all_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% cb_glist_all),"Symbol.x"],allgenes)
	glist_all_mouseVer_covidHigh <- intersect(glist_all_mouseVer, covidHigh_mouseVer_raw)
	print(c(pathway_name, length(glist_all_mouseVer),length(glist_all_mouseVer_covidHigh)))
	PV1 <- wilcox.test(apply(expTable[glist_all_mouseVer, IgGcells],1,mean),
										 apply(expTable[glist_all_mouseVer, Cd8cells],1,mean),alternative="greater")$p.val
	PV2 <- wilcox.test(apply(expTable[glist_all_mouseVer_covidHigh, IgGcells],1,mean),
										 apply(expTable[glist_all_mouseVer_covidHigh, Cd8cells],1,mean),alternative="greater")$p.val
	boxplot(apply(expTable[glist_all_mouseVer, IgGcells],1,mean),
					apply(expTable[glist_all_mouseVer, Cd8cells],1,mean),
					outline=F,ylab="aveExp all pathway genes",las=2,main=paste0(strsplit(pathway_name, "_")[[1]],collapse="\n"),
					names=c("IgG","Cd8a"),col=c("red","blue"),xlab=PV1)
	boxplot(apply(expTable[glist_all_mouseVer_covidHigh, IgGcells],1,mean),
					apply(expTable[glist_all_mouseVer_covidHigh, Cd8cells],1,mean),
					outline=F,ylab="aveExp pathway DEG",las=2,main=paste0(strsplit(pathway_name, "_")[[1]],collapse="\n"),
					names=c("IgG","Cd8a"),col=c("red","blue"),xlab=PV2)
dev.off()




dev.off()






###### project RNA label
cbdataHAR2_prediction <- cbdataHAR2

allen_reference <- readRDS("/gpfs/gpfs0/project/zanglab_project/sh8tv/SunLab/Data/scRNA/AgedMouse_InfectionD60_Lung/01082023_Day61_combined_res2.rds")
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
scRNA_process <- SCTransform(allen_reference) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
DefaultAssay(scRNA_process) <- "SCT"

transfer_anchors <- FindTransferAnchors(
  reference = gse157079,
  query = sobj.combined,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:50
)


predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = gse157079$clusters,
  weight.reduction = sobj.combined[['pca']],
  dims = 1:50
sobj.combined <- AddMetaData(
  object = sobj.combined,
  metadata = predictions
))







projection_CT <- rep("NA",ncol(cbdataHAR2))


projection_assign_name <- function(score_mat){
	cell_max_score <- apply(score_mat, 2, max)
}



cbdataHAR2_prediction[["predictions_Cluster"]] <- predictions_Cluster
cbdataHAR2_prediction[["predictions_VagueCluster"]] <- predictions_VagueCluster


cbdataHAR2_prediction <- AddMetaData(object = cbdataHAR2_prediction, metadata = predictions_VagueCluster)




allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

#anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)














# "orig.ident" = original identity
DimPlot(brain.combined, group.by = "orig.ident")
# "ident" = identity, which are clusters
DimPlot(brain.combined, group.by = "ident", split.by = 'orig.ident')












p1 <- DimPlot(cbdata, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(cbdata, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
    5, 8)), facet.highlight = TRUE, ncol = 3)



de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
    selection.method = "markvariogram")









spot_names <- rownames(S1A@images$slice1@coordinates)
feature <- log10(S1A$nCount_Spatial[spot_names])
norm01_feature <- norm01(feature)
usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
color_feature <- ColorRamp[norm01_feature]



spot_names <- rownames(S1A@images$slice1@coordinates)
feature <- log10(S1A@assays$Spatial@counts["Pecam1",spot_names]+1)#log10(S1A$nCount_Spatial[spot_names])
norm01_feature <- norm01(feature)
usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
color_feature <- ColorRamp[norm01_feature]

x_pos <- S1A@images$slice1@coordinates[spot_names,"row"]
y_pos <- S1A@images$slice1@coordinates[spot_names,"col"]

pdf(file="pecam1_spatialExp.pdf")
par(bg="grey")
plot(x_pos,y_pos,col=color_feature,
	 pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL)
dev.off()


########################## example brain data test
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

pdf(file="brain_default_QCscatter.pdf")
SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")+scale_colour_manual(values = ColorRamp)
dev.off()

spot_names <- rownames(brain@images$anterior1@coordinates)
feature <- (brain$nCount_Spatial[spot_names])
norm01_feature <- norm01(feature)
usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
color_feature <- ColorRamp[norm01_feature]

a <- cbind( feature, norm01_feature,color_feature)

x_pos <- brain@images$anterior1@coordinates[spot_names,"row"]
y_pos <- brain@images$anterior1@coordinates[spot_names,"col"]

pdf(file="brain_test2.pdf")
par(bg="grey")
plot(x_pos,y_pos,col=color_feature,
	 pch=16, cex=0.5 ,xlim=XLYL, ylim=XLYL)
dev.off()


for(i in 1:nrow(S1A@images$slice1@coordinates)){
	spot_name <- rownames(S1A@images$slice1@coordinates)[i]
	x <- S1A@images$slice1@coordinates[spot_name,"row"]
	y <- S1A@images$slice1@coordinates[spot_name,"col"]
	nCount <- norm01_feature[i]
	useimage[x,y] <- nCount
}


pdf(file="default_QCscatter.pdf")
SpatialFeaturePlot(S1A, features = "nCount_Spatial") + theme(legend.position = "right")+scale_colour_manual(values = ColorRamp)
dev.off()

names(130)
plot(S1A@images$slice1@coordinates[,"col"], S1A@images$slice1@coordinates[,"row"], col="black",pch=16,xlim=XLYL,ylim=XLYL)
dev.off()











plot1 <- SpatialFeaturePlot(S1A, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

SpatialFeaturePlot(S1A, features = c("Cd8a", "Cd4", "Cd3d","Pecam1","Ptprc","Epcam"))
dev.off()