library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(harmony)
library(gplots)
library(pheatmap)

cbdataHAR2_scale_allgenes <- readRDS(file="/sfs/weka/scratch/sh8tv/lungST/mouseST/cbdataHAR2_scale_allgenes.rds")

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




human_covidHigh <- read.table("human_condition_covidHigh_DEGname.txt")[,1]
human_ctrlHigh <- read.table("human_condition_ctrlHigh_DEGname.txt")[,1]
human_KRTHigh <- read.table("DEG_KRT8_g_AGER_glist.txt")[,1]
human_AEHigh <- read.table("DEG_KRT8_l_AGER_glist.txt")[,1]
human_KRTHigh_covid <- read.table("DEG_KRT8_g_AGER_covid_glist.txt")[,1]
human_AEHigh_covid <- read.table("DEG_KRT8_l_AGER_covid_glist.txt")[,1]
mouse_KRTHigh <- read.table("mouse_celltype3_KrtHigh_DEGname.txt")[,1]
mouse_AEHigh <- read.table("mouse_celltype3_AEHigh_DEGname.txt")[,1]


human_covidHigh_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% human_covidHigh),"Symbol.x"],allgenes)
human_ctrlHigh_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% human_ctrlHigh),"Symbol.x"],allgenes)
human_KRTHigh_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% human_KRTHigh),"Symbol.x"],allgenes)
human_AEHigh_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% human_AEHigh),"Symbol.x"],allgenes)
human_KRTHigh_covid_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% human_KRTHigh_covid),"Symbol.x"],allgenes)
human_AEHigh_covid_mouseVer <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% human_AEHigh_covid),"Symbol.x"],allgenes)


hKrtHigh_mKrtHigh <- intersect(human_KRTHigh_mouseVer, mouse_KRTHigh)
hKrtHighCovid_mKrtHigh <- intersect(human_KRTHigh_covid_mouseVer, mouse_KRTHigh)
hCovidHigh_mKrtHigh <- intersect(human_covidHigh_mouseVer, mouse_KRTHigh)

hAeHigh_mAeHigh <- intersect(human_AEHigh_mouseVer, mouse_AEHigh)
hAEHighCovid_mAeHigh <- intersect(human_AEHigh_covid_mouseVer, mouse_AEHigh)
hCtrlHigh_mAeHigh <- intersect(human_ctrlHigh_mouseVer, mouse_AEHigh)

write.table(hKrtHigh_mKrtHigh, file="humanMouseKeyGeneSet/hKrtHigh_mKrtHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(hKrtHighCovid_mKrtHigh, file="humanMouseKeyGeneSet/hKrtHighCovid_mKrtHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(hCovidHigh_mKrtHigh, file="humanMouseKeyGeneSet/hCovidHigh_mKrtHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(hAeHigh_mAeHigh, file="humanMouseKeyGeneSet/hAeHigh_mAeHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(hAEHighCovid_mAeHigh, file="humanMouseKeyGeneSet/hAEHighCovid_mAeHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(hCtrlHigh_mAeHigh, file="humanMouseKeyGeneSet/hCtrlHigh_mAeHigh.txt",row.names=F,col.names=F,sep="\t",quote=F)
#Venn(human_KRTHigh_mouseVer,mouse_KRTHigh,"hKRT_mKRT")
#Venn(human_KRTHigh_covid_mouseVer,mouse_KRTHigh,"hKRTcovid_mKRT")
#Venn(human_covidHigh_mouseVer,mouse_KRTHigh,"hCovid_mKRT")
#Venn(human_ctrlHigh_mouseVer,mouse_KRTHigh,"hCtrl_mKRT")

hKrt_mKrt <- rbind(c(length(intersect(human_KRTHigh_mouseVer, mouse_KRTHigh)),
      		         length(intersect(human_AEHigh_mouseVer , mouse_KRTHigh))),
      		       c(length(intersect(human_KRTHigh_mouseVer, mouse_AEHigh )),
      		         length(intersect(human_AEHigh_mouseVer , mouse_AEHigh ))))

hKrtcovid_mKrt <- rbind(c(length(intersect(human_KRTHigh_covid_mouseVer, mouse_KRTHigh)),
      		              length(intersect(human_AEHigh_covid_mouseVer , mouse_KRTHigh))),
      		            c(length(intersect(human_KRTHigh_covid_mouseVer, mouse_AEHigh )),
      		              length(intersect(human_AEHigh_covid_mouseVer , mouse_AEHigh ))))
     
hCovid_mKrt <- rbind(c(length(intersect(human_covidHigh_mouseVer, mouse_KRTHigh)),
      		              length(intersect(human_ctrlHigh_mouseVer , mouse_KRTHigh))),
      		            c(length(intersect(human_covidHigh_mouseVer, mouse_AEHigh )),
      		              length(intersect(human_ctrlHigh_mouseVer , mouse_AEHigh ))))
     
fisher.test(hKrt_mKrt,alternative="greater")$p.val
fisher.test(hKrtcovid_mKrt,alternative="greater")$p.val
fisher.test(hCovid_mKrt,alternative="greater")$p.val


#########3 hKrt & mKrt & mAE
library(VennDiagram)
venn.plot <- draw.triple.venn(
area1 = length(human_KRTHigh_mouseVer),
area2 = length(mouse_KRTHigh),
area3 = length(mouse_AEHigh),
n12 = length(intersect(human_KRTHigh_mouseVer,mouse_KRTHigh)),
n23 = length(intersect(mouse_AEHigh,mouse_KRTHigh)),
n13 = length(intersect(human_KRTHigh_mouseVer,mouse_AEHigh)),
n123 = length(intersect(intersect(human_KRTHigh_mouseVer,mouse_KRTHigh),mouse_AEHigh ) ),
category = c("humanKrtHigh", "mouseKrtHigh", "mouseAeHigh"),
fill = c("red", "blue", "green"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.col = c("red", "blue", "green")
);
pdf(file="Venn_humanKrtHigh_mouseKrtHigh_mouseAeHigh.pdf")
grid.draw(venn.plot);
dev.off()
#
#human_KRTHigh_mouseVer
#mouse_KRTHigh
#mouse_AEHigh
#

Cd8cells <- names(cbdataHAR2_scale_allgenes$sample)[which(cbdataHAR2_scale_allgenes$sample %in% c("S1D","S2D"))]
IgGcells <- names(cbdataHAR2_scale_allgenes$sample)[which(cbdataHAR2_scale_allgenes$sample %in% c("S1A","S2A"))]


box2 <- function(g1,M){
Pv <- wilcox.test(apply(expTable[g1, IgGcells],1,mean),
				  apply(expTable[g1, Cd8cells],1,mean),
				  alternative="two.sided")$p.val
boxplot(apply(expTable[g1, IgGcells],1,mean),
				apply(expTable[g1, Cd8cells],1,mean),
				outline=F,ylab="aveExp",las=2,main=M,
				names=c("IgG","Cd8a"),col=c("red","blue"),xlab=Pv)
}

pdf(file="humanDEG_in_mouseCd8aIgG_box.pdf",width=10,height=4)
par(mar=c(10,4,4,2),mfrow=c(1,6))
box2(human_covidHigh_mouseVer,"covidHigh")
box2(human_ctrlHigh_mouseVer,"ctrlHigh")
box2(human_KRTHigh_mouseVer,"KrtHigh")
box2(human_AEHigh_mouseVer,"AEHigh")
box2(human_KRTHigh_covid_mouseVer,"KrtHigh covid")
box2(human_AEHigh_mouseVer,"AEHigh covid")
dev.off()


celltype3[Krt8on_spots] <- "Krt"
celltype3[AE_cells_new] <- "AE"

box2ct3 <- function(g1,M){
Pv <- wilcox.test(apply(expTable[g1, Krt8on_spots],1,mean),
				  apply(expTable[g1, AE_cells_new],1,mean),
				  alternative="two.sided")$p.val
boxplot(apply(expTable[g1, Krt8on_spots],1,mean),
				apply(expTable[g1, AE_cells_new],1,mean),
				outline=F,ylab="aveExp",las=2,main=M,
				names=c("Krt","AE"),col=c("red","blue"),xlab=Pv)
}

pdf(file="humanDEG_in_mouseKrtAe_box.pdf",width=10,height=4)
par(mar=c(10,4,4,2),mfrow=c(1,6))
box2ct3(human_covidHigh_mouseVer,"covidHigh")
box2ct3(human_ctrlHigh_mouseVer,"ctrlHigh")
box2ct3(human_KRTHigh_mouseVer,"KrtHigh")
box2ct3(human_AEHigh_mouseVer,"AEHigh")
box2ct3(human_KRTHigh_covid_mouseVer,"KrtHigh covid")
box2ct3(human_AEHigh_mouseVer,"AEHigh covid")
dev.off()


transitional_pathway <- read.table("pathway_ref_covidVSctrl/transitional_top20",row.names=1,sep="\t",comment.char = ";")
basal_pathway <- read.table("pathway_ref_covidVSctrl/basal_top20",row.names=1,sep="\t",comment.char = ";")
AT1_pathway<- read.table("pathway_ref_covidVSctrl/AT1_top20",row.names=1,sep="\t",comment.char = ";")
AT2_pathway<- read.table("pathway_ref_covidVSctrl/AT2_top20",row.names=1,sep="\t",comment.char = ";")


human2mouse <- function(data){
	all_genes <- paste(data, collapse = ",")
	# Split the string by commas to get individual gene names
	gene_names <- unlist(strsplit(all_genes, ","))
	# Remove any empty strings
	gene_names <- gene_names[gene_names != ""]
	gene_names <- gsub("[\n\"\\)]", "", gene_names)
	gene_names <- gsub(" ", "", gene_names)
	outdata <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% gene_names),"Symbol.x"],allgenes)
	return(outdata)
}

transitional_mm <- human2mouse(transitional_pathway[,4])
basal_mm <- human2mouse(basal_pathway[,4])
AT1_mm <- human2mouse(AT1_pathway[,4])
AT2_mm <- human2mouse(AT2_pathway[,4])



pdf(file="pathway_in_mouseCd8aIgG_box.pdf",width=7.5,height=4)
par(mar=c(10,4,4,2),mfrow=c(1,4))
box2(transitional_mm,"transitional")
box2(basal_mm,"basal")
box2(AT1_mm,"AT1")
box2(AT2_mm,"AT2")
dev.off()

pdf(file="pathway_in_mouseKrtAe_box.pdf",width=7.5,height=4)
par(mar=c(10,4,4,2),mfrow=c(1,4))
box2ct3(transitional_mm,"transitional")
box2ct3(basal_mm,"basal")
box2ct3(AT1_mm,"AT1")
box2ct3(AT2_mm,"AT2")
dev.off()



human_covid_g_ctrl_GSEA <- read.table("/scratch/sh8tv/lungST/humanST/humanLung_10x_spaceranger/human_covid_g_ctrl_GSEA_geneLists.txt",row.names=1,sep="\t")
allGSEAgenes <- c()
pdf(file="covid_g_ctrl_GSEApathway_mouseSTbox.pdf",width=24,height=10)
par(mfcol=c(2,12),mar=c(4,4,8,2))
for(i in 1:11){
	thisdata<- human_covid_g_ctrl_GSEA[i,]
	thisdata <- thisdata[which(thisdata != "")]
	mmdata <- intersect(mh_data[which(mh_data[,"Symbol.y"] %in% thisdata),"Symbol.x"],allgenes)
	allGSEAgenes <- c(allGSEAgenes, mmdata)
	box2(mmdata,gsub("_","\n",rownames(human_covid_g_ctrl_GSEA)[i]))
	box2ct3(mmdata,gsub("_","\n",rownames(human_covid_g_ctrl_GSEA)[i]))
}
box2(allGSEAgenes,"combine")
box2ct3(allGSEAgenes,"combine")
dev.off()


