library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)
#library(dplyr)
library(AUCell)
library(GSEABase)
library(ggplot2)
library(future)
library(tidyverse)
library(RCurl)
library(AnnotationHub)
library(ggplot2)
library(ArchR)


options(future.globals.maxSize = 8000 * 1024^2)
library(BiocParallel)
register(SerialParam())

plan("multiprocess", workers = 50)


load("glioblastoma_single_cell/1.RNA/1.RNA_seq_seurat/glio.int.all.peaks.chromvar.rdata")

p1 <- DimPlot(glio.int, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(glio.int, reduction = "umap", group.by = "CCA_RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(glio.int, reduction = "wnn.umap", group.by = "CCA_RNA_snn_res.0.2",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

col<-c("#B3CDE3","#FCCDE5","#B2DF8A","#E6AB02","#FDBF6F")

getwd()

DefaultAssay(glio.int) <- "CCA_RNA"
clusters_of_interest <- 0
subset_seurat <- subset(glio.int, idents = clusters_of_interest)
DimPlot(subset_seurat, group.by = "ident")

#subset_seurat2<-subset_seurat[,names(which(subset_seurat$sample != "core1"｜ subset_seurat$sample != "core2"｜ subset_seurat$sample != "core3" ))]

subset_seurat2<-subset_seurat[,!(names(subset_seurat$sample) %in% names(which(subset_seurat$sample == "core1"| subset_seurat$sample == "core2"| subset_seurat$sample == "core3" )))]

subset_seurat2

subset_seurat<-subset_seurat2

DefaultAssay(subset_seurat) <- "CCA_RNA"

subset_seurat<-FindNeighbors(subset_seurat,dims=1:30,reduction="pca",k.param=60,prune.SNN=1/20)
#filtered_seurat<-FindNeighbors(filtered_seurat,dims=1:30)
names(subset_seurat@graphs)



for( res in c(0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1)){
   subset_seurat<-FindClusters(subset_seurat,graph.name='CCA_RNA_snn',resolution=res,algorithm=2)
    #filtered_seurat<-FindClusters(filtered_seurat,resolution=res)
}

library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col<-sample(col_vector, n)

DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.1",cols = col) + ggtitle("louvain_0.1")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.15",cols = col) + ggtitle("louvain_0.1")

DimPlot(subset_seurat,reduction="rna.umap",group.by="CCA_RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.2")
DimPlot(subset_seurat,reduction="rna.umap",group.by="CCA_RNA_snn_res.0.3",cols = col) + ggtitle("louvain_0.3")
DimPlot(subset_seurat,reduction="rna.umap",group.by="CCA_RNA_snn_res.0.4",cols = col) + ggtitle("louvain_0.4")
#DimPlot(alldata.int,reduction="tsne",group.by="CCA_snn_res.0.4",cols = col) + ggtitle("louvain_0.4")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.6",cols = col) + ggtitle("louvain_0.6")
DimPlot(subset_seurat,reduction="rna.umap",group.by="CCA_RNA_snn_res.0.8",cols = col) + ggtitle("louvain_0.8")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.1",cols = col) + ggtitle("louvain_1")





DefaultAssay(subset_seurat) <- "peaks"


subset_seurat<-FindNeighbors(subset_seurat,dims=1:30,k.param=60,reduction = 'harmony',prune.SNN=1/20)
#filtered_seurat<-FindNeighbors(filtered_seurat,dims=1:30)
names(subset_seurat@graphs)

for( res in c(0.1,0.2,0.25,0.3,0.4,0.6,0.8,1)){
    subset_seurat<-FindClusters(subset_seurat,graph.name='peaks_snn',resolution=res,algorithm=1)
    #filtered_seurat<-FindClusters(filtered_seurat,resolution=res)
}

DimPlot(subset_seurat, reduction = "umap", group.by = "peaks_snn_res.0.1",cols = col) + ggtitle("louvain_0.1")
DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.2",cols = col) + ggtitle("louvain_0.2")
DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.25",cols = col) + ggtitle("louvain_0.25")

DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.3",cols = col) + ggtitle("louvain_0.3")
DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.4",cols = col) + ggtitle("louvain_0.4")
#DimPlot(alldata.int,reduction="tsne",group.by="CCA_snn_res.0.4",cols = col) + ggtitle("louvain_0.4")
DimPlot(subset_seurat, reduction = "umap", group.by = "peaks_snn_res.0.6",cols = col) + ggtitle("louvain_0.6")
DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.8",cols = col) + ggtitle("louvain_0.8")
DimPlot(subset_seurat, reduction = "umap", group.by = "peaks_snn_res.1",cols = col) + ggtitle("louvain_1")

DimPlot(subset_seurat, reduction = "umap", group.by = "peaks_snn_res.0.1",cols = c("#005084","#4690C2","#9DA7D1")) + ggtitle("louvain_0.1")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.1")

pdf("edge.sub.clsuters.pdf")
DimPlot(subset_seurat, reduction = "umap", group.by = "peaks_snn_res.0.1",cols = c("#005084","#4690C2","#9DA7D1")) + ggtitle("louvain_0.1")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.1")

dev.off()


library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col<-sample(col_vector, n)

DimPlot(subset_seurat, reduction = "rna.umap", group.by = "peaks_snn_res.0.1",cols = c("#005084","#4690C2","#9DA7D1")) + ggtitle("louvain_0.1")


DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.1")

DimPlot(subset_seurat, reduction = "umap", group.by = "peaks_snn_res.0.1",cols = c("#437A66","#6A926A","#65B097")) + ggtitle("louvain_0.1")

table(subset_seurat$CCA_RNA_snn_res.0.2,subset_seurat$peaks_snn_res.0.1)
col<-c("#005084","#4690C2","#9DA7D1")

DefaultAssay(subset_seurat) <- "CCA_RNA"
FeaturePlot(subset_seurat, features = 'NT5E',reduction="umap")
VlnPlot(subset_seurat, features = "NT5E",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
pdf("subcluster.edge.nt5e.pdf")
FeaturePlot(subset_seurat, features = 'NT5E',reduction="umap")
VlnPlot(subset_seurat, features = "NT5E",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
dev.off()

DefaultAssay(subset_seurat) <- "CCA_RNA"
FeaturePlot(subset_seurat, features = 'OSMR',reduction="umap")
VlnPlot(subset_seurat, features = "OSMR",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
pdf("subcluster.edge.osmr.pdf")
FeaturePlot(subset_seurat, features = 'OSMR',reduction="umap")
VlnPlot(subset_seurat, features = "OSMR",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
dev.off()

FeaturePlot(subset_seurat, features = 'LIFR',reduction="umap")
FeaturePlot(subset_seurat, features = 'LIFR',reduction="umap")
VlnPlot(subset_seurat, features = "LIFR",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
pdf("subcluster.edge.LIFR.pdf")
FeaturePlot(subset_seurat, features = 'LIFR',reduction="umap")
VlnPlot(subset_seurat, features = "LIFR",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
dev.off()

a<-table(subset_seurat$CCA_RNA_snn_res.0.2,subset_seurat$peaks_snn_res.0.1)
new<-t(a)
new2<-data.frame(new)
new2$Freq <- new2$Freq
new
new2

new2$Var1<-factor(new2$Var1)
new2$Var2<-factor(new2$Var2)
new2$Var1
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()
pdf("edge.overlap.barplot.scRNA_scATAC.pdf")
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()
dev.off()

a<-table(subset_seurat$sample,subset_seurat$peaks_snn_res.0.1)
new<-t(a)
new2<-data.frame(new)
new2$Freq <- new2$Freq
new
new2
new2$Var1<-factor(new2$Var1)
new2$Var2<-factor(new2$Var2)
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()
pdf("edge.egde.cluser.absolute.pdf")
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()
dev.off()

a<-table(subset_seurat$sample,subset_seurat$peaks_snn_res.0.1)
a
new<-t(a/rowSums(a))
new2<-data.frame(new)
new2$Freq <- new2$Freq

new
new2
new2$Var1<-factor(new2$Var1)
new2$Var2<-factor(new2$Var2)
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()
pdf("edge.egde.cluser.percentage.pdf")
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()
dev.off()

a<-table(subset_seurat$sample,subset_seurat$CCA_RNA_snn_res.0.2)
new<-t(a)
new2<-data.frame(new)
new2$Freq <- new2$Freq
new
new2
new2$Var1<-factor(new2$Var1)
new2$Var2<-factor(new2$Var2)
ggplot(new2, aes( y=Freq, x=Var2,fill=Var1))+geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_classic()

load('/disk1/pengweixing/lene_sc_glio/02.single_cell/03.Edge_analysis/02.quad_meta_model/Neftel.Rdata')
core_matrix = as.matrix(subset_seurat@assays$RNA@data)
core_rankings <- AUCell_buildRankings(core_matrix, nCores=30, plotStats=FALSE)
rm(core_matrix)

MES2_core_AUC <- AUCell_calcAUC(MES2geneSets, core_rankings)
MES1_core_AUC <- AUCell_calcAUC(MES1geneSets, core_rankings)
NPC1_core_AUC <- AUCell_calcAUC(NPC1geneSets, core_rankings)
NPC2_core_AUC <- AUCell_calcAUC(NPC2geneSets, core_rankings)
AC_core_AUC <- AUCell_calcAUC(ACgeneSets, core_rankings)
OPC_core_AUC <- AUCell_calcAUC(OPCgeneSets, core_rankings)

core_AUC_list = c('MES1_core_AUC','MES2_core_AUC','NPC1_core_AUC',
                              'NPC2_core_AUC','AC_core_AUC','OPC_core_AUC')
for(each in core_AUC_list){
    AUC <- t(as.data.frame(get(each)@assays@data$AUC))
    colnames(AUC) <- paste0(colnames(AUC), "_AUC")
    AUC <- data.frame(AUC)
    #print(AUC)
    subset_seurat <- AddMetaData(subset_seurat, metadata = AUC)
}

core_meta <- subset_seurat@meta.data

core_meta$MES_mean <- apply(cbind(core_meta$MES1_AUC,core_meta$MES2_AUC),1,mean)
core_meta$NPC_mean <- apply(cbind(core_meta$NPC1_AUC,core_meta$NPC2_AUC),1,mean)
a1 <- apply(cbind(core_meta$OPC_AUC, core_meta$NPC_mean), 1, max)
b1 <- apply(cbind(core_meta$AC_AUC, core_meta$MES_mean), 1, max)
core_meta$Y.axis <-  a1 - b1
core_meta$CellClass <- ifelse(core_meta$Y.axis > 0, "OPC.NPC", "AC.MES")

head(core_meta)

cc <- cbind(core_meta$AC_AUC, core_meta$MES_mean, core_meta$NPC_mean, core_meta$OPC_AUC)
colnames(cc) <- c("AC", "MES", "NPC", "OPC")
core_meta$MaxClass <- apply(cc,1,function(x) which(x==max(x)))
core_meta$MaxClass <- gsub(1, "AC", core_meta$MaxClass)
core_meta$MaxClass <- gsub(2, "MES", core_meta$MaxClass)
core_meta$MaxClass <- gsub(3, "NPC", core_meta$MaxClass)
core_meta$MaxClass <- gsub(4, "OPC", core_meta$MaxClass)
    
core_meta$ClassSign_X <- core_meta$MaxClass
core_meta$ClassSign_X <- gsub("AC", -1, core_meta$ClassSign_X )
core_meta$ClassSign_X <- gsub("MES", 1,core_meta$ClassSign_X )
core_meta$ClassSign_X <- gsub("OPC", -1, core_meta$ClassSign_X )
core_meta$ClassSign_X <- gsub("NPC", 1, core_meta$ClassSign_X )

core_meta$OPC.NPC_x <- log(abs(core_meta$OPC_AUC - core_meta$NPC_mean)+1,1.6)
core_meta$AC.MES_x <- log(abs(core_meta$AC_AUC - core_meta$MES_mean)+1,2)
core_meta$X.axis <- NA
core_meta$X.axis[grep("OPC.NPC", core_meta$CellClass)] <- core_meta$OPC.NPC_x[grep("OPC.NPC", core_meta$CellClass)]
core_meta$X.axis[grep("AC.MES", core_meta$CellClass)] <- core_meta$AC.MES_x[grep("AC.MES", core_meta$CellClass)]
#now change the sign based on maxclass for cell
#if AC is the max value, make postive; if MES is the max, make negative
#OPC = +, NPC = -
core_meta$X.axis_Class <- core_meta$X.axis * as.numeric(core_meta$ClassSign_X)

subset_seurat@meta.data

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

library(RColorBrewer)
core_meta$density <- get_density(core_meta$X.axis_Class, core_meta$Y.axis, n = 100)

core_meta<-as.data.frame(core_meta)
#core_meta$X.axis_Class <- factor(core_meta$X.axis_Class)

d <- ggplot(core_meta, aes(x = X.axis_Class, y = Y.axis, color = density)) +
geom_point(aes(X.axis_Class, Y.axis, color = density), alpha = 1, size = 0.7) + geom_density_2d(bins = 3, color = "grey")+scale_colour_gradientn(colours = rev(colorRampPalette(brewer.pal(11,"Greys"))(31))) +
xlab("<-- AC-like      MES-like -->") +
ylab("<-- AC-like  OPC-like -->") + 
geom_vline(xintercept = 0, lty = 2, col = "black") + 
geom_hline(yintercept = 0, lty = 2, col = "black") + 
ggtitle("<-- OPC-like      NPC-like -->")  +
 theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+ylim(-0.3,0.3)+xlim(-0.5,0.5)
d

core_meta2<-data.frame(meta=core_meta$MaxClass)
rownames(core_meta2)<-rownames(core_meta)
subset_seurat <- AddMetaData(subset_seurat, metadata = core_meta2)

col<-c("#CA0020","#F4A582","#92C5DE","#0571B0")
col<-c("#F4A582","#CA0020","#0571B0","#92C5DE")
col<-c("#D64E4C","#BD2329","#33708F","#428DB4")

DimPlot(subset_seurat, reduction = "umap", group.by = "meta",cols = col ) + ggtitle("RNA")


table(core_meta$MaxClass,core_meta$peaks_snn_res.0.1)

common<-table(core_meta$MaxClass,core_meta$peaks_snn_res.0.1)
new<-t(t(common)/colSums(common))
new2<-data.frame(new)
new2$Freq <- new2$Freq*100
new2$Var1<-factor(new2$Var1,levels=c("NPC","OPC","AC","MES"))
col<-c("#CA0020","#F4A582","#92C5DE","#0571B0")
ggplot(new2, aes( y=Freq, x=Var2))+scale_fill_manual(values=c("#0571B0","#92C5DE","#F4A582","#CA0020"))+geom_col(aes(fill = Var1 ), width = 0.7)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_y_continuous(breaks = round(seq(0, 100, by = 20),1))

pdf("subclusters.meta.modules.edge.pdf")
col<-c("#CA0020","#F4A582","#92C5DE","#0571B0")
ggplot(new2, aes( y=Freq, x=Var2))+scale_fill_manual(values=c("#0571B0","#92C5DE","#F4A582","#CA0020"))+geom_col(aes(fill = Var1 ), width = 0.7)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_y_continuous(breaks = round(seq(0, 100, by = 20),1))
dev.off()

save(subset_seurat,file="edge_cells.rdata")

DefaultAssay(subset_seurat) <- "CCA_RNA"
FeaturePlot(subset_seurat, features = 'NT5E',reduction="umap")
VlnPlot(subset_seurat, features = "NT5E",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
pdf("subcluster.edge.nt5e.pdf")
VlnPlot(subset_seurat, features = "NT5E",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
dev.off()

DefaultAssay(subset_seurat) <- "CCA_RNA"
FeaturePlot(subset_seurat, features = 'OSMR',reduction="umap")
VlnPlot(subset_seurat, features = "OSMR",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
pdf("subcluster.edge.osmr.pdf")
VlnPlot(subset_seurat, features = "OSMR",group.by = "peaks_snn_res.0.1",pt.size=0,cols=c("#005084","#4690C2","#9DA7D1"))
dev.off()

DefaultAssay(subset_seurat) <- "CCA_RNA"
FeaturePlot(subset_seurat, features = 'HLA-DPA1',reduction="umap")
VlnPlot(subset_seurat, features = "HLA-DPA1",group.by = "peaks_snn_res.0.1",pt.size=0)

DefaultAssay(subset_seurat) <- "CCA_RNA"
markers_rna <- presto:::wilcoxauc.Seurat(X = subset_seurat, group_by = 'peaks_snn_res.0.1', assay = 'data', seurat_assay = 'CCA_RNA')
DefaultAssay(subset_seurat) <- "peaks"
markers_motifs <- presto:::wilcoxauc.Seurat(X = subset_seurat, group_by = 'peaks_snn_res.0.1', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(subset_seurat, id = motif.names)

DefaultAssay(subset_seurat) <- "peaks"

markers_chromatin <- presto:::wilcoxauc.Seurat(X = subset_seurat, group_by = 'peaks_snn_res.0.1', assay = 'data', seurat_assay = 'peaks')


topTFs_chromatin <- function(peaks_snn_res.0.1, padj.cutoff = 0.05) {
  ctmarkers_chromatin <- dplyr::filter(
    markers_chromatin, group == peaks_snn_res.0.1, padj < 0.05, logFC > 0.3,auc > 0.5) 
    top_tfs <- ctmarkers_chromatin
  return(top_tfs)
}

tf1_chromatin<-topTFs_chromatin("0")
tf2_chromatin<-topTFs_chromatin("1")
tf3_chromatin<-topTFs_chromatin("2")

tf_chromatin<-rbind(tf1_chromatin,tf2_chromatin,tf3_chromatin)


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ChIPseeker)
library(org.Hs.eg.db)
library(annotate)

split1<-strsplit(tf1_chromatin$feature,split = "-")
split1<-t(data.frame(split1))
split2<-data.frame(split1[,1],as.numeric(split1[,2]),as.numeric(split1[,3]))
input<-split2
cluster2_anno<-data.frame(seqnames=input[,1],start=input[,2],end=input[,3],width=input[,3]-input[,2],strand=rep("*",length(input[,2])))
rownames(cluster2_anno)<-paste(input[,1],input[,2],input[,3],sep="_")
cluster_c<-cluster2_anno
cluster_c_g<-makeGRangesFromDataFrame(cluster_c, keep.extra.columns = TRUE)
anno_LiverMinusHindbrain <- annotatePeak(cluster_c_g,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_LiverMinusHindbrain_GRanges <- as.GRanges(anno_LiverMinusHindbrain)
a<-getSYMBOL(anno_LiverMinusHindbrain_GRanges$geneId, data='org.Hs.eg')
anno_LiverMinusHindbrain_GRanges$genes<-a
b<-anno_LiverMinusHindbrain_GRanges[which(anno_LiverMinusHindbrain_GRanges$distanceToTSS < 100000 & anno_LiverMinusHindbrain_GRanges$distanceToTSS > -100000),]
b$genes
kk2<-enrichKEGG(b$geneId,organism='hsa', pvalueCutoff = 0.05)
kk3<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kk4 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
kk4
dotplot(kk4, showCategory=20)

split1<-strsplit(tf2_chromatin$feature,split = "-")
split1<-t(data.frame(split1))
split2<-data.frame(split1[,1],as.numeric(split1[,2]),as.numeric(split1[,3]))
input<-split2
cluster2_anno<-data.frame(seqnames=input[,1],start=input[,2],end=input[,3],width=input[,3]-input[,2],strand=rep("*",length(input[,2])))
rownames(cluster2_anno)<-paste(input[,1],input[,2],input[,3],sep="_")
cluster_c<-cluster2_anno
cluster_c_g<-makeGRangesFromDataFrame(cluster_c, keep.extra.columns = TRUE)
anno_LiverMinusHindbrain <- annotatePeak(cluster_c_g,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_LiverMinusHindbrain_GRanges <- as.GRanges(anno_LiverMinusHindbrain)
a<-getSYMBOL(anno_LiverMinusHindbrain_GRanges$geneId, data='org.Hs.eg')
anno_LiverMinusHindbrain_GRanges$genes<-a
b<-anno_LiverMinusHindbrain_GRanges[which(anno_LiverMinusHindbrain_GRanges$distanceToTSS < 100000 & anno_LiverMinusHindbrain_GRanges$distanceToTSS > -100000),]
kk2<-enrichKEGG(b$geneId,organism='hsa', pvalueCutoff = 0.05)
kk3<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kk4 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
dotplot(kk4, showCategory=20)
pdf("edge.cluster2.pdf")
dotplot(kk4, showCategory=20)
dev.off()

split1<-strsplit(tf3_chromatin$feature,split = "-")
split1<-t(data.frame(split1))
split2<-data.frame(split1[,1],as.numeric(split1[,2]),as.numeric(split1[,3]))
input<-split2
cluster2_anno<-data.frame(seqnames=input[,1],start=input[,2],end=input[,3],width=input[,3]-input[,2],strand=rep("*",length(input[,2])))
rownames(cluster2_anno)<-paste(input[,1],input[,2],input[,3],sep="_")
cluster_c<-cluster2_anno
cluster_c_g<-makeGRangesFromDataFrame(cluster_c, keep.extra.columns = TRUE)
anno_LiverMinusHindbrain <- annotatePeak(cluster_c_g,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_LiverMinusHindbrain_GRanges <- as.GRanges(anno_LiverMinusHindbrain)
a<-getSYMBOL(anno_LiverMinusHindbrain_GRanges$geneId, data='org.Hs.eg')
anno_LiverMinusHindbrain_GRanges$genes<-a
b<-anno_LiverMinusHindbrain_GRanges[which(anno_LiverMinusHindbrain_GRanges$distanceToTSS < 100000 & anno_LiverMinusHindbrain_GRanges$distanceToTSS > -100000),]
b$genes
kk2<-enrichKEGG(b$geneId,organism='hsa', pvalueCutoff = 0.05)
kk3<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kk4 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
dotplot(kk4, showCategory=30)

b$geneId
write.table(b$geneId,file="chromatin.edge.3.xls",quote=F,row.names=F)

getwd()

topTFs <- function(peaks_snn_res.0.1, padj.cutoff = 0.01) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == peaks_snn_res.0.1, RNA.padj < padj.cutoff, RNA.logFC > 0.1,RNA.auc>0.1 )
#    %>% 
#    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == peaks_snn_res.0.1, motif.padj < padj.cutoff, motif.logFC > 0.1,motif.auc > 0.1) 
#    %>% 
#    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

topTFs_rna <- function(peaks_snn_res.0.1, padj.cutoff = 0.01) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == peaks_snn_res.0.1, RNA.padj < 0.01, RNA.logFC > 0.5,RNA.auc>0.6 ) 
  top_tfs <- ctmarkers_rna
  return(top_tfs)
}

topTFs_motif <- function(peaks_snn_res.0.1, padj.cutoff = 0.01) {
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == peaks_snn_res.0.1, motif.padj < padj.cutoff, motif.logFC > 0.5,motif.auc > 0.6) 
  top_tfs <- ctmarkers_motif
  return(top_tfs)
}

tf1_c2<-topTFs("0")
tf2_c2<-topTFs("1")
tf3_c2<-topTFs("2")


tf_c2<-rbind(tf1_c2,tf2_c2,tf3_c2)

tf_c2

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "FOXP2",group.by = "peaks_snn_res.0.1",cols=col,pt.size=0)

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "OSMR",group.by = "peaks_snn_res.0.1",cols=col,pt.size=0)

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "LIFR",group.by = "peaks_snn_res.0.1",cols=col,pt.size=0)

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "IL6ST",group.by = "peaks_snn_res.0.1",cols=col,pt.size=0)

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "NF1",group.by = "peaks_snn_res.0.1",cols=col,pt.size=0)

DefaultAssay(glio.int) <- "CCA_RNA"
#DotPlot(object = subset_seurat, features = unique(tf_c2$gene),group.by = "peaks_snn_res.0.1")+coord_flip() 

tf1_rna<-topTFs_rna("0")
tf2_rna<-topTFs_rna("1")
tf3_rna<-topTFs_rna("2")
tf4_rna<-topTFs_rna("3")
tf5_rna<-topTFs_rna("4")

tf_rna<-rbind(tf1_rna,tf2_rna,tf3_rna,tf4_rna,tf5_rna)

tf1_motif<-topTFs_motif("0")
tf2_motif<-topTFs_motif("1")
tf3_motif<-topTFs_motif("2")

tf_motif<-rbind(tf1_motif,tf2_motif,tf3_motif)

write.table(tf_motif,file="tf.motif.subcluster.without.core.xls ",quote=F,sep="\t",row.names=F)

#tf_chromatin
write.table(tf_chromatin,file="chromatin.subcluster.xls",quote=F,sep="\t",row.names=F)

DefaultAssay(subset_seurat) <- "CCA_RNA"
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
#filtered_seurat<-ScaleData(filtered_seurat,features=as.character(unique(top5$gene)),assay="RNA")
DoHeatmap(subset_seurat, features = as.character(tf_rna$RNA.feature), group.by = "peaks_snn_res.0.1")+scale_fill_gradientn(colours = rev(mapal))+theme(text = element_text(size = 5))

pdf("edge.cluster.rna.heatmap.pdf")

#filtered_seurat<-ScaleData(filtered_seurat,features=as.character(unique(top5$gene)),assay="RNA")
DoHeatmap(subset_seurat, features = as.character(tf_rna$RNA.feature), group.by = "CCA_RNA_snn_res.0.2")+scale_fill_gradientn(colours = rev(mapal))+theme(text = element_text(size = 5))
dev.off()

library(ComplexHeatmap)
library(pheatmap)
library( "RColorBrewer" )
annotation<-data.frame(glio.int.kmeans_3=subset_seurat$`peaks_snn_res.0.1`)
rownames(annotation)<-names(subset_seurat$`peaks_snn_res.0.1`)
#annotation<-NULL
annotation$glio.int.kmeans_3<-factor(subset_seurat$`peaks_snn_res.0.1`)
#color = colorRampPalette(brewer.pal(3,"RdBu"))(21)
color <- colorRampPalette(RColorBrewer::brewer.pal(7,"RdBu"))(21)
breaks<-seq(-1.5,1.5,by=0.15)

annotation3<-data.frame(group=annotation[order(annotation$glio.int.kmeans_3),])
rownames(annotation3)<-rownames(annotation)[order(annotation$glio.int.kmeans_3)]
#annotation3$group<-factor(annotation3$group)


#chromatin 
DefaultAssay(subset_seurat) <- "peaks"

s<-as.matrix(subset_seurat@assays$peaks@data[tf_chromatin$feature,])
s2<-as.matrix(s[,rownames(annotation3)])



pheatmap(s2, show_rownames=T, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="none",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)
fig<-pheatmap(s2, show_rownames=T, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="none",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)
pdf("chromatin.edge.heatmap.pdf")
fig
dev.off()

m<-pheatmap(s2, show_rownames=T, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="none",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)
pdf("heatmap.chroamtin.accessibility.edge.pdf")
m
dev.off()

DefaultAssay(glio.int) <- "CCA_RNA"
markers_rna <- presto:::wilcoxauc.Seurat(X =subset_seurat , group_by = 'peaks_snn_res.0.1', assay = 'data', seurat_assay = 'CCA_RNA')
DefaultAssay(glio.int) <- "peaks"
markers_motifs <- presto:::wilcoxauc.Seurat(X = subset_seurat, group_by = 'peaks_snn_res.0.1', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(glio.int, id = motif.names)

topTFs_motif <- function(peaks_snn_res.0.1, padj.cutoff = 0.01) {
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == peaks_snn_res.0.1, motif.padj < padj.cutoff, motif.logFC > 0.5,motif.auc > 0.6) 
  top_tfs <- ctmarkers_motif
  return(top_tfs)
}

topTFs_rna <- function(peaks_snn_res.0.1, padj.cutoff = 0.01) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == peaks_snn_res.0.1, RNA.padj < 0.01, RNA.logFC > 0.5,RNA.auc>0.5 ) 
  top_tfs <- ctmarkers_rna
  return(top_tfs)
}

tf1_motif<-topTFs_motif("0")
tf2_motif<-topTFs_motif("1")
tf3_motif<-topTFs_motif("2")

tf_motif<-rbind(tf1_motif,tf2_motif,tf3_motif)

tf1_rna<-topTFs_rna("0")
tf2_rna<-topTFs_rna("1")
tf3_rna<-topTFs_rna("2")
tf4_rna<-topTFs_rna("3")
tf5_rna<-topTFs_rna("4")

tf_rna<-rbind(tf1_rna,tf2_rna,tf3_rna,tf4_rna,tf5_rna)

tf_rna

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "CDH6",group.by = "peaks_snn_res.0.1",cols=col,pt.size=0)

cols<-c("#005084","#4690C2","#9DA7D1")
VlnPlot(subset_seurat, features ="GRIA1",group.by = "peaks_snn_res.0.1",cols=cols,pt.size=0)
VlnPlot(subset_seurat, features ="GRIA2",group.by = "peaks_snn_res.0.1",cols=cols,pt.size=0)
VlnPlot(subset_seurat, features ="GRIA3",group.by = "peaks_snn_res.0.1",cols=cols,pt.size=0)
VlnPlot(subset_seurat, features ="GRIA4",group.by = "peaks_snn_res.0.1",cols=cols,pt.size=0)

VlnPlot(subset_seurat, features ="BDNF",group.by = "peaks_snn_res.0.1",cols=cols,pt.size=0)
VlnPlot(subset_seurat, features ="DLG4",group.by = "peaks_snn_res.0.1",pt.size=0)
VlnPlot(subset_seurat, features ="NLGN3",group.by = "peaks_snn_res.0.1",pt.size=0)
VlnPlot(subset_seurat, features ="NTRK2",group.by = "peaks_snn_res.0.1",cols=cols,pt.size=0)

#write.table(tf_motif,file="tf.motif.subcluster.xls",quote=F,sep="\t",row.names=F)

annotation3<-data.frame(group=annotation[order(annotation$glio.int.kmeans_3),])
rownames(annotation3)<-rownames(annotation)[order(annotation$glio.int.kmeans_3)]
#annotation3$group<-factor(annotation3$group)


s<-as.matrix(glio.int@assays$chromvar@data[tf_motif$motif.feature,])
s2<-as.matrix(s[,rownames(annotation3)])



pheatmap(s2, show_rownames=T, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)


tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "0"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
library(dplyr)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc	 >0.6 ,"sig", "non-significant"))

#tf_motif2<-tf_motif
#tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
#tf_motif2<-tf_motif2[which(tf_motif2$motif.group == "0"),]
ggplot(tf_motif2, aes(x=motif.logFC, y=motif.padj)) + geom_point()+theme_classic()

max(tf_motif2$motif.padj)

tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "0"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-306.363256720647+10
library(dplyr)
library(ggrepel)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >2,"sig", "non-significant"))
which(mydata$motif.feature == head(mydata[which(mydata$threshold == "sig"),]$motif.feature,10))

ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+ geom_hline(yintercept = 2,linetype="dashed", color = "black", size=0.5)+geom_vline(xintercept = 0.5,linetype="dashed", color = "black", size=0.5)+ylim(c(0,400))+theme_classic()
a<-ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+ geom_hline(yintercept = 2,linetype="dashed", color = "black", size=0.5)+geom_vline(xintercept = 0.5,linetype="dashed", color = "black", size=0.5)+ylim(c(0,400))+theme_classic()
pdf("volcano.edge.c0.pdf")
a
dev.off()

tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "0"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-306.363256720647+10
library(dplyr)
library(ggrepel)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >2,"sig", "non-significant"))
which(mydata$motif.feature == head(mydata[which(mydata$threshold == "sig"),]$motif.feature,10))

ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+geom_text_repel(data=subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >305),aes(label = gene),max.overlaps = 10)+ylim(c(0,400))+theme_classic()
a<-ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+geom_text_repel(data=subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >305),aes(label = gene),max.overlaps = 10)+ylim(c(0,400))+theme_classic()
pdf("volcano.edge.c0-add.pdf")
a
dev.off()

g<-subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >305)
g2<-g[order(g$motif.padj,decreasing =T),]
g2[order(g2$motif.logFC,decreasing =T),]





tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "1"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
#tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-356.363256720647+10
tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-306.363256720647+10
library(dplyr)
library(ggrepel)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >2,"sig", "non-significant"))
which(mydata$motif.feature == head(mydata[which(mydata$threshold == "sig"),]$motif.feature,10))

ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+ geom_hline(yintercept = 2,linetype="dashed", color = "black", size=0.5)+geom_vline(xintercept = 0.5,linetype="dashed", color = "black", size=0.5)+ylim(c(0,400))+theme_classic()
a<-ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+ geom_hline(yintercept = 2,linetype="dashed", color = "black", size=0.5)+geom_vline(xintercept = 0.5,linetype="dashed", color = "black", size=0.5)+ylim(c(0,400))+theme_classic()
pdf("volcano.edge.c1.pdf")
a
dev.off()

tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "1"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
#tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-356.363256720647+10
tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-306.363256720647+10
library(dplyr)
library(ggrepel)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >2,"sig", "non-significant"))
which(mydata$motif.feature == head(mydata[which(mydata$threshold == "sig"),]$motif.feature,10))

ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=1) +scale_colour_manual(values = c("black","#f2938e"))+geom_text_repel(data=subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >315),aes(label = gene),max.overlaps = 30)+ylim(c(0,400))+theme_classic()
a<-ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=1) +scale_colour_manual(values = c("black","#f2938e"))+geom_text_repel(data=subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >315),aes(label = gene),max.overlaps = 30)+ylim(c(0,400))+theme_classic()
pdf("volcano.edge.c1-add.pdf")
a
dev.off()

g<-subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >315)
g2<-g[order(g$motif.padj,decreasing =T),]
g2[order(g2$motif.logFC,decreasing =T),]

tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "2"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
#tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-356.363256720647+10
library(dplyr)
library(ggrepel)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >2,"sig", "non-significant"))
which(mydata$motif.feature == head(mydata[which(mydata$threshold == "sig"),]$motif.feature,10))

ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+ geom_hline(yintercept = 2,linetype="dashed", color = "black", size=0.5)+geom_vline(xintercept = 0.5,linetype="dashed", color = "black", size=0.5)+ylim(c(0,400))+theme_classic()
a<-ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+ geom_hline(yintercept = 2,linetype="dashed", color = "black", size=0.5)+geom_vline(xintercept = 0.5,linetype="dashed", color = "black", size=0.5)+ylim(c(0,400))+theme_classic()
pdf("volcano.edge.c2.pdf")
a
dev.off()

tf_motif2<-markers_motifs[which(markers_motifs$motif.group == "2"),]
tf_motif2$motif.padj<--log10(tf_motif2$motif.padj)
#tf_motif2$motif.padj[is.infinite(tf_motif2$motif.padj)]<-356.363256720647+10
library(dplyr)
library(ggrepel)
mydata<-tf_motif2%>%mutate(threshold = ifelse(motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >2,"sig", "non-significant"))
which(mydata$motif.feature == head(mydata[which(mydata$threshold == "sig"),]$motif.feature,10))

ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+geom_text_repel(data=subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >122),aes(label = gene),max.overlaps = 30)+ylim(c(0,200))+theme_classic()
a<-ggplot(mydata, aes(x=motif.logFC, y=motif.padj))+geom_point(aes(colour = threshold), size=3) +scale_colour_manual(values = c("black","#f2938e"))+geom_text_repel(data=subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >122),aes(label = gene),max.overlaps = 30)+ylim(c(0,200))+theme_classic()
pdf("volcano.edge.c2-add.pdf")
a
dev.off()

g<-subset(mydata,motif.logFC>= 0.5 & motif.auc>0.6 & motif.padj >122)
g2<-g[order(g$motif.padj,decreasing =T),]
g2[order(g2$motif.logFC,decreasing =T),]

markers_motifs_new<-markers_motifs
markers_motifs_new$motif.pval<--log10(markers_motifs_new$motif.pval)

tf_motif22<-markers_motifs_new[which(markers_motifs_new$motif.group == "2"),]
tf_motif21<-markers_motifs_new[which(markers_motifs_new$motif.group == "1"),]
tf_motif20<-markers_motifs_new[which(markers_motifs_new$motif.group == "0"),]

tf_all<-cbind(tf_motif20[,c(7)],tf_motif21[,c(7)],tf_motif22[,c(7,11)])
colnames(tf_all)<-c("0","1","2")
tf_all2<-tf_all[,c(1,2,3)]
rownames(tf_all2)<-tf_all[,4]
#tf_all2<-tf_all2[is.finite(rowSums(tf_all2)),]
tf_all2[sapply(tf_all2, is.infinite)] <- 306.363256720647+10

library(pheatmap)
cols = colorRampPalette(c("white", "red"))(30)
pheatmap(tf_all2,color = cols,fontsize_row=1,cluster_cols = F)
a<-pheatmap(tf_all2,color = cols,fontsize_row=1,cluster_cols =F)
pdf("pvalue.edge.pdf")
a
dev.off()

log10(6.3059239497785e-148)
log10(4.57988696098355e-146 )

DefaultAssay(subset_seurat) <- "CCA_RNA"
markers_rna <- presto:::wilcoxauc.Seurat(X = subset_seurat, group_by = 'CCA_RNA_snn_res.0.2', assay = 'data', seurat_assay = 'CCA_RNA')
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
markers_rna$gene <- markers_rna$RNA.feature

topTFs_rna <- function(CCA_RNA_snn_res.0.2, padj.cutoff = 0.01) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == CCA_RNA_snn_res.0.2, RNA.padj < 0.01, RNA.logFC > 0.2,RNA.auc>0.6 ) 
  top_tfs <- ctmarkers_rna
  return(top_tfs)
}

tf1_rna<-topTFs_rna("0")
tf2_rna<-topTFs_rna("1")
tf3_rna<-topTFs_rna("2")
tf4_rna<-topTFs_rna("3")
tf5_rna<-topTFs_rna("4")

tf_rna<-rbind(tf1_rna,tf2_rna,tf3_rna)


s<-as.matrix(glio.int@assays$CCA_RNA@scale.data[as.character(tf_rna$RNA.feature),])
s2<-as.matrix(s[,rownames(annotation3)])

pheatmap(s2, show_rownames=T, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=6,fontsize_col=4,breaks=breaks)


DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "NFIA",group.by = "CCA_RNA_snn_res.0.2",cols=col,pt.size=0)

write.table(tf_rna,file="subcluster.rna.xls",quote=F,sep="\t",row.names=F)

FeaturePlot(subset_seurat, features = 'GRIA4',reduction="rna.umap")
FeaturePlot(subset_seurat, features = 'COL6A1',reduction="rna.umap")


DefaultAssay(subset_seurat) <- "CCA_RNA"
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
#filtered_seurat<-ScaleData(filtered_seurat,features=as.character(unique(top5$gene)),assay="RNA")
DoHeatmap(subset_seurat, features = as.character(tf_rna$RNA.feature), group.by = "CCA_RNA_snn_res.0.2")+scale_fill_gradientn(colours = rev(mapal))+theme(text = element_text(size = 5))

DefaultAssay(subset_seurat) <- "CCA_RNA"
DotPlot(object = subset_seurat, features = unique(tf_c2$gene),group.by = "peaks_snn_res.0.1")+coord_flip() 


pdf("edge.doptplot.pdf")
DotPlot(object = subset_seurat, features = unique(tf_c2$gene),group.by = "peaks_snn_res.0.1")+coord_flip() 
dev.off()


#FeaturePlot(subset_seurat, features = 'NFIA',reduction="rna.umap")
FeaturePlot(subset_seurat, features = 'FOXP2',reduction="rna.umap")
#FeaturePlot(subset_seurat, features = 'ETV1',reduction="rna.umap")
#FeaturePlot(subset_seurat, features = 'PRRX1',reduction="rna.umap")


getwd()


DefaultAssay(subset_seurat) <- "peaks"
subset_seurat <- RegionStats(subset_seurat, genome = BSgenome.Hsapiens.UCSC.hg38)


peaklink <- LinkPeaks(
  object = subset_seurat,
  peak.assay = "peaks",
  expression.assay = "CCA_RNA",
  genes.use = c("RUNX2","FOXP2")
)

p1 <- CoveragePlot(
  object = peaklink,
  region = "FOXP2",
  features = "FOXP2",
  expression.assay = "CCA_RNA",
  group.by = "peaks_snn_res.0.1",
#  idents = c("0","1","2"),
  extend.upstream = 500,
  extend.downstream = 10000
)
patchwork::wrap_plots(p1, ncol = 1)

p2 <- CoveragePlot(
  object = peaklink,
  region = "RUNX2",
  features = "RUNX2",
  expression.assay = "CCA_RNA",
  group.by = "peaks_snn_res.0.1",
  extend.upstream = 500,
  extend.downstream = 10000
)
patchwork::wrap_plots(p2, ncol = 1)

pdf("panel.edge.runx2.pdf")
patchwork::wrap_plots(p2, ncol = 1)
dev.off()

pdf("link.scatac.without.filter.pdf")
patchwork::wrap_plots(p1, ncol = 1)
patchwork::wrap_plots(p2, ncol = 1)
dev.off()

DefaultAssay(subset_seurat) <- "CCA_RNA"

variable<-FindVariableFeatures(subset_seurat,method="vst",nfeatures=8000,verbose=T)
gene_name<-VariableFeatures(variable)

peaklink2 <- LinkPeaks(
  object = subset_seurat,
    distance = 500000,
  min.cells = 1000,
  peak.assay = "peaks",
  expression.assay = "CCA_RNA",
 genes.use = gene_name
)

table<-data.frame(peaklink2@assays$peaks@links)


table2<-table[which(table$pvalue < 0.01),]

unique(table2$gene)

num<-as.numeric(table(as.character(table2$gene)))

table(num)

names(table(num))

data.frame(name=c(names(table(num))[1:7],">8"),value=c(as.numeric(table(num))[1:7],18))
table3<-data.frame(name=as.character(c(names(table(num))[1:7],">8")),value=c(as.numeric(table(num))[1:7],18))

table3

library(ggplot2)
ggplot(table3, aes( y=value, x=name))+geom_bar(stat="identity", fill="steelblue")+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("barplot.all.pdf")
ggplot(table3, aes( y=value, x=name))+geom_bar(stat="identity", fill="steelblue")+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



s<-as.matrix(subset_seurat@assays$CCA_RNA@scale.data[as.character(table2$gene),])
s2<-as.matrix(s[,rownames(annotation3)])

pheatmap(s2, show_rownames=T, cluster_rows=T,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)


#chromatin 
DefaultAssay(subset_seurat) <- "peaks"

s<-as.matrix(subset_seurat@assays$peaks@data[tf_chromatin$feature,])
s2<-as.matrix(s[,rownames(annotation3)])

pheatmap(s2, show_rownames=T, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= F,color = rev(color),annotation_col=annotation3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)




library(presto)

pbmc.genes <- wilcoxauc(glio.int, 'CCA_RNA_snn_res.0.2')
dplyr::count(pbmc.genes, group)

library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
#library(tibble)

m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CGP")#我们使用H
#m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
#m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "REACTOME")
#m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP")


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10) 

cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<-as.vector(unlist(cluster0.genes$auc))
names(ranks)<-cluster0.genes$feature
head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>% head(20)

# 显示top20信号通路
ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 20)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色

cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "1") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<-as.vector(unlist(cluster0.genes$auc))
names(ranks)<-cluster0.genes$feature
head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>% head(20)

# 显示top20信号通路
ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 20)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色

cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "2") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<-as.vector(unlist(cluster0.genes$auc))
names(ranks)<-cluster0.genes$feature
head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>% head(20)

# 显示top20信号通路
ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 10)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色

cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "3") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<-as.vector(unlist(cluster0.genes$auc))
names(ranks)<-cluster0.genes$feature
head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>% head(20)

# 显示top20信号通路
ggplot(fgseaResTidy %>% filter(padj < 0.01) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 10)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色




