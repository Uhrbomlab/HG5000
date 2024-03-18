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


load("/disk1/xilu/glioblastoma_single_cell/1.RNA/1.RNA_seq_seurat/glio.int.all.peaks.chromvar.rdata")

p1 <- DimPlot(glio.int, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(glio.int, reduction = "umap", group.by = "CCA_RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(glio.int, reduction = "wnn.umap", group.by = "CCA_RNA_snn_res.0.2",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

col<-c("#B3CDE3","#FCCDE5","#B2DF8A","#E6AB02","#FDBF6F")

getwd()



DefaultAssay(glio.int) <- "CCA_RNA"
clusters_of_interest <- c(0,1)
subset_seurat2 <- subset(glio.int, idents = clusters_of_interest)
DimPlot(subset_seurat2, reduction = "umap", group.by = "CCA_RNA_snn_res.0.1",cols = c("#FDC086","#FCCDE5")) + ggtitle("louvain_0.2")

DimPlot(subset_seurat2, reduction = "umap", group.by = "CCA_RNA_snn_res.0.1",cols = col) + ggtitle("louvain_0.2")


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

DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.3",cols = c("#941B2A","#AF4137","#A37660")) + ggtitle("louvain_0.3")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.1")

pdf("core.sub.clsuters.pdf")
DimPlot(subset_seurat,reduction="umap",group.by="peaks_snn_res.0.3",cols = c("#941B2A","#AF4137","#A37660")) + ggtitle("louvain_0.3")
DimPlot(subset_seurat, reduction = "rna.umap", group.by = "CCA_RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.1")

dev.off()


DimPlot(subset_seurat, reduction = "rna.umap", group.by = "peaks_snn_res.0.3",cols = col) + ggtitle("louvain_0.1")






save(subset_seurat,file="core_cells.rdata")

FeaturePlot(subset_seurat, features = 'NT5E',reduction="umap")
VlnPlot(subset_seurat, features = "NT5E",group.by = "peaks_snn_res.0.3",pt.size=0,cols = c("#941B2A","#AF4137","#A37660"))

FeaturePlot(subset_seurat, features = 'OSMR',reduction="umap")
VlnPlot(subset_seurat, features = "OSMR",group.by = "peaks_snn_res.0.3",pt.size=0,cols = c("#941B2A","#AF4137","#A37660"))

FeaturePlot(subset_seurat, features = 'NRP1',reduction="umap")
VlnPlot(subset_seurat, features = "NRP1",group.by = "peaks_snn_res.0.3",pt.size=0,cols = c("#941B2A","#AF4137","#A37660"))


FeaturePlot(subset_seurat, features = 'NF1',reduction="umap")
VlnPlot(subset_seurat, features = "NF1",group.by = "peaks_snn_res.0.3",pt.size=0)

DefaultAssay(subset_seurat) <- "CCA_RNA"
VlnPlot(subset_seurat, features = "NF1",group.by = "peaks_snn_res.0.3",cols=col,pt.size=0)


