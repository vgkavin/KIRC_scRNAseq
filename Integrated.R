library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(pheatmap)


#load KIRC 1 raw datasets and make a seurat obj- if list from individual analysis is not in global environment
h5_files <- list.files("GSE159115_RAW", pattern = "\\.h5$", full.names = TRUE)

#create an empty list to convert all matrices into seurat and store for ease of editing
kirc.1.seurat_list <- list()

for (file in h5_files) {
  filename <- basename(file)
  
  # Extract text BEFORE "_filtered"
  sample_name <- sub("_filtered.*", "", filename)
  
  kirc.1.seurat_list[[sample_name]] <- CreateSeuratObject(
    Read10X_h5(file),
    project = sample_name
  )
}


#edit metadata for all seurat objects to indicate tissue location, study info(GEO accession, Sample IDs)
for (sample in names(kirc.1.seurat_list)) {
  obj <- kirc.1.seurat_list[[sample]]
  obj$Tissue <- "Primary tumor"
  obj$Study <- "GSE159115"
  obj$Sample_Id <- sub("_.*", "", obj$orig.ident)
  kirc.1.seurat_list[[sample]] <- obj
}



#load KIRC 2 raw datasets and make a seurat obj- if list from individual analysis is not in global environment
#load raw count matrix files
raw_files <- list.files("GSE156632_RAW", pattern = "\\.csv$", full.names = TRUE)


raw_data_list <- list()

for (file in raw_files) {
  raw_data_list[[sub("\\.csv$", "", basename(file))]] <- read.csv(file, stringsAsFactors = FALSE)
}


for (m in names(raw_data_list)) {
  matrix <- raw_data_list[[m]]       
  counts <- matrix[!duplicated(matrix$Symbol), ]
  rownames(counts) <- counts$Symbol
  counts <- counts[, -(1:2)]
  counts <- as.matrix(counts)
  raw_data_list[[m]] <- counts        
}


kirc.2.seurat_list <- list()

for (m in names(raw_data_list)) {
  counts <- raw_data_list[[m]]
  
  kirc.2.seurat_list[[m]] <- CreateSeuratObject(
    counts = counts,
    project = m
  )
}


for (sample in names(kirc.2.seurat_list)) {
  obj <- kirc.2.seurat_list[[sample]]
  obj$Tissue <- "Primary tumor"
  obj$Study <- "GSE156632"
  obj$Sample_Id <- sub("_.*", "", obj$orig.ident)
  kirc.2.seurat_list[[sample]] <- obj
}

#combine datasets from both studies into a single list
kirc.obj.seurat_list <- c(kirc.1.seurat_list, kirc.2.seurat_list)



#integration and QC filtering workflow
for(i in 1:length(kirc.obj.seurat_list)){
  kirc.obj.seurat_list[[i]][["MT"]] <- PercentageFeatureSet(kirc.obj.seurat_list[[i]], pattern = "^MT-")
  kirc.obj.seurat_list[[i]] <- subset(kirc.obj.seurat_list[[i]], subset = nFeature_RNA >= 300 & MT <= 25)
  kirc.obj.seurat_list[[i]] <- NormalizeData(kirc.obj.seurat_list[[i]], verbose = TRUE)
  kirc.obj.seurat_list[[i]] <- FindVariableFeatures(kirc.obj.seurat_list[[i]], nfeatures = 3000, verbose = TRUE )
}

# select integration features 
kirc.features <- SelectIntegrationFeatures(object.list = kirc.obj.seurat_list)

#Find integration anchors
kirc.anchors <- FindIntegrationAnchors(object = kirc.obj.seurat_list,
                                  anchor.features = kirc.features)

#integrate all data into a single seurat object using the anchors 
kirc.int <- IntegrateData(anchorset = kirc.anchors)

#scaling to reduce unwanted noise
kirc.int.genes <- rownames(kirc.int)
kirc.int <- ScaleData(kirc.int, features = kirc.int.genes)

#Linear Dimensionality reduction by identifying Princple components 
kirc.int.DEG <- VariableFeatures(kirc.int)
kirc.int <- RunPCA(kirc.int, features = kirc.int.DEG)
ElbowPlot(kirc.int)

#identify genes with similar expression patterns and make cell clusters at author provided resolution
kirc.int <- FindNeighbors(kirc.int)
kirc.int <- FindClusters(kirc.int, resolution = c(0.05, 0.1, 0.13, 0.15, 0.17, 0.2))
Idents(kirc.int) <- kirc.int$RNA_snn_res.0.2

#Non-linear dimentionality reduction using UMAP
kirc.int <- RunUMAP(kirc.int, dims = 1:20)

UMAP_kirc.int.cluster.0.05 <- DimPlot(kirc.int, reduction = "umap", group.by = "integrated_snn_res.0.05")
UMAP_kirc.int.cluster.0.1 <- DimPlot(kirc.int, reduction = "umap", group.by = "integrated_snn_res.0.1")
UMAP_kirc.int.cluster.0.15 <- DimPlot(kirc.int, reduction = "umap", group.by = "integrated_snn_res.0.15")
UMAP_kirc.int.cluster.0.2 <- DimPlot(kirc.int, reduction = "umap", group.by = "integrated_snn_res.0.2")
UMAP_kirc.int.cluster.0.13 <- DimPlot(kirc.int, reduction = "umap", group.by = "integrated_snn_res.0.13")
UMAP_kirc.int.patient <- DimPlot(kirc.int, reduction = "umap", group.by = "Sample_Id")
UMAP_kirc.int.study <- DimPlot(kirc.int, reduction = "umap", group.by = "Study")



grid.arrange(UMAP_kirc.int.cluster.0.05, UMAP_kirc.int.cluster.0.1, UMAP_kirc.int.cluster.0.13, UMAP_kirc.int.cluster.0.15, 
             UMAP_kirc.int.cluster.0.2, ncol = 3, nrow = 2)

grid.arrange(UMAP_kirc.int.patient, UMAP_kirc.int.study, UMAP_kirc.int.cluster.0.1, ncol = 2)

Idents(kirc.int) <- kirc.int$integrated_snn_res.0.13

#Join assay layers
kirc.int <- JoinLayers(kirc.int, assay = "RNA")

#find mmajor marker genes for each clusters
kirc.int.markers <- FindAllMarkers(
  kirc.int,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

kirc.int.markers.sig <- (kirc.int.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))
write.csv(kirc.int.markers.sig, "KIRC_Markers.csv")



kirc.int.cluster.ids <- c(
  "0" = "Tumor",
  "1" = "Macrophages",
  "2" = "Vascular endothelial Cells",
  "3" = "Myofibroblasts",
  "4" = "T-Cells",
  "5" = "Fibroblasts",
  "6" = "Monocytes",
  "7" = "Venular endothelial Cells",
  "8" = "Proliferating Cells",
  "9" = "Plasma Cells"
)

#match annotated celltypes to cluster idents and add to metadata
kirc.int$annotation <- factor( x = kirc.int$integrated_snn_res.0.13,
                                 levels = names(kirc.int.cluster.ids),
                                 labels = kirc.int.cluster.ids)

#visualize
kirc.int.ann.umap <- DimPlot(kirc.int, reduction = "umap", group.by = "annotation")+ ggtitle("Integrated UMAP")
ggsave(kirc.int.ann.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KIRC_annotated.png")

#save final results 
saveRDS(kirc.int, "KIRC.rds")
kirc.int <- readRDS("KIRC.rds")

#fix raw assay data as defaukt to map gene counts
DefaultAssay(kirc.int) <- "RNA"

#list of target genes
gene_list <- c("C1orf174", "TMEM161B", "ZNF808", 
               "TMEM9B", "ZNF683", "ZNF860", "TMEM31", 
               "C2orf49", "C5orf22")

#visualize gene expression in UMAP, splitting data from both studies side by side
kirc.plot.list <- list()

kirc.int$Study <- factor(kirc.int$Study,
                         levels = c("GSE159115", "GSE156632"))

for (g in seq_along(gene_list)) {
  gene <- gene_list[g]
  plot <- FeaturePlot(kirc.int,
                      features = gene,
                      split.by = "Study",
                      reduction = "umap")
  kirc.plot.list[[g]] <- plot
}

names(kirc.plot.list) <- gene_list


#save images
for (gene in names(kirc.plot.list)) {
  plot <- kirc.plot.list[[gene]]
  ggsave(plot = plot, height = 6, width = 12, dpi = 300,   filename = paste0(gene, "_study_comparisson.png"))
}

