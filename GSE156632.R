library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(dplyr)


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


list2env(kirc.2.seurat_list, envir = .GlobalEnv)
ls()

kirc.2.obj <- merge(GSM4735364_RCC1t, y = c(GSM4735366_RCC2t, GSM4735368_RCC3t, 
                                            GSM4735370_RCC4t, GSM4735372_RCC5t, 
                                            GSM4735374_RCC6t),
                    project = "KIRC2")

#initial downstream workflow to check for batch effects

#Calculate %mitochondrial genes and add to metadata
kirc.2.obj$MT <- PercentageFeatureSet(kirc.2.obj, pattern = "^MT-") 

#QC
#Filtering sets based on criteria provided by author
kirc.2.obj <- subset(kirc.2.obj, subset = nFeature_RNA >= 300 & MT <= 25)

#normalize data
kirc.2.obj <- NormalizeData(kirc.2.obj)

#find differentally expressed genes to segregate clusters
kirc.2.obj <- FindVariableFeatures(kirc.2.obj, nfeatures = 3000, verbose = TRUE)

#scaling to reduce unwanted noise
kirc.2.genes <- rownames(kirc.2.obj)
kirc.2.obj <- ScaleData(kirc.2.obj, features = kirc.2.genes)

#Linear Dimensionality reduction by identifying Princple components 
kirc.2.DEG <- VariableFeatures(kirc.2.obj)
kirc.2.obj <- RunPCA(kirc.2.obj, features = kirc.2.DEG)
ElbowPlot(kirc.2.obj)

#identify genes with similar expression patterns and make cell clusters at author provided resolution
kirc.2.obj <- FindNeighbors(kirc.2.obj)
kirc.2.obj <- FindClusters(kirc.2.obj, resolution = 0.2)
Idents(kirc.2.obj) <- kirc.2.obj$RNA_snn_res.0.2

#Non-linear dimentionality reduction using UMAP
kirc.2.obj <- RunUMAP(kirc.2.obj, dims = 1:20)

#plot UMAPs and comparing patient level and cluster level variation to test for batch effects
kirc.2.patient <- DimPlot(kirc.2.obj, reduction = "umap", group.by = "Sample_Id")
kirc.2.cluster <- DimPlot(kirc.2.obj, reduction = "umap", group.by = "RNA_snn_res.0.2")

kirc.2.unint.umap <- grid.arrange(kirc.2.patient, kirc.2.cluster, ncol = 2)
ggsave(kirc.2.unint.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KIRC_2_Unintegrated.png")

#Btach effect seen when grouping UMAP by patients- Need to integrate

#integration(batch correction) using CCA in seurat


#use the initally created list of seurats for initial QC, normalization and finding DEGs
for(i in 1:length(kirc.2.seurat_list)){
  kirc.2.seurat_list[[i]][["MT"]] <- PercentageFeatureSet(kirc.2.seurat_list[[i]], pattern = "^MT-")
  kirc.2.seurat_list[[i]] <- subset(kirc.2.seurat_list[[i]], subset = nFeature_RNA >= 300 & MT <= 25)
  kirc.2.seurat_list[[i]] <- NormalizeData(kirc.2.seurat_list[[i]], verbose = TRUE)
  kirc.2.seurat_list[[i]] <- FindVariableFeatures(kirc.2.seurat_list[[i]], nfeatures = 3000, verbose = TRUE )
}

# select integration features 
kirc.2.features <- SelectIntegrationFeatures(object.list = kirc.2.seurat_list)

#Find integration anchors
kirc.2.anchors <- FindIntegrationAnchors(object = kirc.2.seurat_list,
                                         anchor.features = kirc.2.features)

#integrate all data into a single seurat object using the anchors 
kirc.2.int <- IntegrateData(anchorset = kirc.2.anchors)

#similar workflow to previous
# Scaling
kirc.2.int.genes <- rownames(kirc.2.int)
kirc.2.int <- ScaleData(kirc.2.int, features = kirc.2.int.genes)

# Linear Dimensionality reduction- Princple component analysis
kirc.2.int.DEG <- VariableFeatures(kirc.2.int)
kirc.2.int <- RunPCA(kirc.2.int, features = kirc.2.int.DEG)

#to view the top 5 positive and negative sources of heterogeniety in all 5 principle components
print(kirc.2.int[["pca"]], dims = 1:20, nfeatures = 5) 

#to visualize the sources of heterogeinety in each cell.
DimHeatmap(kirc.2.int, dims = 1, cells = 500, balanced = TRUE) 
VizDimLoadings(kirc.2.int, dims = 1:2, reduction = "pca")

#to determine the dimensionality and filter PCs based on variance.
ElbowPlot(kirc.2.int) 

#Clustering

#Clustering based on gene expression similarity
kirc.2.int <- FindNeighbors(kirc.2.int, dims = 1:50)
kirc.2.int <- FindClusters(kirc.2.int, resolution = c(0.1, 0.13, 0.15, 0.2)) #to identify best clustering
Idents(kirc.2.int) <- kirc.2.int$integrated_snn_res.0.13


#non linear dimensionality reduction using umap
kirc.2.int <- RunUMAP(kirc.2.int, dims = 1:50)



#plot UMAPs and comparing patient level and cluster level variation to comfirm integration
UMAP_kirc.2.int.patient <- DimPlot(kirc.2.int, reduction = "umap", group.by = "Sample_Id")
UMAP_kirc.2.int.cluster <- DimPlot(kirc.2.int, reduction = "umap", group.by = "integrated_snn_res.0.13")

kirc.2.int.umap <- grid.arrange(UMAP_kirc.2.int.patient, UMAP_kirc.2.int.cluster, ncol = 2)
ggsave(kirc.2.int.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KIRC_2_integrated.png")

#integration successful- clustering based on gene expression patterns and no batch effect found!


###Annotation

#Join assay layers
kirc.2.int <- JoinLayers(kirc.2.int, assay = "RNA")

#find mmajor marker genes for each clusters
kirc.2.int.markers <- FindAllMarkers(
  kirc.2.int,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#subset top 10 marker genes to identify cell types in cluster
kirc.2.int.markers.sig <- (kirc.2.int.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))
write.csv(kirc.2.int.markers.sig, "KIRC_2_Markers.csv")


kirc.2.cluster.ids <- c(
  "0" = "Endothelial Cells",
  "1" = "Tumor",
  "2" = "Macrophages",
  "3" = "Fibroblasts",
  "4" = "T-Cells",
  "5" = "Pericytes",
  "6" = "Monocytes",
  "7" = "Proliferating Cells",
  "8" = "Myeloid Cells",
  "9" = "Mast Cells"
)

#match annotated celltypes to cluster idents and add to metadata
kirc.2.int$annotation <- factor( x = kirc.2.int$integrated_snn_res.0.13,
                                 levels = names(kirc.2.cluster.ids),
                                 labels = kirc.2.cluster.ids)

#visualize
kirc.2.int.ann.umap <- DimPlot(kirc.2.int, reduction = "umap", group.by = "annotation")+ ggtitle("GSE156632")
ggsave(kirc.2.int.ann.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KIRC2_annotated.png")

saveRDS(kirc.2.int, "KIRC_2_Int.rds")
kirc.2.int <- readRDS("KIRC_2_Int.rds")

View(kirc.2.int@meta.data)

#plot expression of genes of interest in UMAP
gene_list <- c("C1orf174", "TMEM161B", "ZNF808", 
               "TMEM9B", "ZNF683", "ZNF860", "TMEM31", 
               "C2orf49", "C5orf22")

#fix raw assay data as defaukt to map gene counts
DefaultAssay(kirc.2.int) <- "RNA"

#visualize expression of target genes in UMAP
kirc.2.GE <- FeaturePlot(kirc.2.int, features = gene_list, reduction = "umap")
kirc.2.GE<- wrap_plots(kirc.2.GE) + 
  plot_annotation(title = "Target gene expression in GSE156632") &
  theme(plot.title = element_text(
    size = 20, hjust = 0.5, face = "bold"     
  )
  )



ggsave(kirc.2.GE, dpi = 300, height = 12, width = 12, filename = "KIRC2_Gene_expression.png" )

