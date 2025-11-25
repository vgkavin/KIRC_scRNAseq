library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(pheatmap)


#load raw count matrix files
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

View(kirc.1.seurat_list$GSM4819725_SI_18854@meta.data)

#separate all seurats 
list2env(kirc.1.seurat_list, envir = .GlobalEnv)

#merge them into a single seurat object
kirc.1.obj <- merge(GSM4819725_SI_18854, y = c(GSM4819727_SI_18855, GSM4819729_SI_19703, GSM4819734_SI_22368, 
                                             GSM4819736_SI_22604, GSM4819737_SI_23459, GSM4819738_SI_23843),
                  project = "KIRC1")


#initial downstream workflow to check for batch effects

#Calculate %mitochondrial genes and add to metadata
kirc.1.obj$MT <- PercentageFeatureSet(kirc.1.obj, pattern = "^MT-") 

#QC
#Filtering sets based on criteria provided by author
kirc.1.obj <- subset(kirc.1.obj, subset = nFeature_RNA >= 300 & MT <= 25)

#normalize data
kirc.1.obj <- NormalizeData(kirc.1.obj)

#find differentally expressed genes to segregate clusters
kirc.1.obj <- FindVariableFeatures(kirc.1.obj, nfeatures = 3000, verbose = TRUE)

#scaling to reduce unwanted noise
kirc.1.genes <- rownames(kirc.1.obj)
kirc.1.obj <- ScaleData(kirc.1.obj, features = kirc.1.genes)

#Linear Dimensionality reduction by identifying Princple components 
kirc.1.DEG <- VariableFeatures(kirc.1.obj)
kirc.1.obj <- RunPCA(kirc.1.obj, features = kirc.1.DEG)
ElbowPlot(kirc.1.obj)

#identify genes with similar expression patterns and make cell clusters at author provided resolution
kirc.1.obj <- FindNeighbors(kirc.1.obj)
kirc.1.obj <- FindClusters(kirc.1.obj, resolution = 0.2)
Idents(kirc.1.obj) <- kirc.1.obj$RNA_snn_res.0.2

#Non-linear dimentionality reduction using UMAP
kirc.1.obj <- RunUMAP(kirc.1.obj, dims = 1:20)


#plot UMAPs and comparing patient level and cluster level variation to test for batch effects
kirc.1.patient <- DimPlot(kirc.1.obj, reduction = "umap", group.by = "Sample_Id")
kirc.1.cluster <- DimPlot(kirc.1.obj, reduction = "umap", group.by = "RNA_snn_res.0.2")

kirc.1.unint.umap <- grid.arrange(kirc.1.patient, kirc.1.cluster, ncol = 2)
ggsave(kirc.1.unint.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KICR_1_Unintegrated.png")


#Btach effect seen when grouping UMAP by patients- Need to integrate

#integration(batch correction) using CCA in seurat


#use the initally created list of seurats for initial QC, normalization and finding DEGs
for(i in 1:length(kirc.1.seurat_list)){
  kirc.1.seurat_list[[i]][["MT"]] <- PercentageFeatureSet(kirc.1.seurat_list[[i]], pattern = "^MT-")
  kirc.1.seurat_list[[i]] <- subset(kirc.1.seurat_list[[i]], subset = nFeature_RNA >= 300 & MT <= 25)
  kirc.1.seurat_list[[i]] <- NormalizeData(kirc.1.seurat_list[[i]], verbose = TRUE)
  kirc.1.seurat_list[[i]] <- FindVariableFeatures(kirc.1.seurat_list[[i]], nfeatures = 3000, verbose = TRUE )
}

# select integration features 
kirc.1.features <- SelectIntegrationFeatures(object.list = kirc.1.seurat_list)

#Find integration anchors
kirc.1.anchors <- FindIntegrationAnchors(object = kirc.1.seurat_list,
                                  anchor.features = kirc.1.features)

#integrate all data into a single seurat object using the anchors 
kirc.1.int <- IntegrateData(anchorset = kirc.1.anchors)

#similar workflow to previous
# Scaling
kirc.1.int.genes <- rownames(kirc.1.int)
kirc.1.int <- ScaleData(kirc.1.int, features = kirc.1.int.genes)

# Linear Dimensionality reduction- Princple component analysis
kirc.1.DEG.int <- VariableFeatures(kirc.1.int)
kirc.1.int <- RunPCA(kirc.1.int, features = kirc.1.DEG.int)

#to view the top 5 positive and negative sources of heterogeniety in principle components
print(kirc.1.int[["pca"]], dims = 1:20, nfeatures = 5) 

#to visualize the sources of heterogeinety in each cell.
DimHeatmap(kirc.1.int, dims = 1, cells = 500, balanced = TRUE) 

#to determine the dimensionality and filter PCs based on variance.
ElbowPlot(kirc.1.int) 


#Clustering based on gene expression similarity
kirc.1.int <- FindNeighbors(kirc.1.int, dims = 1:50)
kirc.1.int <- FindClusters(kirc.1.int, resolution = 0.2)
Idents(kirc.1.int) <- "integrated_snn_res.0.2"

#non linear dimensionality reduction using umap
kirc.1.int <- RunUMAP(kirc.1.int, dims = 1:50)

#plot UMAPs and comparing patient level and cluster level variation to comfirm integration
kirc.1.int.patient.umap <- DimPlot(kirc.1.int, reduction = "umap", group.by = "Sample_Id")
kirc.1.int.cluster.umap <- DimPlot(kirc.1.int, reduction = "umap", group.by = "integrated_snn_res.0.2")

kirc.1.int.umap <- grid.arrange(kirc.1.int.patient.umap, kirc.1.int.cluster.umap, ncol = 2)
ggsave(kirc.1.int.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KIRC1_integrated.png")

###integration successful- clustering based on gene expression patterns and no batch effect found!


###Annotation
#Join assay layers
kirc.1.int <- JoinLayers(kirc.1.int, assay = "RNA")

#find mmajor marker genes for each clusters
kirc.1.int.markers <- FindAllMarkers(
  kirc.1.int,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#subset top 10 marker genes to identify cell types in cluster
kirc.1.int.markers.sig <- (kirc.1.int.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))
write.csv(kirc.1.int.markers, "KIRC_1_Markers.csv")

#create an object with identified celltypes
kirc.1.cluster.ids <- c(
  "0" = "Tumor",
  "1" = "Macrophages",
  "2" = "Vascular Endothelial Cells",
  "3" = "Muscle",
  "4" = "T-Cells",
  "5" = "Pericytes/Fibroblasts",
  "6" = "Proliferating Cells",
  "7" = "Dendritic Cells",
  "8" = "Venular Endothelial Cells",
  "9" = "B Cells",
  "10" = "Mast Cells",
  "11" = "Erythroid Cells"
)

#match annotated celltypes to cluster idents and add to metadata
kirc.1.int$annotation <- factor( x = kirc.1.int$integrated_snn_res.0.2,
                                 levels = names(kirc.1.cluster.ids),
                                 labels = kirc.1.cluster.ids)

#visualize
kirc.1.int.ann.umap <- DimPlot(kirc.1.int, reduction = "umap", group.by = "annotation")+ ggtitle("GSE159115")
ggsave(kirc.1.int.ann.umap, dpi = 300, height = 6, width = 12, filename = "UMAP_KIRC1_annotated.png")

#save final dataset with all analyses and metadata
saveRDS(kirc.1.int, "KIRC_1_Int.rds")
kirc.1.int <- readRDS("KIRC_1_Int.rds")

#plot expression of genes of interest in UMAP
gene_list <- c("C1orf174", "TMEM161B", "ZNF808", 
               "TMEM9B", "ZNF683", "ZNF860", "TMEM31", 
               "C2orf49", "C5orf22")

#fix raw assay data as defaukt to map gene counts
DefaultAssay(kirc.1.int) <- "RNA"

#visualize expression of target genes in UMAP
kirc.1.GE <- FeaturePlot(kirc.1.int, features = gene_list, reduction = "umap")
kirc.1.GE<- wrap_plots(kirc.1.GE) + 
  plot_annotation(title = "Target gene expression in GSE159115") &
  theme(plot.title = element_text(
      size = 20, hjust = 0.5, face = "bold"     
    )
  )

ggsave(kirc.1.GE, dpi = 300, height = 12, width = 12, filename = "KIRC_1_Gene_expression.png" )


