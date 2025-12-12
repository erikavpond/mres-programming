          ################################
          ###        GSE163974         ###
          ################################
# This dataset contains 6 human skin samples: 3 from keloids, and 3 from normal scar tissue
# My goal is to compare the two groups and identify differences in gene expression and presence of cell types
# Potentially looking into pathways that are upregulated/downregulated to explore functional disruption

setwd("~/scratch/module_4_project") #sets working directory appropriate to my HPC setup

library(Seurat)
library(dplyr)
library(ggplot2)



          ################################
          ###      Loading Data        ###
          ################################

# Defining which folder within working directory holds the sample data
data_dir <- "GSE163974/" 
# Creating a list of sample folder names within the data directory
sample_dirs <- list.dirs(data_dir, 
                         full.names = TRUE,
                         recursive = FALSE)

# Creating function called 'load_sample'
load_sample <- function(path) {
  counts <- Read10X(path) #standard method for reading 10x scRNA data
  
  obj <- CreateSeuratObject(counts, project = basename(path))
  if (grepl("^K", basename(path), ignore.case = TRUE)) {
    obj$condition <- "keloid"
  } else {
    obj$condition <- "normal"}
  return(obj)
} #Loads an individual sample, and assigns condition based on known file naming style

# Loads all samples in GSE163974 set by applying paths generated in sample_dirs as input
samples <- lapply(sample_dirs, load_sample) #each sample is made into a seurat object, and added to a list called 'samples'
names(samples) <- basename(sample_dirs) #each object is named according to original folder name e.g. KF1_matrix



          ################################
          ###            QC            ###
          ################################

#Percent mitochondrial transcripts and basic metrics
qual_control <- function(sobj) {
  #Calculates percent.mt
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  #Summary per object
  cat("SAMPLE NAME: ", sobj@project.name, "\n") #\n for new line improving readability
  cat("nCount_RNA: ")
  print(summary(sobj$nCount_RNA))
  cat("nFeature_RNA: ")
  print(summary(sobj$nFeature_RNA))
  cat("% mitochondrial: ")
  print(summary(sobj$percent.mt))
}
samples_qc <- lapply(samples, qual_control) #Runs function and applies to samples

#Ensuring percent.mt added to samples
samples <- lapply(samples, function(sobj) {
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  return(sobj)
})

#Gives a pdf for each sample, with 3 plots for QC in each
for (i in seq_along(samples_qc)) {
  sample_name <- samples[[i]]@project.name
  pdf(paste0("results/QC_", sample_name, "_plots.pdf"), width = 50, height= 30)
  p1 <- VlnPlot(samples[[i]], 
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.005)
  print(p1)
  p2 <- FeatureScatter(samples[[i]],
                       "nCount_RNA",
                       "percent.mt")
  print(p2)
  p3 <- FeatureScatter(samples[[i]],
                       "nCount_RNA",
                       "nFeature_RNA")
  print(p3)
  dev.off()
}



          ################################
          ###        Normalise         ###
          ################################

# Individually normalising each Seurat object within the list
# WARNING: This step takes a long time to complete
# PLEASE NOTE: install glmGamPoi to speed up SCTransform process. This package was not compatible with the R v4 on the Aston HPC
samples <- lapply(samples, function(sobj) {
  sobj <- subset(sobj, subset= nFeature_RNA>200 & percent.mt<15)
  sobj <- SCTransform(sobj, verbose = FALSE) #SCTransform() used rather than NormalizeData() as needed for integration in next step
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  return(sobj)
})

saveRDS(samples, file = "results/seurat_objects.rds")


          ################################
          ###        Integration       ###
          ################################

# Selecting the features (genes) for integration
features <- SelectIntegrationFeatures(object.list = samples)

#Prepare for integration
samples <- PrepSCTIntegration(object.list = samples, 
                              anchor.features = features)

# Find anchors for integration
# WARNING: This step takes a long time to complete
anchors <- FindIntegrationAnchors(object.list = samples,
                                  normalization.method = "SCT",
                                  anchor.features = features)
# Integrate
integrated_samples <- IntegrateData(anchorset = anchors,
                                    normalization.method = "SCT")

# Ensuring Seurat uses the integrated assay
DefaultAssay(integrated_samples) <- "integrated"



          ################################
          ###         PCA, UMAP,       ###
          ###  Clustering, Neighbours  ###
          ################################

# Dimensional reduction
integrated_samples <- RunPCA(integrated_samples, features = VariableFeatures(object = integrated_samples), verbose = FALSE)

pdf(file = "results/PCA_elbowplot.pdf", width = 10, height = 10)
ElbowPlot(integrated_samples)
dev.off()

# Runs UMAP on the first 15 principle components
integrated_samples <- RunUMAP(integrated_samples, dims = 1:15, verbose = FALSE)

# Finding nearest neighbours
integrated_samples <- FindNeighbors(integrated_samples, dims = 1:15, verbose = FALSE)

# Clustering cells
integrated_samples <- FindClusters(integrated_samples, resolution = 0.25) #reduced resolution for time efficiency when finding markers

saveRDS(integrated_samples, file = "results/integrated_GSE163974.rds")

          ################################
          ###       Visualisation      ###
          ################################

# Visualising UMAP by 'condition', which is keloid vs normal here
pdf("results/UMAP_conditions.pdf", width=10, height=10)
DimPlot(integrated_samples, group.by = "condition", reduction = "umap")
dev.off()

#Visualising UMAP by cluster
pdf("results/UMAP_clusters.pdf", width=10, height=10)
DimPlot(integrated_samples, group.by = "seurat_clusters", reduction = "umap", label = TRUE)
dev.off()



          ################################
          ###        Markers/ DE       ###
          ################################

#Switching back to RNA assay
DefaultAssay(integrated_samples) <- "RNA"

# IMPORTANT: RNA assay must be normalised and scaled manually 
integrated_samples <- NormalizeData(integrated_samples)
integrated_samples <- FindVariableFeatures(integrated_samples)
integrated_samples <- ScaleData(integrated_samples)

#Keloid vs normal scar DE
Idents(integrated_samples) <- integrated_samples$condition
# WARNING: This step takes a long time to complete
kvn_markers <- FindMarkers(integrated_samples,
                           ident.1 = "keloid",
                           ident.2 = "normal",
                           group.by = "condition",
                           min.pct = 0.2, #expression in at least 20% of cells in group
                           logfc.threshold = 0.5)
write.csv(kvn_markers, "results/DE_keloid_vs_normal.csv")

#Top 10 upregulated genes in keloid samples
top10_tabl <- kvn_markers |>
  filter(p_val_adj < 0.05) |> #only significant genes kept
  arrange(desc(avg_log2FC)) |> #highest log2fc first
  head(10)  #top10

top10 <- rownames(top10_tabl)
top3 <- head(top10, 3)
  
#Violin plot and feature plot on UMAP for top 10 genes
pdf("results/top_upregulated_keloid.pdf", width=40, height=30)
p4 <- VlnPlot(merged_samples, features = top10, group.by = "condition", pt.size = 0.01, ncol = 5)
print(p4)
DefaultAssay(integrated_samples) <- "integrated"
p5 <- FeaturePlot(integrated_samples, features = top10)
print(p5)
dev.off()
