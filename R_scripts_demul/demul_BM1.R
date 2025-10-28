library('Seurat')
library(Matrix)

data_dir_bc_filtered <- "/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/BM1/GEX_BM1_multi/outs/per_sample_outs/GEX_BM1_multi/count/sample_filtered_feature_bc_matrix/"
data.list <- Read10X(data.dir = data_dir_bc_filtered)

rna_counts <- data.list[['Gene Expression']]
hto_counts <- data.list[["Antibody Capture"]]


# Fraction of zeros in RNA
sparsity_rna <- 1 - (length(rna_counts@x) / (nrow(rna_counts) * ncol(rna_counts)))
sparsity_rna

# Fraction of zeros in HTO
sparsity_hto <- 1 - (length(hto_counts@x) / (nrow(hto_counts) * ncol(hto_counts)))
sparsity_hto

# Keep only cells that exist in both RNA and HTO
joint_bcs <- intersect(colnames(rna_counts), colnames(hto_counts))
rna_counts <- rna_counts[, joint_bcs]
hto_counts <- as.matrix(hto_counts[, joint_bcs])

# Create a Seurat object using RNA counts
seurat_obj <- CreateSeuratObject(counts = rna_counts, project = 'BM1')

# Add HTO data as a new assay
seurat_obj[["HTO"]] <- CreateAssayObject(counts = hto_counts)


# Normalize HTO counts using centered log-ratio (CLR)
seurat_obj <- NormalizeData(seurat_obj, 
                            assay = "HTO", 
                            normalization.method = "CLR",
                            margin = 2)


seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)
table(seurat_obj$HTO_classification.global)


RidgePlot(
  seurat_obj,
  assay = "HTO",
  features = rownames(seurat_obj[["HTO"]]),
  group.by = "HTO_classification.global"
)

#group cells based on their max expressed HTO
Idents(seurat_obj) <- "HTO_maxID"

# Plot each HTO signal across these groups
RidgePlot(
  seurat_obj,
  assay = "HTO",
  features = rownames(seurat_obj[["HTO"]]),
  ncol = 2
)

hto_features<- rownames(seurat_obj[['HTO']])
print(hto_features)

#Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
FeatureScatter(seurat_obj, feature1 = "TotalSeqB-B0301-Hashtag-1", feature2 = "TotalSeqB-B0302-Hashtag-2")


DefaultAssay(seurat_obj) <- "HTO"
hto_names <- rownames(seurat_obj[["HTO"]])
pairs <- combn(hto_names, 2, simplify = FALSE)

#THIS IS SO COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOL 
library(ggplot2)
for (pair in pairs) {
  p <- FeatureScatter(seurat_obj,
                      feature1 = pair[1],
                      feature2 = pair[2],
                      slot = "data") +
    ggtitle(paste(pair[1], "vs", pair[2]))
  print(p)
  readline(prompt = "Press [Enter] to view next plot")
}


#####################################
library(patchwork)

DefaultAssay(seurat_obj) <- "HTO"
hto_names <- rownames(seurat_obj[["HTO"]])
pairs <- combn(hto_names, 2, simplify = FALSE)

plots <- list()

for (pair in pairs) {
  message("Creating ", pair[1], " vs ", pair[2])
  p <- FeatureScatter(seurat_obj,
                      feature1 = pair[1],
                      feature2 = pair[2],
                      slot = "data") +
    ggtitle(paste(pair[1], "vs", pair[2]))
  plots[[paste(pair, collapse = "_vs_")]] <- p
}

# Combine them into a grid layout
combined_plot <- wrap_plots(plots)
combined_plot
