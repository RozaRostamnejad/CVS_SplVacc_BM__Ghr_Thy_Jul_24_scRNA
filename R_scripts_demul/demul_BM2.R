library('Seurat')
library(Matrix)

data_dir_bc_filtered <- "/gpfs/helios/home/rostamne/CVS_SplVacc_BM__Ghr_Thy_Jul_24_scRNA/code/BM2/GEX_BM2_multi/outs/per_sample_outs/GEX_BM2_multi/count/sample_filtered_feature_bc_matrix/"
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

############### Setup Seurat object and add in the HTO data
seurat_obj <- CreateSeuratObject(counts = rna_counts, project = 'BM2')
# Normalize RNA data with log normalization
seurat_obj <- NormalizeData(seurat_obj)
# Find and scale variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Add HTO data as a new assay
seurat_obj[["HTO"]] <- CreateAssayObject(counts = hto_counts)
# Normalize HTO counts using centered log-ratio (CLR)
seurat_obj <- NormalizeData(seurat_obj, 
                            assay = "HTO", 
                            normalization.method = "CLR",
                            margin = 2)


###plots are too long i have to rename the hashtags:

seurat_obj <- RenameAssays(seurat_obj, HTO = "HTO")  # if assay name is long or messy
rownames(seurat_obj)
Assays(seurat_obj)
rownames(seurat_obj[['HTO']])
#this above only change the assay name

#tochange HTOs fill name

############################## HTOdemux:
seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)
############################## 

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

################################
#Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
hto_names <- rownames(seurat_obj[["HTO"]])
short_names <- setNames(paste0("HTO", seq_along(hto_names)), hto_names)
#hto_names[1]
#short_names[hto_names[1]]
#table(seurat_obj$HTO_maxID)
#head(Idents(seurat_obj)) this is why the label still shows the long names

hto_names <- rownames(seurat_obj[["HTO"]])
short_names <- setNames(paste0("HTO", seq_along(hto_names)), hto_names)

FeatureScatter(
  seurat_obj,
  feature1 = hto_names[1],
  feature2 = hto_names[2],
  slot = "data"
) + ggtitle(paste(short_names[hto_names[1]], "vs", short_names[hto_names[2]]))

################################################################

DefaultAssay(seurat_obj) <- "HTO"
hto_names <- rownames(seurat_obj[["HTO"]])
pairs <- combn(hto_names, 2, simplify = FALSE)

#THIS IS SO COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOL for vewing only one. 
library(ggplot2)
for (pair in pairs) {
  p <- FeatureScatter(seurat_obj,
                      feature1 = pair[1],
                      feature2 = pair[2],
                      slot = "data") +
    ggtitle(short_names[pair[1]], "vs", short_names[pair[2]])
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
    ggtitle(paste(short_names[pair[1]], "vs", short_names[pair[2]]))
  plots[[paste(pair, collapse = "_vs_")]] <- p
}

# Combine them into a grid layout
combined_plot <- wrap_plots(plots)
combined_plot

##########################################
#Compare number of UMIs for singlets, doublets and negative cells

Idents(seurat_obj) #each cell barcode assigned to a hashtag
DefaultAssay(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "HTO_classification.global"
VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


##########################################
#Generate a two dimensional tSNE embedding for HTOs. Here we are grouping cells by singlets and doublets for simplicity.
# First, we will remove negative cells from the object
seurat_obj.subset <- subset(seurat_obj, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(seurat_obj.subset) <- "HTO"
seurat_obj.subset <- ScaleData(seurat_obj.subset, features = rownames(seurat_obj.subset),
                                 verbose = FALSE)

#got error, hence removing duplicates: 
hto_matrix <- t(as.matrix(GetAssayData(seurat_obj.subset, assay = "HTO", slot = "data")))
hto_matrix_unique <- unique(hto_matrix)
cells_to_keep <- rownames(hto_matrix_unique)
seurat_obj.subset <- subset(seurat_obj.subset, cells = cells_to_keep)


seurat_obj.subset <- RunPCA(seurat_obj.subset, features = rownames(seurat_obj.subset), approx = FALSE)
seurat_obj.subset <- RunTSNE(seurat_obj.subset, dims = 1:4, perplexity = 100)
DimPlot(seurat_obj.subset)


##########################################
#Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(seurat_obj, assay = "HTO", ncells = 5000)



##########################################
#Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.

