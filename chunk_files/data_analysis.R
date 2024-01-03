# load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(umap)


# .RDS format
rds_obj <- readRDS("D:\\Finished projects\\Scrna_data_analysis\\datasets\\ependymal_cells.rds")

# Display structure of the object
str(rds_obj)

# Display the first few rows of the object (assuming it's a data frame)
rds_obj@raw.data

rds_obj@data

rds_obj@meta.data

# Access expression data
expression_data <- rds_obj@scale.data




# 10X CellRanger .HDF5 format 
hdf5_obj <- Read10X_h5(filename = "D:\\Finished projects\\Scrna_data_analysis\\datasets\\20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)
seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)

seurat_hdf5@meta.data
seurat_hdf5@project.name
seurat_hdf5@assays
seurat_hdf5@commands
seurat_hdf5@active.assay
seurat_hdf5@active.ident
seurat_hdf5@reductions
seurat_hdf5@assays$RNA$counts


# .mtx file
# Use double backslashes
mtx_obj <- ReadMtx(mtx = "D:\\Finished projects\\Scrna_data_analysis\\datasets\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix\\raw_feature_bc_matrix\\matrix.mtx.gz",
                   features = "D:\\Finished projects\\Scrna_data_analysis\\datasets\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix\\raw_feature_bc_matrix\\features.tsv.gz",
                   cells = "D:\\Finished projects\\Scrna_data_analysis\\datasets\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix\\raw_feature_bc_matrix\\barcodes.tsv.gz")

seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

seurat_mtx@meta.data
a <- seurat_mtx@assays$RNA$counts

a@Dimnames

a@i

a@Dim
a@x







# 10X CellRanger .HDF5 format 
nsclc.sparse.m<- Read10X_h5(filename ="D:\\Finished projects\\Scrna_data_analysis\\20k_NSCLC_raw_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)


str(nsclc.sparse.m)
cts <-  nsclc.sparse.m$`Gene Expression`


nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj

#29552 features across 42081 samples within 1 assay 


# 1. QC -------
View(nsclc.seurat.obj@meta.data)


# Access the RNA assay data (gene expression matrix)
gene_expression_matrix <- GetAssayData(nsclc.seurat.obj, assay = "RNA")
# Extract gene names from the row names of the matrix
gene_names <- rownames(gene_expression_matrix)
head(gene_names)
gene_names_df <- data.frame(GeneName = gene_names)
View(gene_names_df)

# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')




# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)


# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)


# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top25 <- head(VariableFeatures(nsclc.seurat.obj), 25)
top25


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top25, repel = TRUE) 


# 5. Scaling -------------  TO remove th batch effects
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)



# 6. Perform Linear dimensionality reduction -------------- To identify the hetrogenicity 
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)




# 7. Clustering ------------ to find the most variation in the dataset 
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)


# understanding resolution  --> lower number lower cluster, higher number higher cluster 
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

# ploting the each cluster 
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)


# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap")


