library(tools)
library(Seurat)
library(dplyr)
library(shinyjs)
library(shiny)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(markdown)
library(ggthemes)
library(plotly)
library(monocle3)
library(dplyr)
library(tidyverse)
library(SeuratWrappers)
library(SingleCellExperiment)




# Read in file and perform validation.
load_seurat_obj <- function(path){
    errors <- c()
    # check file extension
    if (!tolower(tools::file_ext(path)) == "rds") { # ignores case
        errors <- c(errors, "Invalid rds file.")
        return(errors)
    }

    # try to read in file
    tryCatch(
        {
        obj <- readRDS(path)
        },
        error = function(e) {
            errors <- c(errors, "Invalid rds file.")
            return(errors)
        }
    )

    # Validate obj is a seurat object
    if (!inherits(obj, "Seurat")){
        errors <- c(errors, "File is not a seurat object")
        return(errors)
    }

    return(obj)
}


#2d viloen plot 
violen_plot_gene  <- function(obj, gene) {
    if (gene %in% rownames(obj)) {
        vp <- VlnPlot(obj, features = gene)
    } else {
        vp <- ggplot() + 
            theme_void() + 
            VlnPlot(obj, features = gene) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    return(vp)
}


# dimplot
generate_dim_plot <- function(obj, group.by) {
    dim_plot <- DimPlot(obj, reduction = "umap", group.by = group.by,label =TRUE)
    return(dim_plot)
}



#2d viloen plot 
rp_plot_gene  <- function(obj, gene) {
    if (gene %in% rownames(obj)) {
        rp <- RidgePlot(obj, features = gene)
    } else {
        rp <- ggplot() + 
            theme_void() + 
            RidgePlot(obj, features = gene) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    return(rp)
}


#function without trajectory 

create_plot_indent <- function(obj) {
   
   cluster_plot <- plot_cells(obj, color_cells_by = "ident", label_groups_by_cluster = FALSE, group_label_size = 5) +
      theme(legend.position = "right")
      
   return(cluster_plot)
}


create_plot_clusters <- function(obj) {
   
   cluster_plot <- plot_cells(obj, color_cells_by = "seurat_clusters", label_groups_by_cluster = FALSE, group_label_size = 5) +
      theme(legend.position = "right")
      
   return(cluster_plot)
}




# function for trajectory learning 

learning_trajectories <- function(obj) {

    
    obj@meta.data[["ident"]] <- obj@active.ident
    
    cds <- as.cell_data_set(obj, group.by = "ident")
    
    DefaultAssay(obj) <- "RNA"
    
    fData(cds)$gene_short_name <- rownames(fData(cds))
    
    recreate.partion  <- c(rep(1, length(cds@colData@rownames)))
    names(recreate.partion)  <- cds@colData@rownames
    recreate.partion <- as.factor(recreate.partion)
    
    cds@clusters$UMAP$partitions <- recreate.partion
    
    list_cluster <- obj@active.ident
    cds@clusters$UMAP$clusters <- list_cluster
    
    # Assign UMAP coordinates - cell embeddings 
    cds@int_colData@listData$reducedDims$UMAP  <- obj@reductions$umap@cell.embeddings
    
    cds <- learn_graph(cds, use_partition = FALSE)
    
    return(cds)
    }

create_trajectory_plot_clusters <- function(obj) {
   
   cluster_plot <- plot_cells(obj, color_cells_by = "seurat_clusters", label_groups_by_cluster = FALSE, group_label_size = 5) +
      theme(legend.position = "right")
      
   return(cluster_plot)
}



create_trajectory_plot_indent <- function(obj) {
   
   cluster_plot <- plot_cells(obj, color_cells_by = "ident", label_groups_by_cluster = FALSE, group_label_size = 5) +
      theme(legend.position = "right")
      
   return(cluster_plot)
}


# for Pseudotime analysis 
pseudotime_analysis <- function(cds) {
  # Order the cells based on the pseudotime analysis
  root_cells <- colnames(cds[, clusters(cds) == "Adventitia"])
  cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = root_cells)
  return(cds)
}



# getting the pseudataframe for box plot
pseudotime_df <-function(pseudo_obj){
  pseudo_obj$molocle3_pseudotime <- pseudotime(pseudo_obj)
  data.pseudo <- as.data.frame(colData(pseudo_obj))
  return(data.pseudo)
}


# estimate size factors 
estimate_size <- function(cds){
  cds <- estimate_size_factors(cds)
  cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(cds[["RNA"]])
  return(cds)
}


# getting the diffrential genes

df_genes <- function(cds){
    
    Differenial_expressed_genes<- graph_test(cds,neighbor_graph = "principal_graph",cores = 3)
    
    
    differential_genes_filtered <-Differenial_expressed_genes %>%
      arrange(q_value) %>%
      filter(status=="OK")
}



#creating the UMAP 
create_metadata_UMAP <- function(obj, col){
    if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")){
    col_df <- data.frame(obj@reductions$umap@cell.embeddings, data = obj@meta.data[,col])
    umap <- ggplot(data = col_df) +
        geom_point(mapping = aes(umap_1, umap_2, color = log10(data)), size = 0.01) +
        scale_colour_gradientn(colours = rainbow(7))
    } else if (col %in% colnames(obj@meta.data)) {
    umap <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "umap")
    } else {
    umap <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    class_of_umap <- class(umap)
    print(class_of_umap) 
    return(umap)
}



#create feature plot 
create_feature_plot <- function(obj, gene) {
    if (gene %in% rownames(obj)) {
        FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
    } else {
        FP <- ggplot() + 
        theme_void() + 
        geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    return(FP)
}






#create violin plot for genes
create_violen_plot <- function(obj, gene, col) {
    if (gene %in% rownames(obj)) {
        df <- obj@meta.data
        vp <- plot_ly(data = df, y = ~df[[col]], x = obj@assays$RNA_Seq$data[, gene], type = "violin", box = list(visible = TRUE), line = list(color = 'black')) %>%
            layout(title = gene)
    } else {
        vp <- plot_ly() %>%
            layout(title = gene) %>%
            layout(showlegend = FALSE) %>%
            add_trace(type = "violin")
    }
    return(vp)
}













create_boxplot <- function(data, x_col, y_col,fill_v,a = median) {
  ggplot(data, aes_string(x = x_col, y = reorder(y_col, fill =fill_v,FUN = a ))) +
    geom_boxplot() +
    labs(title = paste("Boxplot of", y_col, "by", x_col),
         x = x_col,
         y = y_col) +
    theme_minimal()
}


