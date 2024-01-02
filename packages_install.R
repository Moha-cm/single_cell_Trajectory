
install.packages("tidyverse")
install.packages("ggthemes")
install.packages("DT")



install.packages("dplyr")
install.packages("shinyjs")

install.packages("shinydashboard")
install.packages("shinydashboardPlus")
install.packages("ggplot2")

install.packages("shinybusy")
install.packages("markdown")





if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}


# package installation  for monocle3

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                     'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                     'SummarizedExperiment', 'batchelor', 'HDF5Array',
                     'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')


 # seurat pacakage
install.packages('Seurat')

remotes::install_github('satijalab/seurat-wrappers')


# install searut disk


remotes::install_github("mojaveazure/seurat-disk")

# -----------------------------------------------------------------------------------------------------------------------


install.packages("umap")
 install.packages("shiny")

install.packages("plotly")


library(flexmix)
library(velocyto.R)


library(plotly)
library(shiny)
library(monocle3)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)
library(umap)
 











