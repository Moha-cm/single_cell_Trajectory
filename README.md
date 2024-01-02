# Project Title : single_cell_Trajectory Analysis

Single_cell_Trajectory Analysis: A User-Friendly Tool Using R shiny for datananlysis tool for ScRNaseq gene expression data

## **Overview**
The primary goal of the tool is to facilitate trajectory analysis on scRNA-seq data. Trajectory analysis helps in understanding the developmental paths of individual cells over time, providing insights into cell differentiation, transitions, and heterogeneity.


## Approach

1. **Technology Stack:**
       **R Shiny:** The tool is built using the R Shiny framework, which allows for the creation of interactive web applications directly from R.

2. **Data Input:**
   
   - The tool likely supports the input of scRNA-seq gene expression data, which is a matrix representing the expression levels of genes across individual cells.
Users may be able to upload their seurat data or use create seurat dataset with the help of  markdown code provided.

3. **Data Analysis**
   - Apart form the dashboard, the trajectory workflow markdown code  is provide for the your analysis. Using the dashboard help  you can to visualize the UMAP,Gene expression of each  gene, Trajectory  and Pseudotime Plots  more easily.

## R packages
use **packages_install.R** file ,to download the packages required for the dashboard 


Dowload the source files from repo  and  use the bellow commandas to run

## Run the application using the following command
```
shiny::runApp('Trajectory_app.R')
```
This will launch the Rshipy application.
