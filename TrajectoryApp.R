# shiny::runApp()
source('global.R')





ui <- dashboardPage(
    dashboardHeader(title = "Option Menu "),
    dashboardSidebar(
     tags$head(
        tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
        ),
        # setting  the options Menu
        sidebarMenu( id = 'tab',
        
        #setting the menu item 
        menuItem("Home",tabName = "home",icon = icon("list")),
        useShinyjs(),
                                                               
        
        # setting the condition to the each options  
            menuItem("Trajectroy Analyzer",tabName = "input",icon = icon("edit")),
                conditionalPanel(condition = "input.tab == 'input'",
                    div(
                        fileInput("file","Upload File",multiple = FALSE,accept = c('.rds',".rds.gz",".txt")),
                        actionButton("reset","Reset",icon = icon("undo"),style="color: #fff; background-color: #dc3545; width:87.25%"),
                        actionButton("run","Run",icon = icon("play"),style="color: #fff; background-color: #28a745; width:87.25%")
                 )
            )
        )
    ),
    dashboardBody(
        # set the tabs references from the input tabs
        tabItems(
            #subset of tabs 
            tabItem(tabName = "input",tabsetPanel(
                id = "main_tabs",tabPanel("Instructions",includeMarkdown('markdown\\instructions.md'))
            )
            ), 
            tabItem(tabName = "home",tags$h1(HTML("<b>Welcome Singel Cell RNAseq Trajectory Analysis </b>")))
        )
    )
)


server <- function(input,output,session){
    
    options(shiny.maxRequestSize=3000*1024^2) # maximizing the input size 

   # Disable Run by default
    shinyjs::disable("run")

   observe({
    if(is.null(input$file) != TRUE) {
        shinyjs::enable("run")
    } else {
        shinyjs::disable("run")
    }
    })

    observeEvent(input$reset,{
        shinyjs::reset("file")
        shinyjs::disable("run")
        removeTab("main_tabs","UMAP")
        removeTab("main_tabs","Gene Expression")
    })

    observeEvent(input$run,{
    shinyjs::disable("run")
    show_modal_spinner(text = "Plotting the graphs ")

    obj <- load_seurat_obj(input$file$datapath) # contain the uploaded file path 


# Error show the pop up 
    if(is.vector(obj)){
        showModal(modalDialog(title = "Error",HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
        paste(unlist(obj), collapse = "<br><br>"))))
        shinyjs::enable("run")}
    else {
        
        output$umap <-renderPlotly({create_metadata_UMAP(obj,input$metadata_col)})    # getting the umaps plots 
        output$dim_plot_v <- renderPlotly({generate_dim_plot(obj,input$metadata_col)}) # dimplots


        output$feature_plot <- renderPlotly({create_feature_plot(obj,input$gene1)})    # getting the feature plots



        output$umap <- renderPlotly({
              create_metadata_UMAP(obj, input$metadata_col_umap)
              })

        output$dim_plot_v <- renderPlotly({
              generate_dim_plot(obj, input$metadata_col_dim_plot)
              })
        
       


       # output$violin_plot2 <- renderPlotly({create_violen_plot2(obj,input$gene2,input$metadata_col)})     # getting the violen plots 
        

        output$violin_plot <- renderPlotly({
            gene <- input$gene2
            col <- input$metadata_col1
            if (gene %in% rownames(obj)) {
                         plots <- list()# Create a list to store individual violin plots
            for (g in gene) {   # Loop through each gene and create a violin plot
                plot_gene <- create_violen_plot(obj, g, col)
                plots[[g]] <- plot_gene
            }
        
        # Combine individual plots into a subplot
        subplot(plots)
    } else {
        # Display an empty plot if the gene is not found
        plot_ly() %>%
            layout(showlegend = FALSE) %>%
            add_trace(type = "violin")
    }
        })
        

        output$v_gene_plot <- renderPlotly({    gene <- input$gene6
            if (length(gene) > 0) {
                      violen_plot_gene(obj, gene)
                          }
                        })
                    
        cds <- learning_trajectories(obj)
        print(cds)


        output$indent_traj <- renderPlotly({create_trajectory_plot_indent(cds)})

        output$indent_clusters <- renderPlotly({create_trajectory_plot_clusters(cds)})

  
        pseudo_time <- pseudotime_analysis(cds)
        print( pseudo_time)

        
        output$pseudo_plot_c <- renderPlotly({create_trajectory_plot_clusters(pseudo_time)})
        output$pseudo_plot_in <- renderPlotly({create_trajectory_plot_indent(pseudo_time)})

        pseudo_time_df <- pseudotime_df(pseudo_time)
        

       
        
        size_factors <- estimate_size(pseudo_time)

        print(size_factors)
        differential_genes_df <- df_genes(size_factors)
        print( differential_genes_df )


        output$D_feature_plot <- renderPlot({create_feature_plot(differential_genes_df,input$gene5)}) 
        
        print( head(pseudo_time_df))
        print(colnames(pseudo_time_df))

        output$boxplot <-renderPlot({create_boxplot(pseudo_time_df,colnames(pseudo_time_df),colnames(pseudo_time_df),colnames(pseudo_time_df)
        )}) 
        


       

        insertTab(
            inputId = "main_tabs",tabPanel("UMAP Plots",
        fluidRow(
        column(width = 8, plotlyOutput(outputId ="umap"), downloadButton('download_umap', "Download UMAP")),
        column(width = 4, selectizeInput("metadata_col_umap", "Metadata Column", colnames(obj@meta.data))),
        column(width = 8, plotlyOutput(outputId ="dim_plot_v"), downloadButton('download_dim_plot', "Download Dim Plot")),
        column(width = 4, selectizeInput("metadata_col_dim_plot", "Metadata Column", colnames(obj@meta.data)))
         )
         ))


        insertTab(inputId = "main_tabs",
            tabPanel("Gene Expression",
            fluidRow(
                column(width = 8,plotlyOutput(outputId ="feature_plot"),downloadButton('downloadFeaturePlot',"Download Feature Plot")),
                column(width=4,selectizeInput ( "gene1","Genes",rownames(obj))),
                column(width = 8, plotlyOutput(outputId ="v_gene_plot"), downloadButton('downloadVGenePlot', "Download Violin Gene Plot")),
        column(width = 4, selectizeInput("gene6", "Genes", rownames(obj)))
                
              )
            ))


        insertTab(inputId = "main_tabs",
        tabPanel("Trajectory Plots",
        fluidRow(
      column(width = 6,
        plotlyOutput(outputId ="indent_traj"),
        downloadButton('download_traj_plot', "Download Trajectory Plot")
      ),
      column(width = 6,
        plotlyOutput(outputId ="indent_clusters"),
        downloadButton('download_clusters_plot', "Download Clusters Plot")
      )
    )
  )
)
        insertTab(inputId = "main_tabs",
        tabPanel("Pseudotime Plots",
        fluidRow(
      column(width = 6,
        plotlyOutput(outputId ="pseudo_plot_in"),
        downloadButton('download_traj_plot', "Download Trajectory Plot")
      ),
      column(width = 6,
        plotlyOutput(outputId ="pseudo_plot_c"),
        downloadButton('download_clusters_plot', "Download Clusters Plot")
      )
    )
  )
)
           
         
         
            remove_modal_spinner()
            shinyjs::disable(("run"))
            }
    })

     # Clear all sidebar inputs when 'Reset' button is clicked
    observeEvent(input$reset, {
        shinyjs::reset("file")
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")
        removeTab("main_tabs","Pseudotime Plots")
        removeTab("main_tabs","Trajectory Plots")
        shinyjs::disable("run")
    })
}


shinyApp(ui,server)
