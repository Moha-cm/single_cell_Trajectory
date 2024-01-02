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
        
        menuItem("Home",tabName = "home",icon = icon("list")),
        useShinyjs(),                                                              # to dynamically disenable and enabling while uploading the file 
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

    if(is.vector(obj)){
        showModal(modalDialog(title = "Error",HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
        paste(unlist(obj), collapse = "<br><br>"))))
        shinyjs::enable("run")}
    else {
        
        output$umap <-renderPlot({create_metadata_UMAP(obj,input$metadata_col)})    # getting the umaps plots 
        output$feature_plot <- renderPlot({create_feature_plot(obj,input$gene1)})    # getting the feature plots
           # getting the feature plots

        
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
        
        cds <- learning_trajectories(obj)
        print(cds)
        pseudo_time <- pseudotime_analysis(cds,set_root = "Adventitia")
        print( pseudo_time)

        output$traj_plot  <- renderPlot({create_trajectory_plot(cds)})#
        output$pseudo_plot  <- renderPlot({create_trajectory_plot(pseudo_time)})

        pseudo_time_df <- pseudotime_df(pseudo_time)
        

       
        
        size_factors <- estimate_size(pseudo_time)

        print(size_factors)
        differential_genes_df <- df_genes(size_factors)
        print( differential_genes_df )
        print(row.names(differential_genes_df))
        print(colnames(differential_genes_df))


        output$D_feature_plot <- renderPlot({create_feature_plot(differential_genes_df,input$gene5)}) 
        
        print( head(pseudo_time_df))
        print(colnames(pseudo_time_df))

        output$boxplot <-renderPlot({create_boxplot(pseudo_time_df,colnames(pseudo_time_df),colnames(pseudo_time_df),colnames(pseudo_time_df)
        )}) 

        
        # setting the function for downloading the umap 
        output$download_umap <- downloadHandler(
                filename = function(){
                    paste0(input$metadata_col, '_UMAP', '.png')
                },
                content = function(file){
                    plot <- create_metadata_UMAP(obj, input$metadata_col)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )
    # setting the function for downloading the feature plot 
        output$downloadFeaturePlot <- downloadHandler(
                filename = function(){
                    paste0(input$gene, '_feature_plot', '.png')
                },
                content = function(file){
                    plot <- create_feature_plot(obj, input$gene)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )

        output$downloadviolenPlot <- downloadHandler(
                filename = function(){
                    paste0(input$gene, '_violen_plot', '.png')
                },
                content = function(file){
                    plot <- create_violin_plot(obj, input$gene)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )


      #  setting the slider meny for ploting 
        insertTab(inputId = "main_tabs",
            tabPanel("UMAP",
            fluidRow(
                column(width = 8,plotOutput(outputId ="umap"),downloadButton('download_umap',"Download UMAP")),
                column(width=4,selectizeInput ( "metadata_col","Metadata Column",colnames(obj@meta.data)))
                )
            )
            )

            insertTab(inputId = "main_tabs",
            tabPanel("Gene Expression",
            fluidRow(
                column(width = 8,plotOutput(outputId ="feature_plot"),downloadButton('downloadFeaturePlot',"Download Feature Plot")),
                column(width=4,selectizeInput ( "gene1","Genes",rownames(obj)))
                )  )
            )


            insertTab(inputId = "main_tabs",
            tabPanel("violen Plot",
            fluidRow(
                column(width = 8, plotlyOutput(outputId ="violin_plot"),downloadButton('downloadviolenPlot',"Download violen Plot")),
                column(width=2,selectizeInput ( "gene2","Genes",rownames(obj)), selectInput("metadata_col1", "Meta columns", choices = colnames(obj@meta.data)))
            )
            )
            )

            insertTab(inputId = "main_tabs",
            tabPanel("Trajectory Plot",
            fluidRow(
                column(width = 8, plotOutput(outputId ="traj_plot"),downloadButton('downloadtrajPlot',"Download Trajectory Plot"))
            )
            )
            )

            insertTab(inputId = "main_tabs",
            tabPanel("Pseudo Plot",
            fluidRow(
                column(width = 8, plotOutput(outputId ="pseudo_plot"),downloadButton('downloadPseudoPlot',"Download Pseudo Plot"))
            )
            )
            )

            insertTab(inputId = "main_tabs",
            tabPanel("Pseudo  Box plot ",
            fluidRow(
                column(width = 8, plotOutput(outputId ="boxplot"),downloadButton('downloadPseudoPlot',"Download Pseudo Plot")),
                column(width=2,selectInput ( "x1","X _value",colnames(pseudo_time_df)), selectInput("Y1", "Y_value", choices = colnames(pseudo_time_df),selectInput("Y2","Fill",choices = colnames(pseudo_time_df))))
            )
            )
            )

            
           
            insertTab(inputId = "main_tabs",
            tabPanel("Diffrential Gene Expression",
            fluidRow(
                column(width = 8,plotOutput(outputId ="D_feature_plot"),downloadButton('downloadFeaturePlot',"Download Feature Plot")),
                column(width=4,selectizeInput ( "gene5","Genes",rownames(differential_genes_df)))
                )  )
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
        shinyjs::disable("run")
    })
}


shinyApp(ui,server)
