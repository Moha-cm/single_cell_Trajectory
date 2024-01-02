library(shinyjs)
library(shiny)
library(shinybusy)
library(shinydashboard)
library(shinydashboardPlus)

# setting the UI interference 
ui <- dashboardPage(
    dashboardHeader(title = "Option Menu"),
    # sidebar menu 
    dashboardSidebar(

        tags$style(HTML(".sidebar-menu { margin-bottom: 20px; /* Adjust the margin as needed */ } ")),
        tags$style(HTML(".sidebar-menu .treeview-menu > li {margin-bottom: 10px; /* Adjust the margin as needed */}")),
       # tags$style(HTML(".main-sidebar {width: 250px; /* Adjust the sidebar width as needed */}")),
        
        # setting the menuitem

        sidebarMenu(
            useShinyjs(), 
            id =  "sidebar", menuItem("Home",tabName = "home",icon = icon("list")),
        
            # first menu item
            
            menuItem(text ="Preprocessing Trajectory ",tabName = "file_t",icon = icon("edit"),
            conditionalPanel(condition = "file_t.sidebar == 'file_t'",
             div(
                fileInput("file1","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset1","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:45%"),
                actionButton("run1","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:45%"))
            )),
            
            #second menu item
            menuItem(text ="Analyze  Trajectory ",tabName = "An_traj",icon = icon("edit"),
            conditionalPanel(condition = "An_traj.sidebar == 'An_traj'",
             div(
                fileInput("file2","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset2","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:45%"),
                actionButton("run2","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:45%")))),

            # 3rd menu item
            menuItem(text ="Preprocess seurat",tabName = "file_s",icon = icon("edit"),
            conditionalPanel(condition = "file_s.sidebar == 'file_s'",
             div(
                fileInput("file3","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset3","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:45%"),
                actionButton("run3","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:45%")))),
            
            # 4th menu item
            menuItem(text ="Analyze seurat ",tabName = "An_s",icon = icon("edit"),
            conditionalPanel(condition = "An_s.sidebar == 'An_s'",
             div(
                fileInput("file4","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset4","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:45%"),
                actionButton("run4","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:45%"))))
    )
    ),

    dashboardBody(

        # setting the Home page  Tabs 

        tabItems(
            tabItem(tabName = "home",
            tabsetPanel(id = "main_An_tab",
            tabPanel("Basic Instructions",includeMarkdown("markdown\\basic_instructions.md")),
            tabPanel("Trajectory Data PreProcessing",includeMarkdown("markdown\\preprocessing_trajectory.md")),
            tabPanel("Analyze Trajectory  Data",includeMarkdown("markdown\\trajetory_analysis_markdown.md")),
            tabPanel("seurat Data  PreProcessing",includeMarkdown("markdown\\preprocessing_seurat_markdown.md")),
            tabPanel("Analyze seurat Data",includeMarkdown("markdown\\Analyze_seurat_markdown.md")))
            ),
            tabItem(tabName = "home")
        )
            )     
)


# setting the server which passes the input and output into the ui interference 

server <- function(input,output,session){  
    
    # setting the maximum size 
     options(shiny.maxRequestSize=30000*1024^2)


    # setting the Run button disable 
    shinyjs::disable('run1')
    shinyjs::disable('run2')
    shinyjs::disable('run3')
    shinyjs::disable('run4')
    observe({
    if(is.null(input$file1) != TRUE) {
        shinyjs::enable("run1")
    } else {
        shinyjs::disable("run1")
    }
    })

    observe({
    if(is.null(input$file2) != TRUE) {
        shinyjs::enable("run2")
    } else {
        shinyjs::disable("run2")
    }
    })

    observe({
    if(is.null(input$file3) != TRUE) {
        shinyjs::enable("run3")
    } else {
        shinyjs::disable("run3")
    }
    })

    observe({
    if(is.null(input$file4) != TRUE) {
        shinyjs::enable("run4")
    } else {
        shinyjs::disable("run4")
    }
    })

    # setting the reset button click 

    observeEvent(input$reset1,{
        shinyjs::reset("file1")
        shinyjs::disable("run1")
        #removeTab("main_tabs","UMAP")
        #removeTab("main_tabs","Gene Expression")
    })

    observeEvent(input$reset2,{
        shinyjs::reset("file2")
        shinyjs::disable("run2")
        #removeTab("main_tabs","UMAP")
        #removeTab("main_tabs","Gene Expression")
    })

    observeEvent(input$reset3,{
        shinyjs::reset("file3")
        shinyjs::disable("run3")
       # removeTab("main_tabs","UMAP")
        #removeTab("main_tabs","Gene Expression")
    })

    observeEvent(input$reset4,{
        shinyjs::reset("file4")
        shinyjs::disable("run4")
       # removeTab("main_tabs","UMAP")
       # removeTab("main_tabs","Gene Expression")
    })

}


# compile the UI and server 

shinyApp(ui,server)




# ---------------------------------------- compeleted getting the input and disable the reset  button --------------------------------------------------------
