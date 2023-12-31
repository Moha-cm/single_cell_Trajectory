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
        # setting the menuitem
        sidebarMenu(

            id =  "sidebar", menuItem("Home",tabName = "home",icon = icon("list"))),
        
            #p_seruat_file --> file uploading to Preprocessing Trajectory
            # tar_file -->  file uploding to Analyze  Trajectory
            # hds_file --> file uploading to Preprocessing Seurat 
            # rds_file --> file  uploding to Analyze  Seurat 



            # first menu item
            menuItem(text ="Preprocessing Trajectory File",tabName = "file_t",icon = icon("edit"),
            conditionalPanel(condition = "file_t.sidebar == 'file_t'",
             div(
                fileInput("p_seruat_file","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:87%"),
                actionButton("run","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:87%"))
            )),
            
            #second menu item
            menuItem(text ="Analyze  Trajectory ",tabName = "An_traj",icon = icon("edit"),
            conditionalPanel(condition = "An_traj.sidebar == 'An_traj'",
             div(
                fileInput("tar_file","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:87%"),
                actionButton("run","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:87%")))),

            # 3rd menu item
            menuItem(text ="Preprocess seurat",tabName = "file_s",icon = icon("edit"),
            conditionalPanel(condition = "file_s.sidebar == 'file_s'",
             div(
                fileInput("hds_file","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:87%"),
                actionButton("run","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:87%")))),
            
            # 4th menu item
            menuItem(text ="Analyze seurat ",tabName = "An_s",icon = icon("edit"),
            conditionalPanel(condition = "An_s.sidebar == 'An_s'",
             div(
                fileInput("rds_file","Upload File",multiple = FALSE,accept = c('.rds')),
                actionButton("reset","Reset",icon = icon("undo"),style ="color: #fff; background-color: #dc3545; width:87%"),
                actionButton("run","Run",icon = icon("play"),style ="color: #fff; background-color: #28a745; width:87%"))))
            
    ),
    dashboardBody()
)



# setting the server which passes the input and output into the ui interference 

server <- function(input,output,session){  
}


# compile the UI and server 

shinyApp(ui,server)