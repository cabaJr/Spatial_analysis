pkg_req <- c("shiny", "shinyFiles", "plotly", "fs", "tools", "htmltools")


for (package in pkg_req) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

ui <- fluidPage(
  titlePanel("3D Plot Comparison Dashboard"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("folder1", "Select Folder 1", "Choose first folder"),
      shinyDirButton("folder2", "Select Folder 2", "Choose second folder"),
      checkboxInput("sync_rotation", "Sync 3D Rotation", TRUE),
      width = 3
    ),
    mainPanel(
      uiOutput("plot_ui"),
      width = 9
    )
  )
)

server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), getVolumes()())
  
  shinyDirChoose(input, "folder1", roots = volumes, session = session)
  shinyDirChoose(input, "folder2", roots = volumes, session = session)
  
  folder1_path <- reactive({
    req(input$folder1)
    return(parseDirPath(volumes, input$folder1))
  })
  
  folder2_path <- reactive({
    req(input$folder2)
    return(parseDirPath(volumes, input$folder2))
  })
  
  matched_files <- reactive({
    req(folder1_path(), folder2_path())
    files1 <- dir_ls(folder1_path(), glob = "*.html")
    files2 <- dir_ls(folder2_path(), glob = "*.html")
    
    common_files <- intersect(basename(files1), basename(files2))
    data.frame(name = common_files,
               path1 = file.path(folder1_path(), common_files),
               path2 = file.path(folder2_path(), common_files))
  })
  
  observe({
    req(folder1_path(), folder2_path())
    
    # Ensure both folders are accessible
    addResourcePath("folder1", folder1_path())
    addResourcePath("folder2", folder2_path())
  })
  
  
  output$plot_ui <- renderUI({
    req(matched_files())
    plots <- lapply(seq_len(nrow(matched_files())), function(i) {
      fluidRow(
        column(6, uiOutput(paste0("plot1_", i))),
        column(6, uiOutput(paste0("plot2_", i)))
      )
    })
    do.call(tagList, plots)
  })
  
  observe({
    req(matched_files())
    
    lapply(seq_len(nrow(matched_files())), function(i) {
      output[[paste0("plot1_", i)]] <- renderUI({
        # Check if the HTML file exists and add it to iframe
        tags$iframe(src = file.path("folder1", basename(matched_files()$path1[i])), 
                    width = "100%", height = "600px", frameborder = 0)
      })
      
      output[[paste0("plot2_", i)]] <- renderUI({
        # Check if the HTML file exists and add it to iframe
        tags$iframe(src = file.path("folder2", basename(matched_files()$path2[i])), 
                    width = "100%", height = "600px", frameborder = 0)
      })
    })
  })
}

shinyApp(ui, server)
