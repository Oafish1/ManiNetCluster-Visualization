# NOTE: Crashes every other launch because of retriculate
# Libraries
library(class)
library(cluster)
library(ManiNetCluster)
library(plot3D)
library(plyr)
library(RColorBrewer)
library(shiny)
library(shinyjs)

# CTW returns only two dim for some reason?
alignments = c(
  "cca", 
  "linear manifold", 
  "manifold warping", 
  "nonlinear manifold aln",
  "nonlinear manifold warp"
)
names(alignments) = c(
  "CCA", 
  "Linear Manifold Alignment", 
  "Manifold Warping", 
  "Non-Linear Manifold Alignment",
  "Non-Linear Manifold Warping"
)

# UI
ui <- fluidPage(
  # Meta
  title="ManiNetCluster WebApp",
  
  # Main plot
  # asddf: Add interactive plotting
  useShinyjs(),
  h2("ManiNetCluster Webapp"),
  textOutput("warnings"),
  fluidRow(
    column(6, plotOutput("content1")),
    column(6, plotOutput("content2")),
  ),
  p(strong("Label Transfer Accuracy:"), textOutput("label_transfer_accuracy", inline=T)),
  hr(),
  
  # Upload/Download UI
  fluidRow(
    column(3,
           h2("Alignment"),
           selectInput("method", "Method", choices=names(alignments), selected="Non-Linear Manifold Alignment"),
           sliderInput("d", "Dimensions", min = 1, max = 20, value = 3),
           sliderInput("knn", "Nearest Neighbors", min = 0, max = 20, value = 3), 
           sliderInput("kmed", "Medoids", min = 0, max = 20, value = 6),
           strong("Note:"),
           p("To use medoids, nearest neighbors must be 0.")
    ),
    column(3,
           h2("Visualization"),
           selectInput("cluster_method", "Method", choices = c("None", "K-Means", "PAM")),
           sliderInput("num_clusters", "Clusters", min = 1, max = 7, value = 4),
    ),
    column(6,
      column(6,
        h2("Upload"),
        strong("Row and column labels are expected in all files.  All files are
                also assumed to be ordered similarly (by cell)."),
        br(), div(style = "margin-top: 30px"),
        
        # asddf: Replace div hack
        fileInput("mat1", NULL, buttonLabel="First modality", multiple=FALSE),
        div(style = "margin-top: -30px"),
        fileInput("mat2", NULL, buttonLabel="Second modality", multiple=FALSE),
        div(style = "margin-top: -30px"),
        p("Modality data should be of dimension (cells x features) in CSV format."),
        
        br(),
        fileInput("meta1", NULL, buttonLabel="First modality meta", multiple=FALSE),
        div(style = "margin-top: -30px"),
        fileInput("meta2", NULL, buttonLabel="Second modality meta", multiple=FALSE),
        div(style = "margin-top: -30px"),
        p("Meta data should be ordered similarly to modality data in CSV format.  A",
           em("time"), "column is expected.  Optionally, a", em("label"), "column",
          "may be included for additional error reporting."),
        
        br(),
        fileInput("corr", NULL, buttonLabel="Correspondence", multiple=FALSE),
        div(style = "margin-top: -30px"),
        p("Correspondence matrix (cells1 x cells2) should be of CSV format. ",
          "Meant to represent inter-dataset correspondence and can be calculated",
          "in a variety of ways.  If aligned and no CSV is provided, will default",
          "to the identity matrix."),
      ),
      column(6,
        h2("Download"),
        downloadButton("download", "Download"),
        br(), div(style = "margin-top: 30px"),
      )
    ),
  ),
  
)

# Server
server <- function(input, output, session) {
  output$warnings <- renderText({
    mat1 <- input$mat1
    mat2 <- input$mat2
    corr <- input$corr
    meta1 <- input$meta1
    meta2 <- input$meta2
    
    # Warnings in priority order
    if (is.null(mat1) || is.null(mat2))
      return("Default dataset is being used, please upload data")
    else {
      if (!is.null(corr)) {
        corr <- as.matrix(read.csv(corr$datapath, row.names=1))
      }
      else if (dim(mat1)[1] == dim(mat2)[1]) {
        validate(need(F, "No correspondence matrix provided, identity matrix being used"))
        corr = diag(dim(mat1)[1])
      }
      else {
        validate(need(F, "Could not construct correspondence matrix, please upload a csv"))
        return(NULL)
      }
    }
    if (is.null(meta1) || is.null(meta2))
      return("Default meta is being used, please upload meta")
  })
  
  get_method <- function(method_code) {
    for (i in 1:length(alignments)) {
      if (method_code == names(alignments)[i])
        return(alignments[i])
    }
    return(NULL)
  }
  
  # Change sliders based on selection
  refresh_dimensions <- reactive({
    # -----------------------------------------
    # Noted for adding new methods
    method = get_method(input$method)
    to_show = switch(
      get_method(input$method),
      "cca" = c(),
      "linear manifold" = c("knn"),
      "manifold warping" = c("knn"),
      "nonlinear manifold aln" = c("knn"),
      "nonlinear manifold warp" = c("knn")
    )
    
    # Make all invisible
    for (id in c("knn", "kmed"))
      shinyjs::disable(id=id)
    # Show chosen sliders
    shinyjs::enable(id="d")
    dimensions = input$d
    for (id in to_show)
      shinyjs::enable(id=id)
    if (input$knn == 0 && "knn" %in% to_show) {
      shinyjs::enable(id="kmed")
      shinyjs::disable(id="d")
      dimensions = input$kmed
    }
    return(dimensions)
  })
  observeEvent(input$method, refresh_dimensions())
  
  # Conditional sliders
  observeEvent(input$knn, {
    if (input$knn == 0) {
      shinyjs::enable(id="kmed")
      shinyjs::disable(id="d")
    }
    else {
      shinyjs::disable(id="kmed")
      shinyjs::enable(id="d")
    }
  })
  
  # Generate clusters / Get colors
  get_clusters <- function(df, default_color) {
    working_data = df[,3:(2+input$d)]
    # Max 7 colors
    colors = c("blue","green","yellow","orange","red","pink","purple")
    
    # asdf: Do clusters need to be same?  Even with non-aligned cells?
    if (input$cluster_method == "None")
      return (colorRampPalette(brewer.pal(n=9, default_color))(12))
    else if (input$cluster_method == "K-Means")
      clusters = kmeans(working_data, input$num_clusters)$cluster
    else if (input$cluster_method == "PAM")
      clusters = pam(working_data, input$num_clusters)$clustering
    else
      validate(need(F, "Invalid clustering method"))
    return (colors[clusters])
  }
    
  perform_alignment <- reactive({
    mat1 <- input$mat1
    mat2 <- input$mat2
    corr <- input$corr
    
    if (is.null(mat1) || is.null(mat2)) {
      mat1 = as.matrix(read.csv("data/mat1.csv", row.names=1))
      mat2 = as.matrix(read.csv("data/mat2.csv", row.names=1))
      corr = as.matrix(read.csv("data/corr.csv", row.names=1))
    }
    else {
      mat1 <- as.matrix(read.csv(mat1$datapath, row.names=1))
      mat2 <- as.matrix(read.csv(mat2$datapath, row.names=1))
      
      # Use KNN as correspondence if applicable
      if (!is.null(corr)) {
        corr <- as.matrix(read.csv(corr$datapath, row.names=1))
      }
      else if (dim(mat1)[1] == dim(mat2)[1]) {
        corr = diag(dim(mat1)[1])
      }
      else {
        return(NULL)
      }
    }
    XY_corr = Correspondence(matrix=corr)
    
    # Get vars
    method = get_method(input$method)
    
    # Run NLMA
    ManiNetCluster(
      mat1,mat2,
      nameX='sample1',nameY='sample2',
      corr=XY_corr,
      d=as.integer(input$d),
      method=method,
      k_NN=as.integer(input$knn),
      k_medoids=as.integer(input$kmed)
    )
  })
  
  plot_alignment = function(sample_label, default_color) {
    dimensions = refresh_dimensions()
    validate(need(dimensions >= 3, paste(
      "Must have dimensions \u22653 to plot, currently have", dimensions)))
    
    # Get inputs
    meta1 <- input$meta1
    meta2 <- input$meta2
    
    if (is.null(meta1) || is.null(meta2)) {
      meta1 = read.csv("data/meta1.csv", row.names=1)
      meta2 = read.csv("data/meta2.csv", row.names=1)
    }
    else {
      meta1 = read.csv(meta1$datapath, row.names=1)
      meta2 = read.csv(meta2$datapath, row.names=1)
    }
    
    # Plot code from paper
    df2 <- perform_alignment()
    # Assumes ordering from NLMA
    df2$time = c(meta1$time,meta2$time)
    res = data.frame(df2[df2$data==sample_label,])
    res0 = data.frame(df2)
    time.cols1 = get_clusters(res0, default_color)
    s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],
                   colvar=as.numeric(mapvalues(res$time,names(table(res$time)),c(2:13))),
                   col = c(time.cols1),pch=c(16,17)[as.numeric(as.factor(res$data))],
                   colkey=F, theta = 300, phi = 30,cex=2,
                   xlim=c(min(res0$Val0),max(res0$Val0)),
                   ylim=c(min(res0$Val1),max(res0$Val1)),
                   zlim=c(min(res0$Val2),max(res0$Val2)))
    
    #legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
    #legend("right", legend = levels(as.factor(res$time)), col = c(time.cols1),pch=16,inset = 0.1, xpd = TRUE, horiz = F,cex=2)
  }
  
  # Accuracy metrics
  output$label_transfer_accuracy <- renderText({
    # Get inputs
    meta1 <- input$meta1
    meta2 <- input$meta2
    
    if (is.null(meta1) || is.null(meta2)) {
      meta1 = read.csv("data/meta1.csv", row.names=1)
      meta2 = read.csv("data/meta2.csv", row.names=1)
    }
    else {
      meta1 = read.csv(meta1$datapath, row.names=1)
      meta2 = read.csv(meta2$datapath, row.names=1)
    }
    
    # Check for possibility
    if (!("label" %in% colnames(meta1) && "label" %in% colnames(meta2)))
      return ("Labels not provided")
    
    # Get calculations
    aligned_data = perform_alignment()
    source_data = data.frame(aligned_data[aligned_data$data=='sample1',])[,3:(2+input$d)]
    transfer_data = data.frame(aligned_data[aligned_data$data=='sample2',])[,3:(2+input$d)]
    source_labels = meta1$label
    transfer_labels = meta2$label
    
    predictions = knn(source_data, transfer_data, cl=source_labels, k=5)
    correct = sum(predictions == transfer_labels)
    print(transfer_labels)
    
    return (correct / length(transfer_labels))
  })
  
  # Outputs
  output$content1 <- renderPlot({plot_alignment('sample1', 'Greens')})
  output$content2 <- renderPlot({plot_alignment('sample2', 'Oranges')})
  
  output$download <- downloadHandler(
    "aligned.csv",
    function(fname) {
      write.csv(perform_alignment(), fname)
    },
    "text/csv"
  )
}

# Run
shinyApp(ui=ui, server=server)
