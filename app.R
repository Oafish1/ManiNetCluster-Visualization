# NOTE: Crashes every other launch because of retriculate
# Libraries
library(class)
library(cluster)
library(ManiNetCluster)
library(pdist)
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
  title="ManiNetCluster Alignment Visualization Tool",
  
  # Main plot
  # asddf: Add interactive plotting
  useShinyjs(),
  h2("ManiNetCluster Webapp"),
  textOutput("warnings"),
  fluidRow(
    column(4, plotOutput("content1")),
    column(4, plotOutput("content2")),
    column(4, plotOutput("heatmap")),
  ),
  hr(),
  
  # UI
  fluidRow(
    column(5,
      h2("Upload"),
      strong("Row and column labels are expected in all files.  All files are
        also assumed to be ordered similarly (by sample)."),
      br(), div(style = "margin-top: 30px"),
      
      div(style="display:inline-block;margin-bottom:-30px",
        fluidRow(
        column(6,
          fileInput("mat1", NULL, buttonLabel="First Dataset")
        ),
        column(6,
          fileInput("mat2", NULL, buttonLabel="Second Dataset")
        ),
      )),
      fluidRow(
        column(8, p("Data should be of dimension (samples x features) in CSV format.")),
        column(4, img(src='data.png', align = "left", width=80)),
      ),
      
      br(),
      div(style="display:inline-block;margin-bottom:-30px",
          fluidRow(
          column(6,
            fileInput("meta1", NULL, buttonLabel="First Metadata"),
          ),
          column(6,
            fileInput("meta2", NULL, buttonLabel="Second Metadata"),
          ),
      )),
      fluidRow(
        column(8, p("Metadata should be ordered similarly to dataset data in CSV format.  A",
         em("time"), "column is expected.  All other columns are assumed to be categorical",
         "and are used for error reporting.")),
        column(4, img(src='metadata.png', align = "left", width=80)),
      ),
      
      br(),
      fileInput("corr", NULL, buttonLabel="Correspondence", multiple=FALSE),
      div(style = "margin-top: -30px"),
      fluidRow(
        column(8,
        p("Correspondence matrix (samples1 x samples2) should be of CSV format. ",
          "Meant to represent inter-dataset correspondence and can be calculated",
          "in a variety of ways.  If aligned and no CSV is provided, will default",
          "to the identity matrix.")),
      column(4, img(src='corr.png', align = "left", width=80)),
      ),
    ),
    column(4,
      h2("Alignment"),
      selectInput("method", "Method", choices=names(alignments), selected="Non-Linear Manifold Alignment"),
      sliderInput("d", "Dimensions", min = 1, max = 20, value = 3),
      sliderInput("knn", "Nearest Neighbors", min = 0, max = 20, value = 3), 
      sliderInput("kmed", "Medoids", min = 0, max = 20, value = 6),
      strong("Note:"),
      p("To use medoids, nearest neighbors must be 0."),
      hr(),
      downloadButton("download", "Download Aligned Data"),
      div(style = "margin-top: 15px"),
      p("Download data as a single CSV using the chosen parameters.  Includes clusters, if applicable."),
    ),
    column(3,
      h2("Clustering"),
      checkboxInput("show_clusters", "Show Clusters in Plot", value=F),
      selectInput("cluster_method", "Method", choices=c("K-Means", "PAM"), selected="PAM"),
      sliderInput("num_clusters", "Clusters", min = 1, max = 7, value = 4),
      
      # h2("Gene Information"),
      # em("Placeholder:"),
      # tableOutput("gci"),
      
      hr(),
      plotOutput("label_transfer_accuracy", height=200),
      em("Labels are taken from columns present in both metadata files.  zero values",
        "are assumed to be errant and are excluded.  Case-sensitive.  If no",
        "barplot is shown, no nonzero labels were detected.")
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
  get_clusters <- function(df, default_color="green", normal_colors=F) {
    working_data = df[,3:(2+input$d)]
    # Max 7 colors
    colors = c("blue","green","yellow","orange","red","pink","purple")
    
    if (normal_colors)
      return (list("colors"=colorRampPalette(brewer.pal(n=9, default_color))(12)))
    else if (input$cluster_method == "K-Means")
      clusters = kmeans(working_data, input$num_clusters)$cluster
    else if (input$cluster_method == "PAM")
      clusters = pam(working_data, input$num_clusters)$clustering
    else if (input$cluster_method == "Spectral") {
      s1 <- as.matrix(data.frame(df[df$data=="sample1",])[,3:(2+input$d)])
      s2 <- as.matrix(data.frame(df[df$data=="sample2",])[,3:(2+input$d)])
      dist = as.matrix(pdist(s1, s2))
    }
    else
      validate(need(F, "Invalid clustering method"))
    return (list("clusters"=clusters, "colors"=colors[clusters]))
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
  
  plot_alignment = function(sample_label, default_color, use_legend=F) {
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
    time.cols1 = get_clusters(res0, default_color, !input$show_clusters)$colors
    if (!input$show_clusters)
    {
      par(mar=c(5.1,4.1,4.1,4.1), xpd=T)
      s3d<-scatter3D(main=paste("Dataset", substr(sample_label, nchar(sample_label), nchar(sample_label))),
                     x=res[,3],y=res[,4],z=res[,5],
                     # colvar=as.numeric(mapvalues(res$time,names(table(res$time)),c(2:13))),
                     colvar=as.numeric(res$time),
                     col = c(time.cols1), pch=16, #pch=c(16,17)[as.numeric(as.factor(res$data))],
                     colkey=F, theta = 300, phi = 30,cex=2,
                     xlim=c(min(res0$Val0),max(res0$Val0)),
                     ylim=c(min(res0$Val1),max(res0$Val1)),
                     zlim=c(min(res0$Val2),max(res0$Val2)))
    }
    else {
      s3d<-scatter3D(main=paste("Dataset", substr(sample_label, nchar(sample_label), nchar(sample_label))),
                     x=res[,3],y=res[,4],z=res[,5],
                     col = c(time.cols1), pch=16, #pch=c(16,17)[as.numeric(as.factor(res$data))],
                     colkey=F, theta = 300, phi = 30,cex=2,
                     xlim=c(min(res0$Val0),max(res0$Val0)),
                     ylim=c(min(res0$Val1),max(res0$Val1)),
                     zlim=c(min(res0$Val2),max(res0$Val2)))
    }
    
    if (!input$show_clusters)
      legend("right", title="Time", legend=levels(as.factor(res$time)), col = c(time.cols1), pch=16, inset = c(0,0), horiz = F,cex=1)
    #legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
  }
  
  plot_heatmap <- function() {
    aligned <- perform_alignment()
    clusters <- get_clusters(aligned)
    
    aligned = aligned[sort(clusters$clusters, index.return=T)$ix,]
    clusters$colors = clusters$colors[sort(clusters$clusters, index.return=T)$ix]
    clusters$clusters = clusters$clusters[sort(clusters$clusters, index.return=T)$ix]
    clusters = clusters[sort(clusters$clusters, index.return=T)$ix]
    
    s1 <- as.matrix(data.frame(aligned[aligned$data=="sample1",])[,3:(2+input$d)])
    s2 <- as.matrix(data.frame(aligned[aligned$data=="sample2",])[,3:(2+input$d)])
    rsc = clusters$colors[aligned$data=="sample1"]
    csc = clusters$colors[aligned$data=="sample2"]
    
    s1 = s1[nrow(s1):1,]
    rsc = rsc[length(rsc):1]
    
    dist = as.matrix(pdist(s1, s2))
    
    par(mar=c(5.1,4.1,4.1,4.1), xpd=T)
    ramp = colorRampPalette(rev(brewer.pal(9, "Greys")))
    heatmap(dist, main="Distance Heatmap", RowSideColors=rsc, ColSideColors=csc, Colv=NA, Rowv=NA, labCol=F, labRow=F, col=ramp(256))
    legend("topright", legend=c("close", "mid", "far"), inset = c(0,0), fill=ramp(3))
    legend("bottomright", legend = paste("Cluster ", unique(clusters$clusters), sep=""), col = unique(clusters$colors), pch=16, inset = c(0,0), xpd = TRUE, horiz = F)
  }
  
  # Accuracy metrics
  output$label_transfer_accuracy <- renderPlot({
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
    
    # Prerequisite data
    aligned_data = perform_alignment()
    source_data = data.frame(aligned_data[aligned_data$data=='sample1',])[,3:(2+input$d)]
    transfer_data = data.frame(aligned_data[aligned_data$data=='sample2',])[,3:(2+input$d)]
    
    # Iterate through labels
    labels = c()
    acc = c()
    for (cname in colnames(meta1)) {
      if (!(cname %in% colnames(meta1) && cname %in% colnames(meta2)))
        next
      
      # Get calculations
      source_labels = meta1[cname][,1]
      transfer_labels = meta2[cname][,1]
      
      # Batch transfer
      predictions = knn(source_data, transfer_data, cl=source_labels, k=5)
      correct = sum(predictions == transfer_labels)
      cname_acc = (correct / length(transfer_labels))
      if (cname_acc != 0) {
        labels = append(labels, cname)
        acc = append(acc, cname_acc)
      }
      
      # Reverse
      # predictions = knn(transfer_data, source_data, cl=transfer_labels, k=5)
      # correct = sum(predictions == source_labels)
      # if (cname_acc != 0) {
      #   labels = append(labels, paste("reverse", cname))
      #   acc = append(acc, cname_acc)
      # }
    }
    
    # Barplot
    if (length(acc) < 1) {
      plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=0:1, ylim=0:1)
      return(NULL)
    }
    barplot(acc, main="Label Transfer Accuracy", names.arg=labels, xlim=c(0,1), horiz=T)
  })
  
  # Outputs
  output$content1 <- renderPlot({plot_alignment('sample1', 'Greens')})
  output$content2 <- renderPlot({plot_alignment('sample2', 'Oranges')})
  output$heatmap <- renderPlot({plot_heatmap()})
  # output$gci <- renderTable(data.frame(Gene=c(0,1,2), Contribution=c(0,1,0)))
  
  output$download <- downloadHandler(
    "aligned.csv",
    function(fname) {
      aligned = perform_alignment()
      aligned = cbind(aligned, cluster=get_clusters(aligned)$clusters)
      write.csv(aligned, fname)
    },
    "text/csv"
  )
}

# Run
shinyApp(ui=ui, server=server)

