# NOTE: Crashes every other launch because of retriculate
# Libraries
library(class)
library(cluster)
library(dtw)
library(htmlwidgets)
library(pdist)
library(plot3D)
library(plotly)
library(plyr)
library(purrr)
library(RColorBrewer)
library(shiny)
library(shinyjs)

# Set up reticulate (https://github.com/ranikay/shiny-reticulate-app)
# Requires python 3
PYTHON_DEPENDENCIES = c('pip', 'matplotlib', 'pandas', 'scipy', 'sklearn')
virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
python_path = Sys.getenv('PYTHON_PATH')

# Create virtual env and install dependencies
# https://github.com/namtk/ManiNetCluster
reticulate::virtualenv_create(envname=virtualenv_dir, python=python_path)
reticulate::virtualenv_install(virtualenv_dir, packages=PYTHON_DEPENDENCIES, ignore_installed=FALSE)
reticulate::use_virtualenv(virtualenv_dir, required=T)

# Import python functions
# devtools::install_github("namtk/ManiNetCluster")
reticulate::source_python('./maniNetCluster/pyManifold.py')

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
boma_alignments = c(
  "dtw"
  # asddf: Add KNN
)
names(boma_alignments) = c(
  "Dynamic Time Warping"
)
datasets = c(
  "None",
  # Brain
  "./data/Sample1/",
  "./data/LiScience2018/",
  "./data/NowakowskiScience2017/",
  # "./data/TrevinoCell2021/",
  # "./data/BhaduriNature2020Brain/",
  # Organoid
  "./data/Sample2/",
  "./data/GordonNature2021/",
  "./data/BireyNature2017/"
  # "./data/BhaduriNature2020Org/"
)
names(datasets) = c(
  "Uploaded Dataset (Below)",
  # Brain
  "Sample 1 (490 Bulk x 903 Genes)",
  "Li et al., Science, 2018 (460 Bulk x 559 Genes)",
  "Nowakowski, Science, 2017 (4261 x 56864 Cells)",
  # "Trevino, Cell, 2021 ()",
  # "Bhaduri, Nature, 2020, Brain ()",
  # Organoid
  "Sample 2 (497 Bulk x 903 Genes)",
  "Gordon et al., Nature, 2021 (62 Bulk x 559 Genes)",
  "Birey, Nature, 2017 (11838 Cells x 21092 Genes)"
  # "Bhaduri, Nature, 2020, Organoid ()"
)
brain_datasets = names(datasets)[c(1, 2, 3, 4)]
organoid_datasets = names(datasets)[c(1, 5, 6, 7)]
brain_default = names(datasets)[3]
organoid_default = names(datasets)[6]

# Options
options(shiny.maxRequestSize=0)

# UI
ui <- fluidPage(
  # Meta
  title="BOMA App",
  
  # Main plot
  useShinyjs(),
  h2("Brain Organoid Manifold Alignment"),
  # textOutput("warnings"),
  fluidRow(
    column(4, img(src='figure1.png', align = "left", width="100%")),
    column(4, plotlyOutput("content1")),
    column(4, plotlyOutput("content2")),
  ),
  hr(),
  
  # UI
  fluidRow(
    column(4,
      h2("Upload"),
      strong("Row and column labels are expected in all files.  All files are
        also assumed to be ordered similarly (by sample).  Default datasets and source code can be
        found", a("here", href="https://github.com/Oafish1/ManiNetCluster-Visualization/tree/main/data")),
      
      hr(),
      fluidRow(
        column(6,
          selectInput("custom_mat1", "Use existing brain dataset",
                     choices=brain_datasets, selected=brain_default),
          div(style = "margin-top: -10px"),
          fileInput("mat1", NULL, buttonLabel="First Dataset"),
          div(style = "margin-top: -30px"),
          fileInput("meta1", NULL, buttonLabel="First Metadata"),
        ),
        column(6,
          selectInput("custom_mat2", "Use existing organoid dataset",
                      choices=organoid_datasets, selected=organoid_default),
          div(style = "margin-top: -10px"),
          fileInput("mat2", NULL, buttonLabel="Second Dataset"),
          div(style = "margin-top: -30px"),
          fileInput("meta2", NULL, buttonLabel="Second Metadata"),
        ),
      ),
      fluidRow(
        column(4, p("Datasets should be of dimension (samples x features) in CSV format. ",
                    "Assumed to be pre-normalized. ", strong("For use with BOMA,",
                    "number of features in each dataset should be the same."))),
        column(2, img(src='data.png', align = "left", width="100%")),
        column(4, p("Metadata should be ordered similarly to dataset data in CSV format. ",
                    "Two columns may be chosen as the basis for point coloring in the scatterplot and",
                    "BOMA ordering. ",
                    "All other columns are assumed to be categorical and are used for error reporting.")),
        column(2, img(src='metadata.png', align = "left", width="100%")),
      ),
      
      hr(),
      fileInput("corr", NULL, buttonLabel="Correspondence", multiple=FALSE),
      div(style = "margin-top: -20px"),
      fluidRow(
        column(8,
        p("Correspondence matrix (samples1 x samples2) should be of CSV format. ",
          "Meant to represent inter-dataset correspondence and can be calculated",
          "in a variety of ways.  If aligned and no CSV is provided, will default",
          "to the identity matrix. ", strong("BOMA does not use this."))),
      column(4, img(src='corr.png', align = "left", width="50%")),
      ),
    ),
    column(4,
      h2("Alignment"),
      checkboxInput("use_boma", "Use BOMA", value=T),
      fluidRow(
        column(6,
               selectInput("boma_method", "Global Alignment", choices=names(boma_alignments), selected="DTW")),
        column(6,
               selectInput("boma_col", "Sample Ordering", choices=NULL)),
      ),
      
      hr(),
      selectInput("method", "Local Alignment", choices=names(alignments), selected="Non-Linear Manifold Alignment"),
      fluidRow(
        column(4,
               sliderInput("d", "Dimensions", min = 1, max = 20, value = 3)),
        column(4,
               sliderInput("knn", "Nearest Neighbors", min = 0, max = 20, value = 3), ),
        column(4,
               sliderInput("kmed", "Medoids", min = 0, max = 20, value = 6),),
      ),
      
      hr(),
      fluidRow(
        column(6,
               downloadButton("download", "Download Aligned Data"),),
        column(6,
               p("Download data as a single CSV using the chosen parameters.  Includes clusters."),),
      ),
      
      hr(),
      plotOutput("statistics", height=200),
      em("Zero values are assumed to be errant and are excluded. ",
         "Labels are extracted from columns present in both metadata files. "),
    ),
    column(4,
      h2("Visualization"),
      checkboxInput("show_clusters", "Show Clusters in Plot", value=F),
      fluidRow(
        column(6,
               selectInput("color_col", "Cluster Coloring", choices=NULL)),
        column(6,
               selectInput("cluster_method", "Clustering Method", choices=c("K-Means", "PAM"), selected="PAM")),
      ),
      sliderInput("num_clusters", "Clusters", min = 1, max = 14, value = 4),
      
      # h2("Gene Information"),
      # em("Placeholder:"),
      # tableOutput("gci"),
      
      hr(),
      fluidRow(
        column(9, plotOutput("heatmap")),
        column(3, plotOutput("colorbar")),
      ),
    ),
  ),
  
)

# Server
server <- function(input, output, session) {
  get_method <- function(method_code) {
    for (i in 1:length(alignments)) {
      if (method_code == names(alignments)[i])
        return(alignments[i])
    }
    return(NULL)
  }
  
  
  get_boma_method <- function(method_code) {
    for (i in 1:length(boma_alignments)) {
      if (method_code == names(boma_alignments)[i])
        return(boma_alignments[i])
    }
    return(NULL)
  }
  
  
  get_custom_dataset <- function(method_code) {
    for (i in 1:length(datasets)) {
      if (method_code == names(datasets)[i])
        return(datasets[i])
    }
    return(NULL)
  }
  
  
  get_clusters <- function(df, default_color="green", normal_colors=F, num_colors=8) {
    # Generate clusters / Get colors
    working_data = df[,3:(2+input$d)]
    # Max 7 colors
    colors = c("blue","green","yellow","orange","red","pink","purple")
    
    if (normal_colors)
      return (list("colors"=colorRampPalette(brewer.pal(n=9, default_color))(num_colors)))
    else if (input$cluster_method == "K-Means")
      clusters = kmeans(working_data, input$num_clusters)$cluster
    else if (input$cluster_method == "PAM")
      clusters = pam(working_data, input$num_clusters)$clustering
    else if (input$cluster_method == "Spectral") {
      s1 <- as.matrix(data.frame(df[df$data=="sample1",])[,3:(2+input$d)])
      s2 <- as.matrix(data.frame(df[df$data=="sample2",])[,3:(2+input$d)])
      dist = as.matrix(pdist::pdist(s1, s2))
    }
    else
      validate(need(F, "Invalid clustering method"))
    return (list("clusters"=clusters, "colors"=colors[clusters]))
  }
  
  
  get_raw_clusters <- function(df) {
    # Generate clusters / Get colors
    working_data = df
    # Max 7 colors
    colors = c("blue","green","yellow","orange","red","pink","purple")
    
    if (input$cluster_method == "K-Means")
      clusters = kmeans(working_data, input$num_clusters)$cluster
    else if (input$cluster_method == "PAM")
      clusters = pam(working_data, input$num_clusters)$clustering
    else
      validate(need(F, "Invalid clustering method"))
    return (clusters)
  }
  
  
  plot_alignment = function(sample_label, default_color, use_legend=F) {
    dimensions = refresh_dimensions()
    validate(need(dimensions >= 3, paste(
      "Must have dimensions \u22653 to plot, currently have", dimensions)))
    validate(need(input$color_col != "", "Coloring column must be non-empty"))
    
    # Get inputs
    meta = get_meta()
    if(is.null(meta))
      return(NULL)
    meta1 = meta$meta1
    meta2 = meta$meta2
    
    # Plot code from paper
    df2 <- perform_alignment()
    if(is.null(df2))
      return(NULL)
    # Assumes ordering from NLMA
    if (is.null(tryCatch({meta1[,input$color_col]; meta2[,input$color_col]},
                  error=function(cond) {return(NULL)})))
      return(NULL)
    df2$time = c(meta1[,input$color_col],meta2[,input$color_col])
    res = data.frame(df2[df2$data==sample_label,])
    res0 = data.frame(df2)
    clusters = get_clusters(df2, default_color, !input$show_clusters, length(unique(res$time)))
    time.cols1 = clusters$colors
    if (!input$show_clusters) {
      par(mar=c(5.1,4.1,4.1,4.1), xpd=T)
      time.cols1 = time.cols1[ as.numeric(mapvalues(res$time,unique(res$time),c(1:length(unique(res$time))))) ]
      data = data.frame(x=res[,3], y=res[,4], z=res[,5], color=time.cols1, label=as.factor(res$time))
      
      legend_title = input$color_col
    } else {
      time.cols1 = time.cols1[df2$data==sample_label]
      labels = paste("Cluster ", clusters$clusters[df2$data==sample_label], sep="")
      data = data.frame(x=res[,3], y=res[,4], z=res[,5], color=time.cols1, label=labels)
      legend_title = "Clusters"
    }
    data = data[order(data$label), ]
    
    s3d = plot_ly(
      type = "scatter3d",
      mode = "markers",
      data = data,
      x = ~x, y = ~y, z = ~z,
      marker = list(bgcolor="#e5e5e5", color=~color),
      split=~label
    )
    # https://stackoverflow.com/a/66117098 for continually rotating
    s3d %>% layout(
      title=paste("Dataset", substr(sample_label, nchar(sample_label), nchar(sample_label))),
      legend=list(title=list(text=legend_title)),
      scene=list(
        xaxis=list(range=c(min(res0$Val0), max(res0$Val0))),
        yaxis=list(range=c(min(res0$Val1), max(res0$Val1))),
        zaxis=list(range=c(min(res0$Val2), max(res0$Val2))),
        aspectmode='cube',
        camera=list(
          eye = list(
            x = 1.25,
            y = 1.25,
            z = 1.25),
          center = list(
            x = 0,
            y = 0,
            z = 0
        ))
    ))
  }
  
  
  get_dist <- function() {
    aligned <- perform_alignment()
    if(is.null(aligned))
      return(NULL)
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
    dist = as.matrix(pdist::pdist(s1, s2))
    
    return(list("dist"=dist, "rsc"=rsc, "csc"=csc, "clusters"=clusters))
  }
  
  
  plot_heatmap <- function() {
    distance <- get_dist()
    if(is.null(distance))
      return(NULL)
    csc = distance$csc
    rsc = distance$rsc
    dist = distance$dist
    clusters = distance$clusters
    
    # layout(t(1:2),widths=c(6,1))
    # par(mfrow=c(1,2), mar=c(4, 4, 1, .5))
    ramp = colorRampPalette(rev(brewer.pal(9, "Greys")))
    heatmap(dist, main="Alignment Distance", xlab="Dataset 2", ylab="Dataset 1", RowSideColors=rsc, ColSideColors=csc, Colv=NA, Rowv=NA, labCol=F, labRow=F, col=ramp(256))
    # legend("topright", legend=c("close", "mid", "far"), inset = c(0,0), fill=ramp(3))
    legend("bottomright", legend = paste("Cluster ", unique(clusters$clusters), sep=""), col = unique(clusters$colors), pch=16, inset = c(0,0), xpd = TRUE, horiz = F)

    # par(mar=c(5, 1, 5, 2.5))
    # n = 20
    # color=ramp(n)  # rainbow(max(dist), s = 1, v = 1, start = 0, end = 1, alpha = 1)
    # image(y=0:n,z=t(0:n), col=color[1:9], axes=FALSE, main="Slope", cex.main=.8)
    # axis(4,cex.axis=0.8,mgp=c(0,.5,0))
  }
  
  
  plot_colorbar <- function() {  
    distance <- get_dist()
    if(is.null(distance))
      return(NULL)
    csc = distance$csc
    rsc = distance$rsc
    dist = distance$dist
    
    n =30
    color = colorRampPalette(rev(brewer.pal(9, "Greys")))(n)
    par(mar=c(7,0,7,4))  # b, l, t, r
    ax_range = seq(from=min(dist), to=max(dist), by=(max(dist) - min(dist))/n)
    image(y=ax_range,z=t(ax_range), col=color[1:n], axes=FALSE, main="Distance\nLegend", ylab="", cex.main=.8)
    axis(4,cex.axis=0.8,mgp=c(0,.5,0))
  }
  
  
  get_meta <- reactive({
    meta1 <- input$meta1
    meta2 <- input$meta2
    
    # Load uploaded datasets
    if (!is.null(meta1))
      meta1 = read.csv(meta1$datapath, row.names=1)
    if (!is.null(meta2))
      meta2 = read.csv(meta2$datapath, row.names=1)
    
    # Load custom datasets
    custom_dataset1 = get_custom_dataset(input$custom_mat1)
    if (custom_dataset1 != "None")
      meta1 = read.csv(paste(custom_dataset1, "meta.csv", sep=""), row.names=1)
    custom_dataset2 = get_custom_dataset(input$custom_mat2)
    if (custom_dataset2 != "None")
      meta2 = read.csv(paste(custom_dataset2, "meta.csv", sep=""), row.names=1)
    
    data = get_data()
    if(is.null(data))
      return(NULL)
    if( (dim(data$mat1)[1] == dim(meta1)[1]) && (dim(data$mat2)[1] == dim(meta2)[1]) )
      return(list("meta1"=meta1, "meta2"=meta2))
    showModal(modalDialog("Length of meta and data must match", footer=NULL, easyClose=T))
    return(NULL)
  })
  
  
  get_data <- reactive({
    mat1 <- input$mat1
    mat2 <- input$mat2
    corr <- input$corr
    
    # Load uploaded datasets
    if (!is.null(mat1))
      mat1 <- as.matrix(read.csv(mat1$datapath, row.names=1))
    if (!is.null(mat2))
      mat2 <- as.matrix(read.csv(mat2$datapath, row.names=1))
    
    # Load custom datasets
    custom_dataset1 = get_custom_dataset(input$custom_mat1)
    if (custom_dataset1 != "None")
      mat1 = read.csv(paste(custom_dataset1, "mat.csv", sep=""), row.names=1)
    custom_dataset2 = get_custom_dataset(input$custom_mat2)
    if (custom_dataset2 != "None")
      mat2 = read.csv(paste(custom_dataset2, "mat.csv", sep=""), row.names=1)
    
    # Exception handling
    if (is.null(mat1) || is.null(mat2)) {
      removeModal()
      showModal(modalDialog("No dataset uploaded", footer=NULL, easyClose=T))
      return(NULL)
    }
    
    # Make sure datasets are compatible
    if (dim(mat1)[2] != dim(mat2)[2]) {
      removeModal()
      showModal(modalDialog(paste("Datasets do not have the same number of features,", dim(mat1)[2], "vs", dim(mat2)[2]),
                            footer=NULL, easyClose=T))
      return(NULL)
    }
      
    # Calculate correspondence if not provided
    if (!is.null(corr))
      corr <- as.matrix(read.csv(corr$datapath, row.names=1))
    else if (dim(mat1)[1] == dim(mat2)[1])
      corr = diag(dim(mat1)[1])
    else if (!input$use_boma) {
      removeModal()
      showModal(modalDialog("Correspondence matrix must be provided unless BOMA is used",
                            footer=NULL, easyClose=T))
      return(NULL)
    } else
      corr = NULL
    
    return(list("mat1"=mat1, "mat2"=mat2, "corr"=corr))
  })
  
  
  refresh_dimensions <- reactive({
    # Change sliders based on selection
    method = get_method(input$method)
    to_show = switch(
      get_method(input$method),
      "cca" = c(),
      "linear manifold" = c("knn"),
      "manifold warping" = c("knn"),
      "nonlinear manifold aln" = c("knn", "kmed"),
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
    
    # Disable corr upload with BOMA
    # if (input$use_boma)
    #   shinyjs::disable(id="corr")
    # else
    #   shinyjs::enable(id="corr")
    
    return(dimensions)
  })
  observeEvent(input$method, refresh_dimensions())
  
  
  conditional_sliders <- reactive({
    # Conditional sliders
    if (input$knn == 0) {
      shinyjs::enable(id="kmed")
      shinyjs::disable(id="d")
    }
    else {
      shinyjs::disable(id="kmed")
      shinyjs::enable(id="d")
    }
  })
  # observeEvent(input$knn, conditional_sliders())
  
  
  meta_col_choices <- reactive({
    meta = get_meta()
    if(is.null(meta))
      return(NULL)
    meta1 = meta$meta1
    meta2 = meta$meta2
    
    color_cols = c()
    for(col in colnames(meta1))
      if(col %in% colnames(meta2))
        color_cols = append(color_cols, col)
    return(color_cols)
  })
  meta_col_choices_none <- reactive({
    meta = get_meta()
    if(is.null(meta))
      return(NULL)
    meta1 = meta$meta1
    meta2 = meta$meta2
    
    color_cols = c("None")
    for(col in colnames(meta1))
      if(col %in% colnames(meta2))
        color_cols = append(color_cols, col)
    return(color_cols)
  })
  default_col_selection <- reactive({
    # Pretty hacky, but inconsequential
    if(identical(meta_col_choices(), c("Days", "time", "PAM.Cluster"))
       || identical(meta_col_choices(), c("None", "Days", "time", "PAM.Cluster")))
      return("time")
    return(NULL)
  })
  observeEvent(meta_col_choices(), updateSelectInput(session, "color_col", choices=meta_col_choices(), selected=default_col_selection()))
  observeEvent(meta_col_choices_none(), updateSelectInput(session, "boma_col", choices=meta_col_choices_none(), selected=default_col_selection()))
  
  
  conditional_colors <- reactive({
    if (input$show_clusters) {
      shinyjs::disable(id="color_col")
    }
    else {
      shinyjs::enable(id="color_col")
    }
  })
  observeEvent(input$show_clusters, conditional_colors())
  
  
  conditional_boma_clustering <- reactive({
    if (input$use_boma) {
      shinyjs::enable(id="boma_method")
      shinyjs::enable(id="boma_col")
    }
    else {
      shinyjs::disable(id="boma_method")
      shinyjs::disable(id="boma_col")
    }
  })
  observeEvent(input$use_boma, conditional_boma_clustering())
    
  
  perform_alignment <- reactive({
    # https://stackoverflow.com/a/52741787
    data = get_data()
    if(is.null(data))
      return(NULL)
    mat1 = data$mat1
    mat2 = data$mat2
    corr = data$corr
    
    # Get vars
    method = get_method(input$method)
    boma_method = get_boma_method(input$boma_method)
    
    # Perform Alignment
    showModal(modalDialog("Calculating alignment, could take a few minutes...", footer=NULL, easyClose=T))
    if(input$use_boma) {
      # Load metadata
      meta = get_meta()
      if(is.null(meta))
        return(NULL)
      meta1 = meta$meta1
      meta2 = meta$meta2
      
      # Reorder data
      if (input$boma_col == "")
        return(NULL)
      else if (input$boma_col != "None") {
        mat1 = mat1[order(meta1[,input$boma_col]),]
        mat2 = mat2[order(meta2[,input$boma_col]),]
      }
      
      # Global alignment
      if (boma_method == "dtw") {
        corr = cor(t(mat1), t(mat2))
        dist = 1 / (1+corr)
        # Asymmetric step means that switching mat1 <-> mat2 might not provide same result
        dtw_out = dtw::dtw(dist, step.pattern=dtw::asymmetric, open.begin=T, open.end=T)
        pairs_x = dtw_out$index1
        pairs_y = dtw_out$index2
        w = matrix(0, nrow(mat1), nrow(mat2))
        for (i in 1:length(pairs_x))
          w[pairs_x[i], pairs_y[i]] = 1
      }
      else {
        removeModal()
        showModal(modalDialog("Unknown global alignment type selected", footer=NULL, easyClose=T))
        return(NULL)
      }
      
      # Local alignment
      aligned = ManiNetCluster(
        mat1,mat2,
        nameX='sample1',nameY='sample2',
        corr=Correspondence(matrix=w),
        d=as.integer(input$d),
        method=method,
        k_NN=as.integer(input$knn),
        k_medoids=as.integer(input$kmed)
      )
      removeModal()
    } else {
      XY_corr = Correspondence(matrix=corr)
      aligned = ManiNetCluster(
        mat1,mat2,
        nameX='sample1',nameY='sample2',
        corr=XY_corr,
        d=as.integer(input$d),
        method=method,
        k_NN=as.integer(input$knn),
        k_medoids=as.integer(3),
        k_medoids=as.integer(input$kmed)
      )
      removeModal()
    }
    
    return(aligned)
  })
  
  
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
  
  
  output$statistics <- renderPlot({
    # Accuracy metrics
    # Get inputs
    data = get_data()
    if(is.null(data))
      return(NULL)
    corr = data$corr
    
    meta = get_meta()
    if(is.null(meta))
      return(NULL)
    meta1 = meta$meta1
    meta2 = meta$meta2
    
    # Prerequisite data
    aligned_data = perform_alignment()
    if(is.null(aligned_data))
      return(NULL)
    source_data = data.frame(aligned_data[aligned_data$data=='sample1',])[,3:(2+input$d)]
    transfer_data = data.frame(aligned_data[aligned_data$data=='sample2',])[,3:(2+input$d)]
    
    # Iterate through labels
    labels = c()
    acc = c()
    for (cname in colnames(meta1)) {
      if (!(cname %in% colnames(meta1) && cname %in% colnames(meta2)))
        next
      
      # Get calculations
      source_labels = meta1[,cname]
      transfer_labels = meta2[,cname]
      
      # Batch transfer
      predictions = tryCatch({knn(source_data, transfer_data, cl=source_labels, k=5)}, error=function(cond) {NULL})
      if (is.null(predictions))
        next
      correct = sum(predictions == transfer_labels)
      cname_acc = (correct / length(transfer_labels))
      if (cname_acc != 0) {
        labels = append(labels, paste("Label Transfer Acc.\n", cname))
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
    
    
    # Pairwise Distance (Marin)
    if (!is.null(corr)) {
      distance_sum = 0
      for (i in 1:nrow(source_data)) {
        for (j in 1:nrow(transfer_data)) {
          if(corr[i,j] == 1) {
            distance_sum = distance_sum + sqrt(sum( (source_data[i,] - transfer_data[j,])^2 ))
          }
        }
      }
      
      if (sum(corr == 1) > 0) {
        pairwise_dist = distance_sum / sum(corr == 1)
        labels = append(labels, 'Pairwise Distance')
        acc = append(acc, pairwise_dist)
      }
    }
    
    # Barplot
    if (length(acc) < 1) {
      plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=0:1, ylim=0:1)
      return(NULL)
    }
    par(mai=c(1,2,1,1))
    bar = barplot(acc, main="Alignment Statistics", names.arg=labels, xlim=c(0,1), horiz=T, las=1)
    # https://stackoverflow.com/a/59658285
    text(max(acc - .3, .2), bar, round(acc,3), font=2)
  })
  
  
  output$content1 <- renderPlotly({plot_alignment('sample1', 'Greens')})
  output$content2 <- renderPlotly({plot_alignment('sample2', 'Oranges')})
  output$heatmap <- renderPlot({plot_heatmap()})
  output$colorbar <- renderPlot({plot_colorbar()})
  # output$gci <- renderTable(data.frame(Gene=c(0,1,2), Contribution=c(0,1,0)))
  
  
  output$download <- downloadHandler("aligned.csv", function(fname) {
      aligned = perform_alignment()
      aligned = cbind(aligned, cluster=get_clusters(aligned)$clusters)
      write.csv(aligned, fname)
    }, "text/csv")
  # traceback()
}

# Run
shinyApp(ui=ui, server=server)

