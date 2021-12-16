# NOTE: Crashes every other launch because of retriculate
# Libraries
library(ManiNetCluster)
library(plot3D)
library(plyr)
library(RColorBrewer)
library(shiny)

# UI
ui <- fluidPage(
  # Meta
  title="NLMA WebApp",
  
  # Main plot
  # asddf: Add interactive plotting
  h2("Result"),
  fluidRow(
    column(6, plotOutput("content1")),
    column(6, plotOutput("content2")),
  ),
  hr(),
  
  # Upload/Download UI
  # asdf: Multiple uploads at once
  # asdf: Label upload
  fluidRow(
    column(4,
      fileInput("mat1", NULL, buttonLabel="Matrix 1", multiple=FALSE),
      fileInput("mat2", NULL, buttonLabel="Matrix 2", multiple=FALSE),
      fileInput("corr", NULL, buttonLabel="Correspondence", multiple=FALSE),
    ),
    column(4,
      sliderInput("d", "d", min = 1, max = 20, value = 3),
      sliderInput("knn", "K_NN", min = 1, max = 20, value = 3), 
      sliderInput("kmed", "K_MED", min = 1, max = 20, value = 6),  
    ),
    column(4,
      downloadButton("download", "Download"),
    ),
  ),
  
)

# Server
server <- function(input, output, session) {
  output$files <- renderTable(input$files)
  
  # asddf: Make less redundant
  perform_alignment = reactive({
    mat1 <- input$mat1
    mat2 <- input$mat2
    corr <- input$corr
    
    if (is.null(mat1) || is.null(mat2)) {
      return(NULL)
    }
    
    mat1 <- as.matrix(read.csv(mat1$datapath, row.names=1))
    mat2 <- as.matrix(read.csv(mat2$datapath, row.names=1))
    
    # Use KNN as correspondence if applicable
    if (!is.null(corr)) {
      corr <- as.matrix(read.csv(corr$datapath, row.names=1))
      XY_corr = Correspondence(matrix=corr)
    }
    else if (dim(mat1)[1] == dim(mat2)[1]) {
      XY_corr = Correspondence(matrix=diag(dim(mat1)[1]))
    }
    else {
      return(NULL)
    }
    
    # Run NLMA
    ManiNetCluster(
      mat1,mat2,
      nameX='sample1',nameY='sample2',
      corr=XY_corr,
      d=as.integer(input$d),
      method='nonlinear manifold aln',
      # TODO: Dropdown menu
      # cca, lm
      k_NN=as.integer(input$knn),
      k_medoids=as.integer(input$kmed)
    )
  })
  
  plot_alignment = function(sample_label, color_str) {
    input$mat1
    df2 <- perform_alignment()
    df2$time = c(sel.meta1$time,sel.meta2$Days) #
    time.cols1 = colorRampPalette(brewer.pal(n=12,color_str))(12) #
    res = data.frame(df2[df2$data==sample_label,]) # Change this line
    res0 = data.frame(df2)
    s3d<-scatter3D(x=res[,3],y=res[,4],z=res[,5],
                   colvar=as.numeric(mapvalues(res$time,names(table(res$time)),c(2:13))),
                   col = c(time.cols1),pch=c(16,17)[as.numeric(as.factor(res$data))],
                   colkey=F, theta = 300, phi = 30,cex=2,
                   xlim=c(min(res0$Val0),max(res0$Val0)),
                   ylim=c(min(res0$Val1),max(res0$Val1)),
                   zlim=c(min(res0$Val2),max(res0$Val2)))
    # asddf: Revise legend
    #legend("top", legend = levels(as.factor(res$data)), pch = c(16, 17),inset = -0.1, xpd = TRUE, horiz = TRUE)
    #legend("right", legend = levels(as.factor(res$time)), col = c(time.cols1),pch=16,inset = 0.1, xpd = TRUE, horiz = F,cex=2)
  }
  
  output$content1 <- renderPlot({plot_alignment('sample1', 'Greens')})
  output$content2 <- renderPlot({plot_alignment('sample2', 'Oranges')})
  
  output$download <- downloadHandler(
    "aligned.csv",
    function(fname) {
      write.csv(perform_alignment(), fname, row.names=F)
    },
    "text/csv"
  )
}

# Run
shinyApp(ui=ui, server=server)
