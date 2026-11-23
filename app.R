library(shiny)
library(tidyverse)
library(ggplot2)
library(umap)
library(NMF)
library(mclust)

# Load data
TBnanostring <- readRDS("TBnanostring.rds")
TBnanostring$TB_Status <- as.factor(TBnanostring$TB_Status)

gene_data <- TBnanostring %>%
  select(-TB_Status) %>%
  as.data.frame()

# PCA coordinates
pca_fit <- prcomp(gene_data, center = TRUE, scale. = FALSE)
pca_data <- as.data.frame(pca_fit$x[, 1:2])
names(pca_data) <- c("Dim1", "Dim2")
pca_data$TB_Status <- TBnanostring$TB_Status

# UMAP coordinates
set.seed(123)
umap_fit <- umap(
  gene_data,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean"
)
umap_data <- as.data.frame(umap_fit$layout)
names(umap_data) <- c("Dim1", "Dim2")
umap_data$TB_Status <- TBnanostring$TB_Status

# NMF coordinates
gene_matrix <- as.matrix(gene_data)
gene_minima <- apply(gene_matrix, 2, min, na.rm = TRUE)
gene_matrix_shifted <- sweep(gene_matrix, 2, gene_minima, FUN = "-")

# NMF input: genes in rows and samples in columns
set.seed(123)
nmf_fit <- nmf(
  t(gene_matrix_shifted),
  rank = 2,
  method = "brunet",
  nrun = 10,
  seed = 123
)

nmf_coordinates <- as.data.frame(t(coef(nmf_fit)))
names(nmf_coordinates) <- c("Dim1", "Dim2")
nmf_coordinates$TB_Status <- TBnanostring$TB_Status

# User interface
ui <- fluidPage(
  titlePanel("TB NanoString: Dimension Reduction and K-means"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "method",
        label = "Choose dimension-reduction method:",
        choices = c("PCA", "UMAP", "NMF"),
        selected = "PCA"
      ),
      sliderInput(
        inputId = "k",
        label = "Number of clusters (k):",
        min = 2,
        max = 10,
        value = 2,
        step = 1
      )
    ),
    mainPanel(
      plotOutput("cluster_plot"),
      h4("Cluster membership compared with TB status"),
      tableOutput("cluster_table"),
      h4("Adjusted Rand Index"),
      tableOutput("ari_table")
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  selected_data <- reactive({
    if (input$method == "PCA") {
      pca_data
    } else if (input$method == "UMAP") {
      umap_data
    } else {
      nmf_coordinates
    }
  })
  
  clustering <- reactive({
    coordinates <- selected_data() %>%
      select(Dim1, Dim2)
    
    set.seed(123)
    kmeans(
      coordinates,
      centers = input$k,
      nstart = 25
    )
  })
  
  output$cluster_plot <- renderPlot({
    current_data <- selected_data()
    km <- clustering()
    
    plot_data <- current_data %>%
      mutate(Cluster = factor(km$cluster))
    
    centers <- as.data.frame(km$centers)
    names(centers) <- c("Dim1", "Dim2")
    
    ggplot(
      plot_data,
      aes(x = Dim1, y = Dim2, color = Cluster, shape = TB_Status)
    ) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_point(
        data = centers,
        aes(x = Dim1, y = Dim2),
        inherit.aes = FALSE,
        shape = 4,
        size = 5,
        stroke = 1.5,
        color = "black"
      ) +
      labs(
        title = paste(input$method, "with K-means Clustering"),
        x = paste(input$method, "Dimension 1"),
        y = paste(input$method, "Dimension 2"),
        color = "Cluster",
        shape = "TB Status"
      ) +
      theme_minimal()
  })
  
  output$cluster_table <- renderTable({
    table(
      Cluster = clustering()$cluster,
      TB_Status = selected_data()$TB_Status
    )
  })
  
  output$ari_table <- renderTable({
    data.frame(
      Method = input$method,
      k = input$k,
      Adjusted_Rand_Index = adjustedRandIndex(
        clustering()$cluster,
        selected_data()$TB_Status
      )
    )
  })
}

shinyApp(ui = ui, server = server)