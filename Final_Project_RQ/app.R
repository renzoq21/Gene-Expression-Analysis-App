# Increase file upload limit
options(shiny.maxRequestSize = 50 * 1024^2)  # Allow uploads up to 50 MB
# Load required libraries
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("Huntington's Disease RNASeq Analysis App"),
  
  # Tabs for each component
  tabsetPanel(
    
    # Tab 1: Sample Information Exploration
    tabPanel("Sample Information",
             sidebarLayout(
               sidebarPanel(
                 fileInput("sample_file", "Upload Sample Information (CSV or TXT):", accept = c(".csv", ".txt"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", verbatimTextOutput("sample_summary")),
                   tabPanel("Table", DTOutput("sample_table")),
                   tabPanel("Plots", 
                            uiOutput("plot_controls"),
                            plotOutput("sample_plot"))
                 )
               )
             )
    ),
    
    # Tab 2: Counts Matrix Exploration
    tabPanel("Counts Matrix",
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_file", "Upload Counts Matrix (CSV or TXT):", accept = c(".csv", ".txt")),
                 sliderInput("var_filter", "Minimum Variance Percentile:", min = 0, max = 100, value = 50),
                 sliderInput("zero_filter", "Minimum Non-Zero Samples:", min = 0, max = 100, value = 10),
                 actionButton("apply_filters", "Apply Filters")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", verbatimTextOutput("counts_summary")),
                   tabPanel("Heatmap", plotOutput("heatmap_plot")),
                   tabPanel("PCA", plotOutput("pca_plot"))
                 )
               )
             )
    ),
    
    # Tab 3: Differential Expression
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 fileInput("de_file", "Upload DE Results (CSV or TXT):", accept = c(".csv", ".txt"))
               ),
               mainPanel(
                 DTOutput("de_table"),
                 plotOutput("volcano_plot")
               )
             )
    ),
    
    # Gene Expression Tab
    tabPanel("Gene Expression",
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_file_gene", "Upload Counts Matrix (CSV or TXT):", accept = c(".csv", ".txt")),
                 fileInput("sample_file_gene", "Upload Sample Information (CSV or TXT):", accept = c(".csv", ".txt")),
                 selectInput("plot_gene", "Select Gene:", choices = NULL),
                 selectInput("grouping_var", "Select Grouping Variable:", choices = NULL),
                 actionButton("plot_button", "Generate Plot")
               ),
               mainPanel(
                 plotOutput("gene_plot")
               )
             )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Helper function to read both CSV and TXT files
  read_file <- function(file_path) {
    if (grepl("\\.csv$", file_path)) {
      read.csv(file_path, row.names = 1)
    } else if (grepl("\\.txt$", file_path)) {
      read.delim(file_path, row.names = 1)
    } else {
      stop("Unsupported file format")
    }
  }
  
  ### Tab 1: Sample Information Exploration ###
  sample_data <- reactive({
    req(input$sample_file)
    read_file(input$sample_file$datapath)
  })
  
  output$sample_summary <- renderPrint({
    req(sample_data())
    summary(sample_data())
  })
  
  output$sample_table <- renderDT({
    req(sample_data())
    datatable(sample_data())
  })
  
  output$plot_controls <- renderUI({
    req(sample_data())
    valid_columns <- setdiff(names(sample_data()), c("X", "symbol"))
    selectInput("plot_var", "Select Variable to Plot:", choices = valid_columns)
  })
  
  output$sample_plot <- renderPlot({
    req(sample_data(), input$plot_var)
    data <- sample_data()[[input$plot_var]]
    if (is.numeric(data)) {
      ggplot(sample_data(), aes_string(x = input$plot_var)) + 
        geom_histogram(bins = 20, fill = "steelblue", color = "black") +
        theme_minimal()
    } else {
      ggplot(sample_data(), aes_string(x = input$plot_var)) + 
        geom_bar(fill = "steelblue", color = "black") +
        theme_minimal()
    }
  })
  
  ### Tab 2: Counts Matrix Exploration ###
  counts_data <- reactive({
    req(input$counts_file)
    read_file(input$counts_file$datapath)
  })
  
  filtered_data <- eventReactive(input$apply_filters, {
    req(counts_data())
    numeric_data <- counts_data()[, sapply(counts_data(), is.numeric)]
    variances <- apply(numeric_data, 1, var)
    non_zero_counts <- rowSums(numeric_data > 0)
    numeric_data[variances >= quantile(variances, input$var_filter / 100) & 
                   non_zero_counts >= input$zero_filter, ]
  })
  
  output$counts_summary <- renderPrint({
    req(filtered_data())
    cat("Samples:", ncol(filtered_data()), "\nGenes:", nrow(filtered_data()))
  })
  
  output$heatmap_plot <- renderPlot({
    req(filtered_data())
    heatmap(as.matrix(filtered_data()[1:50, ]), scale = "row", main = "Filtered Gene Heatmap")
  })
  
  output$pca_plot <- renderPlot({
    req(filtered_data())
    pca <- prcomp(t(filtered_data()))
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
    ggplot(pca_df, aes(PC1, PC2)) + 
      geom_point(color = "blue") + 
      labs(title = "PCA of Filtered Gene Counts", x = "PC1", y = "PC2") +
      theme_minimal()
  })
  
  ### Tab 3: Differential Expression ###
  de_data <- reactive({
    req(input$de_file)
    read_file(input$de_file$datapath)
  })
  
  output$de_table <- renderDT({
    req(de_data())
    datatable(de_data())
  })
  
  output$volcano_plot <- renderPlot({
    req(de_data())
    ggplot(de_data(), aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = padj < 0.05)) +
      labs(title = "Volcano Plot") + theme_minimal()
  })
  
  ### Tab 4: Gene Expression Visualization ###
  gene_counts_data <- reactive({
    req(input$counts_file_gene)
    counts <- clean_counts_matrix(input$counts_file_gene$datapath)
    print("Cleaned Counts Matrix Row Names:")
    print(head(rownames(counts)))
    return(counts)
  })
  
  gene_sample_data <- reactive({
    req(input$sample_file_gene)
    sample_info <- clean_sample_info(input$sample_file_gene$datapath)
    print("Cleaned Sample Information SampleID:")
    print(head(sample_info$SampleID))
    return(sample_info)
  })
  
  observe({
    req(gene_counts_data())
    updateSelectInput(session, "plot_gene", choices = rownames(gene_counts_data()))
  })
  
  observe({
    req(gene_sample_data())
    updateSelectInput(session, "grouping_var", choices = colnames(gene_sample_data()))
  })
  
  output$gene_plot <- renderPlot({
    req(input$plot_button, gene_counts_data(), gene_sample_data())
    
    # Load cleaned data
    counts <- gene_counts_data()
    sample_info <- gene_sample_data()
    selected_gene <- trimws(tolower(input$plot_gene))  # Ensure selected gene is clean
    
    # Debugging: Print selected gene and check matching
    print("Selected Gene:")
    print(selected_gene)
    print("Row Names in Counts Matrix:")
    print(head(rownames(counts)))
    print("SampleID in Sample Information:")
    print(head(sample_info$SampleID))
    
    # Check if the selected gene exists in both datasets
    if (!(selected_gene %in% rownames(counts))) {
      stop("Error: Selected gene not found in counts matrix.")
    }
    if (!(selected_gene %in% sample_info$SampleID)) {
      stop("Error: Selected gene not found in sample information file.")
    }
    
    # Extract expression values for the selected gene
    expression_values <- as.numeric(counts[selected_gene, ])
    
    # Match selected gene in the sample information
    gene_info <- sample_info[sample_info$SampleID == selected_gene, ]
    
    # Prepare plot data
    plot_data <- data.frame(
      Sample = colnames(counts),
      Expression = expression_values,
      Group = rep(gene_info[[input$grouping_var]], length(expression_values))
    )
    
    # Generate Boxplot
    ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot() +
      labs(
        title = paste("Expression of Gene:", selected_gene),
        x = input$grouping_var, y = "Expression Level"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

# Run the app
shinyApp(ui, server)
