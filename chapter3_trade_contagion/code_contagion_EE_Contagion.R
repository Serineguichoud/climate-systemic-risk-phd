contagion_model <- function(trade_matrix, node_weights, shocks, buffers) {
  # Validate inputs (all in one line to avoid line break issues)
  stopifnot(is.matrix(trade_matrix), length(node_weights)==nrow(trade_matrix), length(shocks)==nrow(trade_matrix), length(buffers)==nrow(trade_matrix), all(trade_matrix >= 0), all(node_weights >= 0))
  
  n <- nrow(trade_matrix)  # Simple assignment
  effective_shocks <- pmax(shocks - buffers, 0)
  loss_matrix <- matrix(0, n, n)
  
  for (i in 1:n) {
    if (effective_shocks[i] > 0) {
      loss_matrix[,i] <- trade_matrix[,i] * node_weights[i] * pmin(effective_shocks[i], 1)
    }
  }
  
  return(list(
    loss_matrix = loss_matrix,
    total_impact = rowSums(loss_matrix),
    surviving_buffers = pmax(buffers - shocks, 0)
  ))
}            
     
plot_large_contagion <- function(result, trade_matrix = NULL, top_n = 20) {
  # Set up graphics
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 1), mar = c(4, 6, 2, 1))
  
  # 1. Top-N Impact Summary (Barplot)
  top_impact <- head(sort(result$total_impact, decreasing = TRUE), top_n)
  barplot(top_impact,
          horiz = TRUE,
          las = 1,
          col = colorRampPalette(c("lightblue", "darkred"))(top_n),
          main = paste("Top", top_n, "Most Impacted Nodes"),
          xlab = "Total Contagion Impact")
  
  # 2. Aggregate Sector View (if you have sector labels)
  if (!is.null(rownames(trade_matrix))) {
    sectors <- substr(rownames(trade_matrix), 1, 3)  # Extract first 3 chars as sector
    sector_impact <- tapply(result$total_impact, sectors, sum)
    barplot(sort(sector_impact, decreasing = TRUE),
            cex.names = 0.7,
            main = "Impact by Sector",
            ylab = "Total Impact")
  }
  
  # 3. Simplified Network (Optional for large matrices)
  if (!is.null(trade_matrix) && requireNamespace("igraph", quietly = TRUE) && top_n <= 50) {
    # Create subnet with top nodes
    top_nodes <- names(top_impact)
    subnet <- igraph::graph_from_adjacency_matrix(
      trade_matrix[top_nodes, top_nodes],
      weighted = TRUE,
      mode = "directed"
    )
    # Plot
    set.seed(123)
    plot(subnet,
         edge.arrow.size = 0.3,
         vertex.size = sqrt(result$total_impact[top_nodes]) * 2,
         vertex.label.cex = 0.7,
         vertex.color = scales::alpha("red", 0.7),
         main = paste("Contagion Subnetwork (Top", top_n, "Nodes)"))
  }
}

###########

trade_net <- as.matrix(MatrixAll[-1])
dim(md_mat2)

rownames(trade_net) <- colnames(trade_net) <- MatrixAll$Flows


# Example data
result <- contagion_model(
  trade_matrix = trade_net,
  node_weights = Raw_data_with_all_damages$weights,
  shocks = Raw_data_with_all_damages$shock,
  buffers = Raw_data_with_all_damages$buffer
)

# Plot top 30 nodes and sectors
plot_large_contagion(result, large_trade_matrix, top_n = 30)


# Plot only non-zero losses
non_zero_losses <- result$loss_matrix[result$loss_matrix > 0]
hist(log10(non_zero_losses), main = "Distribution of Loss Magnitudes")

# Export for external tools
write.csv(result$loss_matrix, "loss_matrix.csv")

getwd()

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)
  
####last version

  plot_contagion_effects <- function(result, initial_shocks, node_names = NULL, top_n = 30) {
    # Validate inputs
    stopifnot(length(result$total_impact) == length(initial_shocks))
    
    # Get or create node names
    if (is.null(node_names)) {
      node_names <- if(!is.null(names(result$total_impact))) {
        names(result$total_impact)
      } else if(!is.null(rownames(result$loss_matrix))) {
        rownames(result$loss_matrix)
      } else {
        paste0("Node ", 1:length(result$total_impact))
      }
    }
    
    # Create data frame with guaranteed names
    df <- data.frame(
      Node = factor(node_names, levels = node_names[order(result$total_impact, decreasing = TRUE)]),
      Initial_Shock = initial_shocks,
      Additional_Impact = pmax(result$total_impact - initial_shocks, 0),
      Total_Impact = result$total_impact,
      stringsAsFactors = FALSE
    )
    
    # Select top nodes
    df_top <- head(df[order(-df$Total_Impact), ], top_n)
    
    # Create plot with improved labeling
    library(ggplot2)
    p <- ggplot(df_top, aes(x = reorder(Node, Total_Impact), y = Total_Impact)) +
      geom_col(aes(y = Initial_Shock, fill = "Initial Shock")) +
      geom_col(aes(y = Additional_Impact, fill = "Contagion Effect"), 
               position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = c("Initial Shock" = "#E41A1C", 
                                   "Contagion Effect" = "#377EB8")) +
      labs(title = "Contagion Impact Breakdown",
           x = "Node",
           y = "Impact Value",
           fill = "Impact Type") +
      coord_flip() +  # Horizontal bars for better name visibility
      theme_minimal() +
      theme(legend.position = "top",
            axis.text.y = element_text(size = 8))  # Adjust text size
    
    # Dynamic height based on number of nodes
    plot_height <- max(5, top_n * 0.3)
    ggsave("contagion_plot.png", p, width = 10, height = plot_height, units = "in")
    
    return(p)
  }

# Calculate initial effective shocks
# With explicit node names (recommended for 169 nodes)
node_names <- MatrixAll$Flows # Your 169 names

head(node_names)

# Or use existing names from your matrix
# node_names <- rownames(your_trade_matrix)

initial_effective_shocks=Raw_data_with_all_damages$shock*100
head(initial_effective_shocks)

# Generate plot (top 20 nodes)
plot_contagion_effects(result, initial_effective_shocks,node_names, top_n = 20)



