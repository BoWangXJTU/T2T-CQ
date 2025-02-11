#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# Load required libraries (ensure they are manually installed)
library(ggplot2)
library(optparse)
library(dplyr)

# Command-line argument parsing
option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Path to the default input data file [e.g., data.txt]", metavar = "file"),
  make_option(c("--input_A"), type = "character", default = NULL,
              help = "Path to the first input data file (A) [e.g., data_A.txt]", metavar = "file"),
  make_option(c("--input_B"), type = "character", default = NULL,
              help = "Path to the second input data file (B) [e.g., data_B.txt]", metavar = "file")
)

# Add script usage example
description_text <- "Example usage:\n Rscript violin_plot.R --input data.txt --input_A data_A.txt --input_B data_B.txt"

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

# Check if at least one input file is provided
if (is.null(opt$input) && is.null(opt$input_A) && is.null(opt$input_B)) {
  print_help(opt_parser)
  stop("Error: Please provide at least one input file using --input, --input_A, or --input_B arguments.\n")
}

# Read data
data_list <- list()

if (!is.null(opt$input)) {
  data <- read.table(opt$input, header = FALSE, sep = "\t", col.names = c("ID", "Length"))
  data$Group <- "Input"
  data_list <- append(data_list, list(data))
}

if (!is.null(opt$input_A)) {
  data_A <- read.table(opt$input_A, header = FALSE, sep = "\t", col.names = c("ID", "Length"))
  data_A$Group <- "Input_A"
  data_list <- append(data_list, list(data_A))
}

if (!is.null(opt$input_B)) {
  data_B <- read.table(opt$input_B, header = FALSE, sep = "\t", col.names = c("ID", "Length"))
  data_B$Group <- "Input_B"
  data_list <- append(data_list, list(data_B))
}

# Merge all data
combined_data <- bind_rows(data_list)

# Set ID column as a factor and sort IDs naturally
combined_data$ID <- factor(combined_data$ID, levels = unique(combined_data$ID))

# Convert Length to Mb
combined_data$Length <- combined_data$Length / 1e6

# Create violin plot
p <- ggplot(combined_data, aes(x = ID, y = Length, fill = ID)) +
  geom_violin(trim = FALSE, alpha = 0.8) +  # Adjust transparency of violin plot
  
  # Input data points (small, transparent, same color as violin plot)
  geom_jitter(data = subset(combined_data, Group == "Input"), 
              aes(color = Group), width = 0.2, size = 0.8, alpha = 0.3) +
  
  # Input_A data points (large and dark red)
  geom_jitter(data = subset(combined_data, Group == "Input_A"), 
              aes(color = Group), width = 0.2, size = 2, alpha = 0.7) +
  
  # Input_B data points (large and dark blue)
  geom_jitter(data = subset(combined_data, Group == "Input_B"), 
              aes(color = Group), width = 0.2, size = 2, alpha = 0.7) +
  
  # Custom colors and legend labels
  scale_color_manual(
    values = c("Input_A" = "darkred", "Input_B" = "darkblue", "Input" = "gray50"),
    labels = c("Input_A" = "cq_mat", "Input_B" = "cq_pat", "Input" = "102 samples"),
    guide = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  
  # Mean line: thinner and shorter
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.5) +
  
  # Set Y-axis scale with 1-unit increments
  scale_y_continuous(breaks = seq(0, max(combined_data$Length, na.rm = TRUE) + 1, by = 1)) +
  
  # Legend and axis settings
  labs(title = "",
       x = "Chromosome",
       y = "Length of HOR arrays (Mb)",
       color = "Group") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),  # X-axis text size 15
    axis.text.y = element_text(size = 15),  # Y-axis text size 15
    axis.title.x = element_text(size = 15),  # X-axis label size 15
    axis.title.y = element_text(size = 15),
    axis.line = element_line(color = "black"),  # Solid axis lines
    axis.line.y.right = element_blank(),       # Remove right border
    axis.line.x.top = element_blank(),         # Remove top border
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 15),
    legend.key = element_blank(),              # Remove legend background
    panel.border = element_blank(),            # Remove panel border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  # Remove ID-related legend
  guides(
    fill = "none",
    color = guide_legend(order = 1)
  )

# Save as PDF
output_file <- "violin_plot_combined.pdf"
ggsave(output_file, plot = p, width = 12, height = 6, dpi = 300, device = "pdf")

# Display the plot
print(p)

cat(sprintf("Violin plot saved to %s\n", output_file))
