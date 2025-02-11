#!/usr/bin/env Rscript
# This script is revised from https://github.com/logsdon-lab/CenMAP/blob/main/workflow/scripts/plot_HOR_length.R

# -------------------------
# 📚 Load Required Libraries
# -------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(argparser)
library(cowplot)

# -------------------------
# 📌 Argument Parsing
# -------------------------
p <- arg_parser("Plot STV from formatted HumAS-HMMER output for 3 datasets.")
p <- add_argument(p, "--input", help = "First input HumAS-HMMER formatted output.", type = "character")
p <- add_argument(p, "--input1", help = "Second input HumAS-HMMER formatted output.", type = "character")
p <- add_argument(p, "--input2", help = "Third input HumAS-HMMER formatted output.", type = "character")
p <- add_argument(p, "--chr", help = "Chromosome to plot. ex. chrX", type = "character")
p <- add_argument(p, "--output", help = "Output comparison plot file. Defaults to {chr}_comparison.png.", type = "character", default = NA)
p <- add_argument(p, "--hor_filter", help = "Filter for HORs that occur at least 20 times.", type = "numeric", default = 0)
p <- add_argument(p, "--mer_order", help = "Reorder data to put 'large'-r or 'small'-r -mers on top.", type = "character", default = "large")
p <- add_argument(p, "--plot_width", help = "Plot width", default = 10)
p <- add_argument(p, "--plot_height", help = "Plot height per plot", default = 5)
args <- parse_args(p)

# Set default output filename
args$output <- ifelse(
  is.na(args$output),
  paste0(args$chr, "_comparison.png"),
  args$output
)

# -------------------------
# 📊 Data Processing Function
# -------------------------
process_input <- function(input_file, chr, hor_filter, mer_order) {
  cols_to_take <- seq(5)
  monomer_len <- 170
  cols <- c("chr", "start", "stop", "hor", "strand")
  
  samples <- fread(input_file,
    sep = "\t",
    stringsAsFactors = TRUE,
    fill = TRUE, quote = "",
    header = FALSE, select = cols_to_take,
    col.names = cols
  )
  
  samples$length <- samples$stop - samples$start
  samples$mer <- as.numeric(round(samples$length / monomer_len))
  
  samples <- switch(chr,
    "chr10" = subset(samples, as.numeric(mer) >= 5),
    "chr20" = subset(samples, as.numeric(mer) >= 5),
    "chrY" = subset(samples, as.numeric(mer) >= 30),
    "chr17" = subset(samples, as.numeric(mer) >= 4),
    samples
  )
  
  df_mer <- samples %>%
    group_by(mer) %>%
    filter(n() > hor_filter)
  
  df_final <- df_mer %>%
    group_by(chr) %>%
    mutate(stop = stop - min(start)) %>%
    mutate(start = start - min(start)) %>%
    mutate(chr_order = sprintf("%02d", as.numeric(str_extract(chr, "(?<=\\d).*?(?=\\))")))) %>%
    mutate(chr_order = factor(chr_order, levels = sort(unique(chr_order), decreasing = FALSE))) %>%
    arrange(chr_order) %>%
    mutate(chr = factor(chr, levels = unique(chr)))
  
  df_final <- switch(mer_order,
    "large" = df_final %>% arrange(mer),
    "small" = df_final %>% arrange(-mer),
    stop(paste("Invalid mer reordering option:", mer_order))
  )
  
  # -------------------------
  # 🎨 Fixed Color Mapping
  # -------------------------
  myColors <- c(
    "#D66C54", "#C6625D", "#F4DC78", "#7EC0B3", "#A8275C", "#8CC49F", "#893F89", "#6565AA", "#9AC78A",
    "#cc8fc1", "#3997C6", "#8882c4", "#8abdd6", "#38A49B", "#45B4CE", "#afa7d8", "#a874b5", "#3F66A0", "#BFDD97",
    "#af5d87", "#E5E57A", "#ed975d", "#F9E193", "#E5D1A1", "#A1B5E5", "#9F68A5", "#81B25B"
  )

  names(myColors) <- c(
    "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
    "21", "22", "24", "26", "30", "32", "34", "35"
  )
  
  # -------------------------
  # 📊 Plot Generation
  # -------------------------
  p <- ggplot(df_final) +
    geom_segment(aes(x = start / 1e6, xend = stop / 1e6 + 0.002, y = chr, yend = chr, color = as.factor(mer)), linewidth = 9) +
    scale_color_manual(values = myColors) +
    theme_classic(base_size = 15) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 15)
    ) +
    xlab("Position (Mb)") +
    ylab("HOR size") +
    guides(color = guide_legend(override.aes = list(size = 1)))  # Adjust legend icon size to 1
  return(p)
}

# -------------------------
# 📊 Generate Three Plots
# -------------------------
plot1 <- process_input(args$input, args$chr, args$hor_filter, args$mer_order)
plot2 <- process_input(args$input1, args$chr, args$hor_filter, args$mer_order) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())
plot3 <- process_input(args$input2, args$chr, args$hor_filter, args$mer_order) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

# -------------------------
# 📊 Combine Plots and Adjust Legend
# -------------------------
legend <- get_legend(plot1 + theme(legend.position = "right"))

final_plot <- plot_grid(
  plot1 + theme(legend.position = "none"),
  plot2 + theme(legend.position = "none"),
  plot3 + theme(legend.position = "none"),
  legend,
  ncol = 4,
  rel_widths = c(1, 1, 1, 0.3)
)

# -------------------------
# 💾 Save Plot
# -------------------------
ggsave(args$output, plot = final_plot, width = args$plot_width * 3.5, height = args$plot_height, limitsize = FALSE)
