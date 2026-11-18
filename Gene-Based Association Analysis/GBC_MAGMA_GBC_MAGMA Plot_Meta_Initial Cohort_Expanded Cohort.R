
# MAGMA GENE-BASED MANHATTAN PLOT-GBC_MAGMA_Gene_Meta(Initial Cohort + Expanded Cohort)
# Author: Shubham G Chavan


# SECTION 1: LOAD REQUIRED PACKAGES

library(data.table)
library(ggplot2)


# SECTION 2: DEFINE INPUT AND OUTPUT PATHS

# Input MAGMA gene-based results file
input_file <- "D:/Shubham G Chavan/SSD change new work/GBC_ARTICLE/GBC-MAGMA-Files_All-3/GBC-Meta/Main File/magma.genes.out.csv"

# Output directory for Manhattan plot
output_dir <- "D:/Shubham G Chavan/SSD change new work/GBC_ARTICLE/GBC-MAGMA-Files_All-3/GBC-Meta/Graphs & Plot"

# Output PNG filename
output_file <- file.path(
  output_dir,
  "GBC_MAGMA_Gene_Meta(Initial Cohort + Expanded Cohort).png"
)

# Create output directory if it does not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# SECTION 3: READ AND CLEAN MAGMA DATA

clean_magma_data <- function(filePath) {
  
  # Read MAGMA output file
  df <- fread(filePath)
  
  # Convert column names to uppercase
  setnames(df, names(df), toupper(names(df)))
  
  # Check required SYMBOL column
  if (!"SYMBOL" %in% names(df)) {
    stop("ERROR: SYMBOL column not found in MAGMA output.")
  }
  
  # Define required columns
  required_cols <- c("CHR", "START", "STOP", "P")
  
  # Check missing columns
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Convert required columns to numeric
  df[, CHR := as.numeric(CHR)]
  df[, START := as.numeric(START)]
  df[, STOP := as.numeric(STOP)]
  df[, P := as.numeric(P)]
  
  # Retain autosomal chromosomes
  df <- df[CHR >= 1 & CHR <= 22]
  
  # Calculate gene midpoint
  df[, BP := (START + STOP) / 2]
  
  # Remove missing and invalid P-values
  df <- df[!is.na(P) & P > 0]
  
  return(df)
}


# SECTION 4: CALCULATE GENOMIC POSITIONS

calculate_gene_positions <- function(df) {
  
  # Calculate maximum gene position per chromosome
  chr_len <- df[, .(maxBP = max(BP)), by = CHR]
  
  # Arrange chromosomes numerically
  setorder(chr_len, CHR)
  
  # Calculate cumulative genomic positions
  chr_len[, cum := cumsum(c(0, head(maxBP, -1)))]
  
  # Merge chromosome information with gene data
  df <- merge(df, chr_len, by = "CHR")
  
  # Calculate continuous genome-wide positions
  df[, pos := BP + cum]
  
  # Calculate chromosome midpoint positions
  axis_df <- chr_len[, .(CHR, center = cum + maxBP / 2)]
  
  return(list(data = df, axis = axis_df))
}


# SECTION 5: GENERATE MAGMA GENE-BASED MANHATTAN PLOT

plot_magma_manhattan <- function(df, axis_df, title, outFile) {
  
  # Define genome-wide significance threshold
  gw_p <- 2.5e-6
  
  # Calculate -log10(P)
  df[, LOGP := -log10(P)]
  
  # Identify significant genes
  df[, SIGNIFICANT := P <= gw_p]
  
  # Separate significant and non-significant genes
  bg_df <- df[SIGNIFICANT == FALSE]
  sig_df <- df[SIGNIFICANT == TRUE]
  
  # Define maximum Y-axis value
  ymax <- ceiling(max(df$LOGP)) + 1
  
  # Initialize Manhattan plot
  p <- ggplot() +
    
    # Plot non-significant genes
    geom_point(
      data = bg_df,
      aes(x = pos, y = LOGP, color = as.factor(CHR %% 2)),
      size = 2,
      alpha = 0.85
    ) +
    
    # Highlight significant genes in green
    geom_point(
      data = sig_df,
      aes(x = pos, y = LOGP),
      shape = 21,
      fill = "green",
      color = "black",
      size = 2,
      stroke = 0.5
    ) +
    
    # Add genome-wide significance threshold
    geom_hline(
      yintercept = -log10(gw_p),
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    )
  
  
  # Select genes above significance threshold for labelling
  label_df <- df[LOGP > -log10(gw_p) & !is.na(SYMBOL)]
  
  if (nrow(label_df) > 0) {
    
    # Add default gene labels
    p <- p +
      geom_text(
        data = label_df[!SYMBOL %in% c("HSPA1L", "HLA-DQA1")],
        aes(x = pos, y = LOGP, label = SYMBOL),
        vjust = -1,
        hjust = 0.5,
        size = 3.4,
        fontface = "bold",
        color = "black",
        check_overlap = TRUE
      )
    
    # Manually position HSPA1L label
    p <- p +
      geom_text(
        data = label_df[SYMBOL == "HSPA1L"],
        aes(x = pos, y = LOGP, label = SYMBOL),
        vjust = 1.2,
        hjust = 1.2,
        size = 3.4,
        fontface = "bold",
        color = "black"
      )
    
    # Manually position HLA-DQA1 label
    p <- p +
      geom_text(
        data = label_df[SYMBOL == "HLA-DQA1"],
        aes(x = pos, y = LOGP, label = SYMBOL),
        vjust = -0.3,
        hjust = -0.2,
        size = 3.4,
        fontface = "bold",
        color = "black"
      )
  }
  
  
  # Add chromosome colours, axes and publication formatting
  p <- p +
    
    # Define alternating chromosome colours
    scale_color_manual(
      values = c("#420338", "#875d80")
    ) +
    
    # Display chromosome numbers
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$CHR
    ) +
    
    # Define Y-axis limits and intervals
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = seq(0, ymax, by = 2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    # Define plot title and axis labels
    labs(
      title = title,
      x = "Chromosome",
      y = expression(-log[10](p))
    ) +
    
    # Apply black-and-white theme
    theme_bw(base_size = 15) +
    
    # Customize plot appearance
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.border = element_blank(),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 16
      ),
      
      axis.text = element_text(
        size = 14,
        color = "black"
      )
    )
  
  
  # Save final Manhattan plot as high-resolution PNG
  ggsave(
    filename = outFile,
    plot = p,
    width = 15,
    height = 6.5,
    dpi = 300
  )
}


# SECTION 6: RUN MAGMA ANALYSIS AND GENERATE PLOT

# Read and clean MAGMA gene-based results
df <- clean_magma_data(input_file)

# Calculate continuous genomic positions
pos <- calculate_gene_positions(df)

# Generate and save final Manhattan plot
plot_magma_manhattan(
  pos$data,
  pos$axis,
  "MAGMA Gene-Based Manhattan Plot of Gallbladder Cancer-Meta-analysis (Initial Cohort + Expanded Cohort)",
  output_file
)


# SECTION 7: DISPLAY COMPLETION MESSAGE

cat("\nMAGMA Gene-Based Manhattan Plot generated successfully!\n")
cat("Input file: ", input_file, "\n")
cat("Output file: ", output_file, "\n")
