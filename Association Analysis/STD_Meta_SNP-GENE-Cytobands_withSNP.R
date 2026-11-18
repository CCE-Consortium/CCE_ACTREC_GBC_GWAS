# Title: Gallbladder Cancer GWAS Meta-Analysis
# Analysis: Manhattan Plot and Quantile-Quantile (QQ) Plot
# Author: Shubham G Chavan

#
# Required R package:
#   data.table,qqman


# Load required packages
library(qqman)
library(data.table)

# Define the input GWAS meta-analysis file
filePath <- paste0("Path of CSV main fine")

# Define the folder for saving final plots
OutputPath <- paste0("Path of Output results")

# Create the output folder if it does not exist
if (!dir.exists(OutputPath)) dir.create(OutputPath, recursive = TRUE)

# Define chromosome colors for the Manhattan plot
colors_new <- c(
  "#1874CD99",
  "#79CDCD99",
  "#00BFFF75"
)

# Define the maximum Y-axis limit
YLim <- 20

# Define Lead SNP rsIDs
SNPsOfInterest_RSID <- c(
  "rs6733452",
  "rs15734",
  "rs67274221",
  "rs35604230",
  "rs60541816",
  "rs77155094",
  "rs2312959",
  "rs2131242"
)

# Define gene names for Lead SNPs
gene_names <- list(
  rs6733452 = "ABCG8",
  rs15734 = "LMCD1-AS1",
  rs67274221 = "ARHGEF3",
  rs35604230 = "LSM2",
  rs60541816 = "ABCB1",
  rs77155094 = "PIGCP2",
  rs2312959 = "ABAT",
  rs2131242 = "ANO1"
)

# Define cytoband locations for Lead SNPs
cytoband_names <- list(
  rs6733452 = "2p21",
  rs15734 = "3p26.1",
  rs67274221 = "3p14.3",
  rs35604230 = "6p21.33",
  rs60541816 = "7q21.12",
  rs77155094 = "7q31.1",
  rs2312959 = "16p13.2",
  rs2131242 = "11q13.3"
)

# Read, clean, validate, and select the strongest record for each SNP
clean_data <- function(filePath) {
  
  # Read only the required GWAS columns
  df <- fread(
    file = filePath,
    select = c("CHR", "BP", "P", "SNP")
  )
  
  # Convert chromosome values to numeric
  df[, CHR := suppressWarnings(as.numeric(CHR))]
  
  # Convert base-pair positions to numeric
  df[, BP := suppressWarnings(as.numeric(BP))]
  
  # Convert P-values to numeric
  df[, P := suppressWarnings(as.numeric(P))]
  
  # Convert SNP IDs to character and remove leading or trailing spaces
  df[, SNP := trimws(as.character(SNP))]
  
  # Keep only valid chromosome, position, SNP, and P-value records
  df <- df[
    !is.na(CHR) &
      !is.na(BP) &
      !is.na(P) &
      !is.na(SNP) &
      SNP != "" &
      P > 0 &
      P <= 1
  ]
  
  # For duplicated SNPs, retain the record with the smallest P-value
  df <- df[
    order(P),
    .SD[1],
    by = SNP
  ]
  
  # Remove the grouping column created by data.table
  df <- df[, .(
    SNP,
    CHR,
    BP,
    P
  )]
  
  # Sort the cleaned dataset by chromosome and base-pair position
  setorder(df, CHR, BP)
  
  # Return the cleaned dataset
  return(df)
}

# Calculate chromosome lengths and chromosome center positions
calculate_x_positions <- function(df) {
  
  # Calculate the maximum base-pair position for each chromosome
  chr_lengths <- tapply(df$BP, df$CHR, max)
  
  # Sort chromosomes numerically
  chr_lengths <- chr_lengths[
    order(as.numeric(names(chr_lengths)))
  ]
  
  # Store chromosome names
  chr_names <- names(chr_lengths)
  
  # Calculate the starting position of each chromosome
  chr_start <- c(
    0,
    cumsum(chr_lengths)[-length(chr_lengths)]
  )
  
  # Calculate the center position of each chromosome
  chr_center <- chr_start + chr_lengths / 2
  
  # Assign chromosome names to starting positions
  names(chr_start) <- chr_names
  
  # Assign chromosome names to center positions
  names(chr_center) <- chr_names
  
  # Return chromosome position information
  return(list(
    lengths = chr_lengths,
    start = chr_start,
    center = chr_center
  ))
}

# Calculate the Manhattan X-axis position
calculate_x_pos <- function(chr, bp, chr_positions) {
  
  # Add chromosome start position to the base-pair position
  return(
    chr_positions$start[as.character(chr)] + bp
  )
}

# Create the publication-quality Manhattan plot
plot_manhattan_labels <- function(
    df,
    colors,
    snp_list_rsid,
    gene_names,
    cytoband_names,
    outFile,
    ylim) {
  
  # Open high-resolution PNG graphics device
  png(
    filename = outFile,
    width = 6000,
    height = 2800,
    res = 300,
    type = "cairo-png"
  )
  
  # Calculate chromosome positions
  chr_positions <- calculate_x_positions(df)
  
  # Set plot margins and text sizes
  par(
    mar = c(8, 5, 5, 2),
    xaxs = "i",
    cex.lab = 1.8,
    cex.axis = 1.8
  )
  
  # Draw the Manhattan plot
  manhattan(
    df,
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    main = "",
    ylim = c(0, ylim),
    xlab = "Chromosomes",
    col = colors,
    annotateTop = FALSE,
    genomewideline = -log10(5e-8),
    suggestiveline = FALSE,
    cex = 1.0,
    xaxt = "n"
  )
  
  # Add Manhattan plot title at the upper-left
  mtext(
    text = "Manhattan Plot of Gallbladder Cancer GWAS Meta-Analysis (Initial Cohort + Expanded Cohort)",
    side = 3,
    line = 2.2,
    adj = 0.02,
    cex = 1.6,
    font = 2
  )
  
  # Add manually calculated chromosome labels
  axis(
    side = 1,
    at = chr_positions$center,
    labels = names(chr_positions$center),
    cex.axis = 1.8,
    las = 1
  )
  
  # Add Y-axis labels every five units
  axis(
    side = 2,
    at = seq(0, ylim, 5),
    las = 1,
    cex.axis = 1.8
  )
  
  # Define vertical distance above the Lead SNP
  label_up_shift <- 1.5
  
  # Define distance between SNP and gene labels
  snp_to_gene_gap <- 0.65
  
  # Define distance between gene and cytoband labels
  gene_to_cytoband_gap <- 0.65
  
  # Define Lead SNP label size
  label_cex <- 1.5
  
  # Define Lead SNP label color
  label_colour <- "black"
  
  # Define bold font for SNP and cytoband labels
  label_font <- 2
  
  # Loop through all Lead SNPs
  for (i in seq_along(snp_list_rsid)) {
    
    # Select the current Lead SNP
    rsid <- snp_list_rsid[i]
    
    # Retrieve the SNP directly from the cleaned dataset
    snp_row <- df[SNP == rsid]
    
    # Skip the SNP if it is not found
    if (nrow(snp_row) == 0) next
    
    # Use the first retained record for the SNP
    snp_row <- snp_row[1]
    
    # Extract chromosome from the dataset
    this_chr <- snp_row$CHR[1]
    
    # Extract base-pair position from the dataset
    this_bp <- snp_row$BP[1]
    
    # Extract the P-value from the dataset
    this_p <- snp_row$P[1]
    
    # Calculate the Manhattan X-axis position
    x_pos <- calculate_x_pos(
      this_chr,
      this_bp,
      chr_positions
    )
    
    # Calculate negative log10 P-value
    y_pos <- -log10(this_p)
    
    # Highlight the Lead SNP with green fill and black border
    points(
      x = x_pos,
      y = y_pos,
      pch = 21,
      bg = "green3",
      col = "black",
      cex = 0.85,
      lwd = 1.5
    )
    
    # Calculate vertical position of the SNP label
    snp_y <- y_pos + label_up_shift
    
    # Calculate vertical position of the gene label
    gene_y <- snp_y - snp_to_gene_gap
    
    # Calculate vertical position of the cytoband label
    cytoband_y <- gene_y - gene_to_cytoband_gap
    
    # Set default X-axis position for labels
    x_label <- x_pos
    
    # Set default SNP label position
    snp_y_label <- snp_y
    
    # Set default gene label position
    gene_y_label <- gene_y
    
    # Set default cytoband label position
    cytoband_y_label <- cytoband_y
    
    # Move rs15734 slightly left and downward
    if (rsid == "rs15734") {
      x_label <- x_pos - 0.02 * diff(par("usr")[1:2])
      snp_y_label <- snp_y - 0.25
      gene_y_label <- gene_y - 0.25
      cytoband_y_label <- cytoband_y - 0.25
    }
    
    # Move rs67274221 slightly to the right
    if (rsid == "rs67274221") {
      x_label <- x_pos + 0.02 * diff(par("usr")[1:2])
    }
    
    # Add the Lead SNP rsID
    text(
      x = x_label,
      y = snp_y_label,
      labels = rsid,
      pos = 3,
      cex = label_cex,
      font = label_font,
      col = label_colour
    )
    
    # Retrieve the corresponding gene
    gene <- gene_names[[rsid]]
    
    # Add the gene label
    if (!is.null(gene)) {
      text(
        x = x_label,
        y = gene_y_label,
        labels = gene,
        pos = 3,
        cex = label_cex,
        font = 4,
        col = label_colour
      )
      
      # Retrieve the corresponding cytoband
      cytoband <- cytoband_names[[rsid]]
      
      # Add the cytoband label
      if (!is.null(cytoband)) {
        text(
          x = x_label,
          y = cytoband_y_label,
          labels = cytoband,
          pos = 3,
          cex = label_cex,
          font = label_font,
          col = label_colour
        )
      }
    }
  }
  
  # Close the PNG graphics device
  dev.off()
}

# Read and clean the GWAS meta-analysis dataset
df <- clean_data(filePath)

# Check the exact rs2131242 result used in the Manhattan plot
print(
  df[
    SNP == "rs2131242",
    .(
      SNP,
      CHR,
      BP,
      P,
      NegLog10P = -log10(P)
    )
  ]
)

# Generate the Manhattan plot
plot_manhattan_labels(
  df = df,
  colors = colors_new,
  snp_list_rsid = SNPsOfInterest_RSID,
  gene_names = gene_names,
  cytoband_names = cytoband_names,
  outFile = paste0(
    OutputPath,
    "/Manhattan_RSID_GENE_Cytobands.png"
  ),
  ylim = YLim
)



# ============================================================
# SECTION 2: QQ PLOT
# ============================================================

# Create the publication-quality QQ plot
plot_qq <- function(df, outFile) {
  
  # Extract valid P-values
  p_values <- df$P[
    !is.na(df$P) &
      df$P > 0 &
      df$P <= 1
  ]
  
  # Convert P-values to chi-square statistics
  chi_square <- qchisq(
    1 - p_values,
    df = 1
  )
  
  # Calculate genomic inflation factor Lambda
  lambda_value <- median(
    chi_square,
    na.rm = TRUE
  ) / qchisq(
    0.5,
    df = 1
  )
  
  # Open high-resolution PNG graphics device
  png(
    filename = outFile,
    width = 3000,
    height = 3000,
    res = 300,
    type = "cairo-png"
  )
  
  # Set QQ plot margins
  par(
    mar = c(7, 7, 5, 3)
  )
  
  # Calculate number of P-values
  n <- length(p_values)
  
  # Calculate expected -log10 P-values
  expected <- -log10(
    (seq_len(n) - 0.5) / n
  )
  
  # Calculate observed -log10 P-values
  observed <- -log10(
    sort(p_values)
  )
  
  # Draw the QQ plot
  plot(
    expected,
    observed,
    pch = 20,
    col = "#701705",
    cex = 1.0,
    xlab = "Expected -log10(p)",
    ylab = "Observed -log10(p)",
    main = "",
    font.lab = 1,
    font.axis = 1,
    cex.lab = 1.8,
    cex.axis = 1.6
  )
  
  # Add the QQ plot title
  mtext(
    text = "Q–Q Plot of Gallbladder Cancer GWAS Meta-Analysis (Initial Cohort + Expanded Cohort)",
    side = 3,
    line = 1,
    adj = 0.70,
    cex = 1.4,
    font = 2
  )
  
  # Add the theoretical null diagonal
  abline(
    a = 0,
    b = 1,
    lwd = 2,
    col = "#FF0000"
  )
  
  # Format Lambda to four decimal places
  lambda_text <- paste0(
    "\u03bb = ",
    format(
      round(lambda_value, 4),
      nsmall = 4
    )
  )
  
  # Add Lambda value to the plot
  text(
    x = min(expected) + 0.05 * diff(range(expected)),
    y = max(observed) - 0.05 * diff(range(observed)),
    labels = lambda_text,
    adj = c(0, 1),
    cex = 1.8,
    font = 2,
    col = "#2B0C05"
  )
  
  # Close the PNG graphics device
  dev.off()
}

# Generate the QQ plot
plot_qq(
  df = df,
  outFile = paste0(
    OutputPath,
    "/QQ_Plot.png"
  )
)