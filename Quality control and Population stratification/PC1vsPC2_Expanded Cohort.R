
# PCA PLOTS: PC1 vs PC2 THROUGH PC9 vs PC10 - GBC PHASE 2

options(bitmapType = "cairo")

library(ggplot2)
library(grid)
library(readxl)

# SECTION 1: READ PCA DATA

input_file <- "D:/Shubham G Chavan/SSD change new work/GBC_ARTICLE/PCA PLOT/PCA-Phase2/Main File/P2_COVAR_FILE_7101Samples_PCA-maf-0.05_regular_10PCs_Age_Gen AM SC -Removed-9.xlsx"

pca <- as.data.frame(read_excel(input_file))

initial_samples <- nrow(pca)

# SECTION 2: CHECK REQUIRED COLUMNS

pc_cols <- paste0("PC", 1:10)

required_cols <- c(
  pc_cols,
  "Birth_res_Zone_genetic",
  "Case_control"
)

missing_cols <- setdiff(required_cols, colnames(pca))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

cat("\nInitial total samples:", initial_samples, "\n")

# SECTION 3: CHECK AND FILTER PCA VALUES

# Identify samples containing -9 in any PC column

minus9_rows <- apply(
  pca[, pc_cols, drop = FALSE],
  1,
  function(x) any(x == -9, na.rm = TRUE)
)

cat("\nSamples containing -9:", sum(minus9_rows), "\n")

# Remove confirmed -9 missing-value codes

pca <- pca[!minus9_rows, , drop = FALSE]

# Remove missing and non-finite PCA values

valid_rows <- complete.cases(pca[, pc_cols, drop = FALSE])

for (pc in pc_cols) {
  valid_rows <- valid_rows & is.finite(pca[[pc]])
}

additional_removed <- sum(!valid_rows)

pca <- pca[valid_rows, , drop = FALSE]

cat("\nAdditional samples removed due to missing/non-finite PC values:",
    additional_removed, "\n")

cat("\nFinal samples available for PCA plots:", nrow(pca), "\n")

if (nrow(pca) == 0) {
  stop("No valid PCA samples remain after filtering.")
}

# SECTION 4: DEFINE PCA COMBINATIONS

plot_pairs <- list(
  c("PC1", "PC2")
  , c("PC2", "PC3"), c("PC3", "PC4"),c("PC4", "PC5"),  c("PC5", "PC6"),c("PC6", "PC7"),  c("PC7", "PC8"),
  c("PC8", "PC9"),  c("PC9", "PC10"))
# SECTION 5: FORMAT GROUPING VARIABLES

pca$Case_control <- factor(
  pca$Case_control,
  levels = c("Case", "Control")
)

pca$Birth_res_Zone_genetic <- factor(
  pca$Birth_res_Zone_genetic,
  levels = sort(unique(na.omit(pca$Birth_res_Zone_genetic)))
)

# SECTION 6: DEFINE COLOURS AUTOMATICALLY

# Original nine colours retained


colours_regions <- c(
  "#F2F200", "#81C784", "#FF0000",
  "#7F00FF", "#1837B5", "#08E4FF",
  "#FF7D47", "#770737", "#04EB00"
)


# Identify all available State Code regions

region_levels <- levels(
  droplevels(pca$Birth_res_Zone_genetic)
)

n_regions <- length(region_levels)

# Automatically generate additional colours if required

if (n_regions > length(colours_regions)) {
  
  additional_colours <- grDevices::hcl.colors(
    n_regions - length(colours_regions),
    palette = "Dynamic"
  )
  
  colours_regions <- c(
    colours_regions,
    additional_colours
  )
}

# Assign a unique colour to each region

region_colours <- setNames(
  colours_regions[seq_len(n_regions)],
  region_levels
)

cat("\nTotal State Code regions:", n_regions, "\n")
cat("Total colours available:", length(region_colours), "\n")

cat("\nRegion-wise colour allocation:\n")
print(region_colours)

# SECTION 7: DEFINE OUTPUT FOLDER

out_dir <- "D:/Shubham G Chavan/SSD change new work/GBC_ARTICLE/PCA PLOT/PCA-Phase2/Graph&Plot/Final_Plots_Clr"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# SECTION 8: GENERATE AND SAVE PCA PLOTS

for (pair in plot_pairs) {
  
  x_pc <- pair[1]
  y_pc <- pair[2]
  
  p <- ggplot(
    pca,
    aes(
      x = .data[[x_pc]],
      y = .data[[y_pc]],
      colour = Birth_res_Zone_genetic,
      shape = Case_control
    )
  ) +
    
    geom_point(
      size = 1.8,
      alpha = 0.85,
      stroke = 0.6,
      na.rm = TRUE
    ) +
    
    scale_shape_manual(
      values = c(
        "Case" = 17,
        "Control" = 16
      ),
      name = "Case & Control",
      drop = TRUE
    ) +
    
    scale_colour_manual(
      values = region_colours,
      name = "State Code",
      drop = TRUE
    ) +
    
    labs(
      title = paste("GBC PCA PLOT:", x_pc, "vs", y_pc),
      x = x_pc,
      y = y_pc
    ) +
    
    theme_classic() +
    
    theme(
      axis.text = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.spacing.x = unit(0.4, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      legend.background = element_rect(
        fill = "transparent",
        colour = "transparent"
      ),
      legend.key = element_rect(
        fill = "transparent",
        colour = "transparent"
      ),
      
      plot.margin = margin(
        t = 5, r = 5, b = 5, l = 5,
        unit = "pt"
      ),
      
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      panel.grid.minor = element_blank()
    )
  
  # SAVE EACH GRAPH WITH ITS CORRECT PC COMBINATION NAME
  
  outfile <- file.path(
    out_dir,
    paste0(x_pc, "_vs_", y_pc, ".png")
  )
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  cat("\nPCA plot successfully saved:", outfile, "\n")
}

# SECTION 9: COMPLETION SUMMARY

cat("\n========================================")
cat("\nGBC PHASE 2 PCA PLOT SUMMARY")
cat("\n========================================")

cat("\nInitial samples:", initial_samples)
cat("\nSamples containing -9:", sum(minus9_rows))
cat("\nAdditional invalid samples removed:", additional_removed)
cat("\nFinal samples used:", nrow(pca))
cat("\nTotal State Code regions:", n_regions)
cat("\nTotal PCA plots generated:", length(plot_pairs))
cat("\nOutput folder:", out_dir)

cat("\n\nAll 9 PCA plots successfully generated!\n")
