# ================================================================
# Draw CCA Variable Circle Plots (from pre-computed coordinates)
# ================================================================
# - Input dir: 'cca_coords_exports/' with files named:
#              *_features_labeled_cutNN.csv (Block, Code, x, y, radius, ...)
# - Output dir: 'circle_from_coords_plots/' with one PDF per input file
# - Dependencies: ggplot2, ggrepel, tools
# ================================================================

library(ggplot2)
library(ggrepel)
library(tools)

## ---------------- Config ----------------
INPUT_DIR  <- "cca_coords_exports"
OUT_DIR    <- "circle_from_coords_plots"
LABEL_SIZE <- 2.8   # label text size
POINT_SIZE <- 1.0   # point size for features
INNER_R    <- 0.50  # optional inner cutoff ring

# Maximum number of labels per block (top by radius)
N_PER_BLOCK <- c(
  "MetaT"  = 50,
  "RNAseq" = 50,
  "Bile"   = 50,
  "16S"    = 50,
  "MetaG"  = 50
)

# Colors for each omics block
COL_BLOCK <- c(
  "Bile"   = "#6c7faa",
  "16S"    = "#7fa634",
  "MetaT"  = "#e25f32",
  "MetaG"  = "#c3579c",
  "RNAseq" = "#3fa181"
)

dir.create(OUT_DIR, showWarnings = FALSE)

# Helper to generate a circle for plotting
circle_df <- function(r = 1, n = 1000) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = r * cos(t), y = r * sin(t))
}

## --------------- Per-file plotter ---------------
make_plot <- function(csv_path) {
  df <- tryCatch(
    read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) {
    message("Could not read: ", csv_path)
    return(invisible(NULL))
  }

  # Require these columns
  needed <- c("Block", "Code", "x", "y")
  if (!all(needed %in% names(df))) {
    message("Skipping (missing columns): ", csv_path)
    return(invisible(NULL))
  }

  # Plot title from filename, e.g. "MetaT_VS_Bile_features_labeled_cut70.csv"
  base <- basename(csv_path)
  title_txt <- base
  title_txt <- sub("_features_labeled_cut\\d+\\.csv$", "", title_txt)
  title_txt <- gsub("_", " ", title_txt)
  title_txt <- gsub(" VS ", " vs ", title_txt, fixed = TRUE)
  title_txt <- paste("CCA Variable Circle:", title_txt)

  # Restrict palette to blocks present
  blocks_present <- intersect(unique(df$Block), names(COL_BLOCK))
  pal <- COL_BLOCK[blocks_present]

  # --- Select top-N features per block (by radius) ---
  df$radius <- sqrt(df$x^2 + df$y^2)
  df$to_label <- FALSE
  for (blk in unique(df$Block)) {
    idx <- which(df$Block == blk)
    n   <- if (!is.na(N_PER_BLOCK[blk])) N_PER_BLOCK[blk] else 50
    ord <- idx[order(df$radius[idx], decreasing = TRUE)]
    keep <- head(ord, min(n, length(ord)))
    df$to_label[keep] <- TRUE
  }
  df_lab <- subset(df, to_label)

  # Circles
  circ_outer <- circle_df(1.0)
  circ_inner <- if (!is.na(INNER_R) && INNER_R > 0) circle_df(INNER_R) else NULL

  # ---- Plot ----
  p <- ggplot() +
    geom_path(data = circ_outer, aes(x, y), color = "grey60", linewidth = 0.6) +
    { if (!is.null(circ_inner))
        geom_path(data = circ_inner, aes(x, y),
                  color = "grey80", linewidth = 0.4, linetype = 2)
      else NULL } +
    geom_hline(yintercept = 0, color = "grey85", linetype = 2) +
    geom_vline(xintercept = 0, color = "grey85", linetype = 2) +
    geom_point(data = df, aes(x, y, color = Block), size = POINT_SIZE) +
    geom_text_repel(
      data = df_lab,
      aes(x, y, label = Code, color = Block),
      size = LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.6,
      point.padding = 0.5,
      force = 2.5,
      min.segment.length = 0,
      segment.size = 0.2,
      segment.color = "grey40",
      segment.alpha = 0.9,
      seed = 42
    ) +
    coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(title = title_txt, x = "Comp 1", y = "Comp 2", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")

  out_pdf <- file.path(OUT_DIR, paste0(file_path_sans_ext(base), "_circle_from_coords.pdf"))
  ggsave(out_pdf, p, width = 7, height = 6)
  message("Saved: ", out_pdf)
}

## --------------- Run all ---------------
if (!dir.exists(INPUT_DIR)) {
  message("Input dir not found: ", INPUT_DIR)
} else {
  files <- list.files(INPUT_DIR, pattern = "_features_labeled_cut\\d+\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("No *_features_labeled_cutNN.csv files found in ", INPUT_DIR)
  } else {
    invisible(lapply(files, make_plot))
  }
}