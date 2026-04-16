# ======================================================================
# rCCA on full features + export circle coordinates
# - Fits rCCA for every dataset pair
# - Computes variable circle coordinates (correlations w/ canonical variates)
# - Exports:
#     * *_features_all.csv      (Code, FullName, Block)
#     * *_features_labeled_cut<NN>.csv (… plus x, y, radius), for each cutoff
# ======================================================================

suppressPackageStartupMessages(library(mixOmics))

## ----------------- Config -----------------
OUT_DIR   <- "cca_coords_exports"
dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(42)

# Canonical components and tuning grid
NCMP      <- 2
GRID      <- seq(0.1, 0.5, by = 0.1)

# Labeling cutoffs (radius on variable circle). One labeled CSV per cutoff.
CUTOFFS   <- c(0.70, 0.75, 0.80)

## ----------------- Helpers -----------------
sanitize_df <- function(df){
  cn <- colnames(df)
  cn <- iconv(cn, to = "ASCII//TRANSLIT"); cn[is.na(cn)] <- "unnamed"
  cn <- gsub("[^[:alnum:]_. -]", "_", cn, perl = TRUE)
  cn <- gsub("_+", "_", cn); cn <- trimws(cn)
  cn[cn == ""] <- "unnamed"
  colnames(df) <- cn
  df
}

align_blocks <- function(X, Y){
  common <- intersect(rownames(X), rownames(Y))
  list(
    X = X[common, , drop = FALSE],
    Y = Y[common, , drop = FALSE]
  )
}

# Short codes per block for compact plotting labels
make_codes <- function(names_vec, block_name){
  prefix <- switch(block_name,
                   "Bile"   = "B",
                   "16S"    = "T",
                   "MetaG"  = "M",
                   "MetaT"  = "X",
                   "RNAseq" = "R",
                   "")
  if (prefix == "") {
    data.frame(Code = names_vec, FullName = names_vec, stringsAsFactors = FALSE)
  } else {
    data.frame(Code = paste0(prefix, seq_along(names_vec)),
               FullName = names_vec, stringsAsFactors = FALSE)
  }
}

# Compute variable-circle coordinates (correlations to canonical variates)
# Returns a list with coords/radius for X and Y blocks.
compute_coords <- function(fit, X, Y){
  VX <- fit$variates$X[, 1:2, drop = FALSE]
  VY <- fit$variates$Y[, 1:2, drop = FALSE]

  cx1 <- as.numeric(cor(X, VX[,1])); names(cx1) <- colnames(X)
  cx2 <- as.numeric(cor(X, VX[,2])); names(cx2) <- colnames(X)
  cy1 <- as.numeric(cor(Y, VY[,1])); names(cy1) <- colnames(Y)
  cy2 <- as.numeric(cor(Y, VY[,2])); names(cy2) <- colnames(Y)

  list(
    X = data.frame(var = colnames(X), x = cx1, y = cx2,
                   radius = sqrt(cx1^2 + cx2^2), stringsAsFactors = FALSE),
    Y = data.frame(var = colnames(Y), x = cy1, y = cy2,
                   radius = sqrt(cy1^2 + cy2^2), stringsAsFactors = FALSE)
  )
}

# Label indices by radius cutoff; guarantee >=1 per side
label_indices <- function(radii, cutoff){
  idx <- which(radii >= cutoff)
  if (length(idx) == 0) idx <- which.max(radii)
  idx
}

## ---------------- Runner for one pair ----------------
run_pair_export <- function(x_mat, y_mat, nameX, nameY){
  # 0) clean + align
  X <- sanitize_df(x_mat); Y <- sanitize_df(y_mat)
  al <- align_blocks(X, Y); X <- al$X; Y <- al$Y
  if (nrow(X) < 3) {
    message("Skip ", nameX, " vs ", nameY, ": <3 common samples")
    return(invisible(NULL))
  }
  X[is.na(X)] <- 0; Y[is.na(Y)] <- 0

  # 1) tune + fit on FULL blocks (no pre-filter)
  tune <- tune.rcc(X, Y, grid1 = GRID, grid2 = GRID, validation = "loo")
  fit  <- rcc(X, Y, ncomp = min(NCMP, ncol(X), ncol(Y)),
              lambda1 = tune$opt.lambda1, lambda2 = tune$opt.lambda2)

  # 2) build code maps (short codes only for Bile/16S/MetaG/MetaT/RNAseq)
  mapX <- make_codes(colnames(X), nameX)
  mapY <- make_codes(colnames(Y), nameY)
  codesX <- setNames(mapX$Code, mapX$FullName)
  codesY <- setNames(mapY$Code, mapY$FullName)

  # 3) compute coordinates and radius
  coords <- compute_coords(fit, X, Y)
  cX <- coords$X; cY <- coords$Y

  # 4) write "all features" CSV (no cutoff)
  all_out <- rbind(
    data.frame(Block = nameX,
               Code = codesX[cX$var],
               FullName = cX$var,
               stringsAsFactors = FALSE),
    data.frame(Block = nameY,
               Code = codesY[cY$var],
               FullName = cY$var,
               stringsAsFactors = FALSE)
  )
  all_path <- file.path(OUT_DIR, sprintf("%s_VS_%s_features_all.csv", nameX, nameY))
  write.csv(all_out, all_path, row.names = FALSE)

  # 5) for each cutoff, write "labeled" CSV (subset by radius)
  for (cut in CUTOFFS){
    ix <- label_indices(cX$radius, cut)
    iy <- label_indices(cY$radius, cut)

    labX <- cX[ix, , drop = FALSE]
    labY <- cY[iy, , drop = FALSE]

    labeled_out <- rbind(
      data.frame(Block = nameX,
                 Code = codesX[labX$var],
                 FullName = labX$var,
                 x = labX$x, y = labX$y, radius = labX$radius,
                 stringsAsFactors = FALSE),
      data.frame(Block = nameY,
                 Code = codesY[labY$var],
                 FullName = labY$var,
                 x = labY$x, y = labY$y, radius = labY$radius,
                 stringsAsFactors = FALSE)
    )

    lab_path <- file.path(
      OUT_DIR, sprintf("%s_VS_%s_features_labeled_cut%02d.csv", nameX, nameY, round(cut*100))
    )
    write.csv(labeled_out, lab_path, row.names = FALSE)
  }

  message(sprintf("Exported: %s vs %s  (all + labeled for %s)",
                  nameX, nameY, paste0(sprintf("%.2f", CUTOFFS), collapse = ", ")))
}

## ----------------- Load matrices -----------------
MetaT  <- read.csv("MetaT_CCA.csv",  row.names = 1, check.names = FALSE)
Bile   <- read.csv("Bile_CCA.csv",   row.names = 1, check.names = FALSE)
MetaG  <- read.csv("MetaG_CCA.csv",  row.names = 1, check.names = FALSE)
RNAseq <- read.csv("RNAseq_CCA.csv", row.names = 1, check.names = FALSE)
S16    <- read.csv("16S_CCA.csv",    row.names = 1, check.names = FALSE)

blocks <- list(MetaT = MetaT, Bile = Bile, MetaG = MetaG, RNAseq = RNAseq, `16S` = S16)

## ----------------- Run all 10 pairs -----------------
pairs <- combn(names(blocks), 2, simplify = FALSE)
for (pr in pairs){
  xname <- pr[1]; yname <- pr[2]
  run_pair_export(blocks[[xname]], blocks[[yname]], xname, yname)
}

## ----------------- Save run info -----------------
sink(file.path(OUT_DIR, "run_info.txt"))
cat("rCCA coordinate export\n")
cat("Date: ", format(Sys.time()), "\n\n", sep = "")
cat("NCMP: ", NCMP, "\n", sep = "")
cat("GRID: ", paste(GRID, collapse = ", "), "\n", sep = "")
cat("CUTOFFS: ", paste(CUTOFFS, collapse = ", "), "\n\n", sep = "")
print(sessionInfo())
sink()