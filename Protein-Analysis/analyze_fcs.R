# FCS Analysis Script in R
# Automates analysis of TwistBio Expression Test data
# Generates RidgePlots, DotPlots, and Summary Statistics

# ---- 1. SETUP & DEPENDENCIES ----

# Function to automatically install missing packages
ensure_package <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

# Install critical dependencies
ensure_package("tidyverse")
ensure_package("ggridges")
ensure_package("scales")
ensure_package("hexbin")
ensure_package("flowCore", bioc = TRUE)
ensure_package("flowStats", bioc = TRUE)
ensure_package("ggcyto", bioc = TRUE)
ensure_package("openxlsx")

# Load Libraries
suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
  library(ggcyto)
  library(tidyverse)
  library(ggridges)
  library(scales)
  library(openxlsx)
})

# ---- 2. CONFIGURATION ----

# Paths
# Note: R usually runs from the script directory or project root. 
# Adjusting based on user context: "DNA Analysis" folder
base_dir <- getwd()
data_dir <- file.path(base_dir, "Local", "TwistBio Expression Test")
output_dir <- file.path(base_dir, "Local", "FCS_Analysis_Results_R")
iso_dir <- file.path(output_dir, "Individual_Plots")

if (!dir.exists(data_dir)) {
  stop(paste("Data directory not found:", data_dir))
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(iso_dir)) dir.create(iso_dir, recursive = TRUE)

# Channels will be resolved dynamically after loading data
CH_FSC_A <- "FSC-A"
CH_SSC_A <- "SSC-A"
CH_FSC_H <- "FSC-H"
# Placeholders
CH_FITC  <- NULL 
CH_APC   <- NULL

# ---- 3. DATA LOADING & TRANSFORMATION ----

message("Searching for FCS files...")
fcs_files <- list.files(data_dir, pattern = "\\.fcs$", full.names = TRUE)

if (length(fcs_files) == 0) {
  stop("No FCS files found in directory!")
}

message(paste("Found", length(fcs_files), "files. Loading..."))

# Load into a flowSet
fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)

# Clean up sample names (remove path and extension)
sampleNames(fs) <- gsub(".fcs", "", basename(sampleNames(fs)))

# ---- RESOLVE CHANNELS EXPERIMENTALLY ----
# Check what columns we actually have
cols <- colnames(fs[[1]])
params <- pData(parameters(fs[[1]]))

# Helper to find channel by string matching
find_ch <- function(targets, fallback) {
  # 1. Check direct match in colnames
  for(t in targets) {
    if(t %in% cols) return(t)
  }
  # 2. Check descriptions (formatted as "Desc" or "desc")
  # Some flowFrames use "desc", some "Desc"
  desc_col <- if("desc" %in% colnames(params)) "desc" else "Desc"
  if(!is.null(desc_col)) {
    match_idx <- grep(paste(targets, collapse="|"), params[[desc_col]], ignore.case=TRUE)
    if(length(match_idx) > 0) return(params$name[match_idx[1]])
  }
  # 3. Fallback
  if(fallback %in% cols) return(fallback)
  return(NULL)
}

CH_FITC <- find_ch(c("FITC-A", "FITC"), "FL1-A")
CH_APC  <- find_ch(c("APC-A", "APC", "R660-A", "RED-A"), "FL4-A") # FL4 is standard for APC/Red on Accuri/FacsCalibur

if(is.null(CH_FITC) || is.null(CH_APC)) {
  message("Available Columns: ", paste(cols, collapse=", "))
  stop("Could not automatically resolve FITC or APC channels! Please update script manually.")
}

message(paste("Resolved Channels -> FITC:", CH_FITC, "| APC:", CH_APC))


# Transform (Logicle is better than simple Log for flow data)
# estimateLogicle failed due to missing keywords in this dataset ($P3N)
# Use standard logicle transform instead
message("Applying Logicle Transformation...")
lgcl <- logicleTransform(w = 0.5, t = 262144, m = 4.5, a = 0)
trans <- transformList(c(CH_FITC, CH_APC), lgcl)
fs_trans <- transform(fs, trans)

# ---- 4. GATING STRATEGY ----

# A. Debris Gate (FSC-A vs SSC-A)
# We can use a static rectangle or a data-driven gate. 
# Replicating Python logic: FSC > 10000 && SSC > 2000 roughly
# But using a polygon gate is cleaner in flowCore
g_debris <- rectangleGate(filterId="DebrisFilter", 
                          "FSC-A" = c(10000, 262144), 
                          "SSC-A" = c(2000, 262144))

# B. Singlet Gate (FSC-A vs FSC-H)
# Linear relationship expected. 
# We'll create a polygon gate around the diagonal
# Or simpler: ratio filter
singlet_filter <- function(fr) {
  fsc_a <- exprs(fr)[, CH_FSC_A]
  fsc_h <- exprs(fr)[, CH_FSC_H]
  ratio <- fsc_h / fsc_a
  
  # Standard Singlet Logic: H ~= A, so ratio ~= 1 (or specific factor based on instrument scaling)
  # Python script used Median Ratio +/- 20%
  med_r <- median(ratio, na.rm = TRUE)
  return(ratio > (med_r * 0.8) & ratio < (med_r * 1.2))
}

# Apply Gates
message("Gating data (Debris -> Singlets)...")
gs <- GatingSet(fs_trans) # Move to GatingSet for easier handling

# 1. Add Debris Gate
gs_pop_add(gs, g_debris, parent = "root")

# 2. Add Singlet Gate (Custom Filter per sample for robustness)
# Simple rectangle on the Ratio is hard in GatingSet structure without derived params
# We'll just filter the flowSet manually for analysis simplicity or use a generic polygon
# Let's use flowStats logic for singlets if possible, or sticking to the manual manual logic for speed:

singlet_list <- fsApply(fs_trans, function(fr) {
  # 1. Debris
  fr_deb <- Subset(fr, g_debris)
  
  # 2. Singlets
  if (nrow(fr_deb) < 100) return(fr_deb) # Too few events
  mask <- singlet_filter(fr_deb)
  return(fr_deb[mask, ])
})
fs_singlets <- as(singlet_list, "flowSet")

# Recalculate GatingSet on clean data for plotting convenience
gs_clean <- GatingSet(fs_singlets)


# ---- 5. DETERMINE THRESHOLDS (CONTROLS) ----

# Identify Negative Controls
neg_pattern <- "Negative Control"
is_neg <- grepl(neg_pattern, sampleNames(fs_singlets))
fs_neg <- fs_singlets[is_neg]

thresh_fitc <- 100 # Default fallback
thresh_apc <- 100

if (length(fs_neg) > 0) {
  message("Calculating thresholds from Negative Controls...")
  # Concatenate all negative control events securely
  message(paste("Negative Control Frames:", length(fs_neg)))
  expr_neg_list <- fsApply(fs_neg, exprs, simplify = FALSE)
  expr_neg <- do.call(rbind, expr_neg_list)
  
  message(paste("Combined Neg Events:", nrow(expr_neg)))
  message(paste("Columns available:", paste(colnames(expr_neg), collapse=", ")))

  if (nrow(expr_neg) < 10) {
      warning("Not enough negative events for robust stats!")
  }
  
  # Calculate 99th percentile
  thresh_fitc <- quantile(expr_neg[, CH_FITC], 0.99)
  thresh_apc <- quantile(expr_neg[, CH_APC], 0.99)
  
  # Stats for Fold Change
  neg_mfi_fitc <- median(expr_neg[, CH_FITC])
  neg_mfi_apc  <- median(expr_neg[, CH_APC])
  neg_sd_fitc  <- mad(expr_neg[, CH_FITC]) # Robust SD
} else {
  warning("No Negative Control files found! Using arbitrary thresholds.")
  neg_mfi_fitc <- 1
  neg_mfi_apc <- 1
  neg_sd_fitc <- 1
}

# Add Quadrant Gates to GatingSet
# Use do.call to pass dynamic names
qg_list <- list(thresh_fitc, thresh_apc)
names(qg_list) <- c(CH_FITC, CH_APC)
qg_list$filterId <- "Quad"
g_quad <- do.call(quadGate, qg_list)
gs_pop_add(gs_clean, g_quad)
recompute(gs_clean)

# ---- 6. STATISTICAL ANALYSIS ----

message("Calculating Summary Statistics...")

# Function to get stats for a single frame
calc_stats <- function(fr) {
  d <- exprs(fr)
  total <- nrow(d)
  
  # Percentages
  p_fitc <- sum(d[, CH_FITC] > thresh_fitc) / total * 100
  p_apc <- sum(d[, CH_APC] > thresh_apc) / total * 100
  p_double <- sum(d[, CH_FITC] > thresh_fitc & d[, CH_APC] > thresh_apc) / total * 100
  
  # MFIs
  mfi_fitc <- median(d[, CH_FITC])
  mfi_apc <- median(d[, CH_APC])
  
  # Metrics
  fc_fitc <- mfi_fitc / neg_mfi_fitc
  si_fitc <- (mfi_fitc - neg_mfi_fitc) / (2 * neg_sd_fitc)
  
  return(c(
    TotalEvents = total,
    Pct_FITC = p_fitc,
    Pct_APC = p_apc,
    Pct_Double = p_double,
    MFI_FITC = mfi_fitc,
    MFI_APC = mfi_apc,
    FoldChange_FITC = fc_fitc,
    StainIndex_FITC = si_fitc
  ))
}

stats_list <- fsApply(fs_singlets, calc_stats)
df_stats <- as.data.frame(stats_list)
df_stats$Sample <- rownames(df_stats)
df_stats <- df_stats %>% select(Sample, everything())

# Save Stats
write.xlsx(df_stats, file.path(output_dir, "Summary_Statistics.xlsx"))
write.csv(df_stats, file.path(output_dir, "Summary_Statistics.csv"))


# ---- 7. VISUALIZATION ----

# Convert flowSet to tidy dataframe for ggplot (downsample if huge)
message("Preparing data for plots...")
df_all <- fsApply(fs_singlets, function(fr) {
  d <- exprs(fr)
  # Downsample for speed in heavy plots if needed (e.g. max 5k events per file)
  if(nrow(d) > 5000) d <- d[sample(nrow(d), 5000), ]
  df <- as.data.frame(d)
  df$Sample <- identifier(fr)
  return(df)
})
df_all <- do.call(rbind, df_all)

# Extract Group Name (First word)
df_all$Group <- sapply(strsplit(df_all$Sample, " "), `[`, 1)

# A. Aggregate Ridge Plot (Binding)
message("Generating Ridgeline Plots...")

# Filter outliers for cleaner plots (0.5% - 99.5%)
# Use tidy eval for dynamic column names
fitc_sym <- sym(CH_FITC)

df_ridge <- df_all %>%
  group_by(Sample) %>%
  filter(!!fitc_sym > quantile(!!fitc_sym, 0.01) & !!fitc_sym < quantile(!!fitc_sym, 0.99)) %>%
  ungroup()

# Sort samples by MFI for the plot
order_mfi <- df_stats %>% arrange(MFI_FITC) %>% pull(Sample)
df_ridge$Sample <- factor(df_ridge$Sample, levels = order_mfi)

p_ridge <- ggplot(df_ridge, aes(x = !!fitc_sym, y = Sample, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_x_continuous(name = paste("Binding Logicle (", CH_FITC, ")", sep="")) +
  scale_fill_viridis_c(name = "Intensity", option = "C") +
  geom_vline(xintercept = thresh_fitc, linetype = "dashed", color = "red") +
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(title = "Aggregate Binding Distribution", subtitle = "Sorted by MFI intensity")

ggsave(file.path(output_dir, "Aggregate_Ridgeline_Binding.png"), p_ridge, width = 10, height = max(6, length(unique(df_ridge$Sample)) * 0.3))


# B. Individual Dot Plots (Expression vs Binding)
message("Generating Individual Plots...")

# Using ggcyto for efficiency
# We iterate to save individual files
for (smp in sampleNames(fs_singlets)) {
  
  # Grab stats
  st <- df_stats[df_stats$Sample == smp, ]
  fr <- fs_singlets[[smp]]
  events <- nrow(fr)
  
  message(paste("Plotting", smp, "| Events:", events))
  
  p <- ggcyto(fr, aes(x = !!sym(CH_APC), y = !!sym(CH_FITC)))
  
  # Adaptive plotting based on event count
  if (events < 2000) {
    p <- p + geom_point(alpha = 0.6, size = 0.8)
  } else {
    p <- p + geom_hex(bins = 100)
  }
  
  p <- p + 
    geom_gate(g_quad) + 
    geom_stats(type = "percent", adjust = 0.7) +
    # Data is already transformed, so we use linear scale on the logicle values
    scale_x_continuous() +
    scale_y_continuous() +
    labs(title = smp,
         subtitle = paste0("Double+: ", round(st$Pct_Double, 1), "% | FC: ", round(st$FoldChange_FITC, 1)),
         x = "Expression (APC)", y = "Binding (FITC)") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(iso_dir, paste0(smp, "_dotplot.png")), p, width = 6, height = 6)
}

# C. Summary Bar Plots
message("Generating Summary Charts...")

# Fold Change
p_fc <- ggplot(df_stats, aes(x = reorder(Sample, FoldChange_FITC), y = FoldChange_FITC, fill = FoldChange_FITC)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c() +
  geom_hline(yintercept = 1, linetype="dashed") +
  labs(title = "Binding Fold Change (vs Neg)", y = "Fold Change", x = "Sample") +
  theme_minimal()

ggsave(file.path(output_dir, "Summary_FoldChange.png"), p_fc, width = 8, height = max(6, nrow(df_stats)*0.2))

# Double Positive %
p_dp <- ggplot(df_stats, aes(x = reorder(Sample, Pct_Double), y = Pct_Double, fill = Pct_Double)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Double Positive Population", y = "% Double Positive", x = "Sample") +
  theme_minimal()

ggsave(file.path(output_dir, "Summary_DoublePos.png"), p_dp, width = 8, height = max(6, nrow(df_stats)*0.2))

message("Analysis Complete! Check results in '../Local/FCS_Analysis_Results_R'")
