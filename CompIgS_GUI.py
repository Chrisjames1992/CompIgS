#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
================================================================================
BCR ANALYSIS SUITE - Integrated CompIgS + FastBCR Analysis Platform
================================================================================

This application combines two independent BCR analysis methods:

1. CompIgS Method (Python/VDJ Matching):
   - Uses V+D+J gene matching to identify clonotypes
   - Requires IMGT output files (StatClonotype, mutation stats, AA sequences)
   - Identifies shared/unique clones based on exact VDJ matching

2. FastBCR Method (R/VJ + Junction Similarity):
   - Clusters sequences within each isotype using V+J genes and k-mer seeding
   - Identifies shared clusters between isotypes using V+J genes and junction AA similarity
   - Requires AIRR-formatted TSV files
   - Identifies shared clones based on configurable junction similarity threshold

Both methods can identify divergent clones with configurable thresholds.
================================================================================
"""

import sys
import os
import subprocess
import threading
import json
import logging
import time
import traceback
from pathlib import Path
from datetime import datetime
from io import StringIO

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QLineEdit, QFileDialog, QProgressBar, QWidget, QFrame, QMessageBox, 
    QGroupBox, QComboBox, QCheckBox, QTextEdit, QTabWidget, QSplitter,
    QListWidget, QListWidgetItem, QDialog, QPlainTextEdit, QSpinBox,
    QDoubleSpinBox, QRadioButton, QButtonGroup, QScrollArea, QGridLayout
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer, QProcess
from PyQt5.QtGui import QFont, QPalette, QColor, QLinearGradient, QBrush, QPixmap

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False


# Scientific computing and data analysis
import pandas as pd
import numpy as np
import glob

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Bioinformatics
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# Peptide analysis
try:
    from peptides import Peptide
except ImportError:
    # Fallback if peptides package not installed
    class Peptide:
        def __init__(self, sequence):
            self.sequence = sequence
        def charge(self, pH=7.0, pKscale="Lehninger"):
            return 0  # Returns 0 if package unavailable

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Configuration settings for the application"""
    APP_NAME = "BCR Analysis Suite"
    VERSION = "2.0.0"
    
    # Default parameters
    DEFAULT_JUNCTION_SIMILARITY = 80
    DEFAULT_MIN_AA_DIFF = 2
    DEFAULT_MIN_CLUSTER_DEPTH = 3
    DEFAULT_CDR3_THRESHOLD = 2
    DEFAULT_CDR3CDR2_THRESHOLD = 2
    
    @staticmethod
    def check_r_installation():
        """Check if R is installed and available"""
        try:
            # Hide console window on Windows
            if sys.platform == "win32":
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                startupinfo.wShowWindow = subprocess.SW_HIDE
                creationflags = subprocess.CREATE_NO_WINDOW
            else:
                startupinfo = None
                creationflags = 0
            
            result = subprocess.run(
                ["Rscript", "--version"], 
                capture_output=True, 
                text=True,
                timeout=10,
                startupinfo=startupinfo,
                creationflags=creationflags
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    @staticmethod
    def check_fastbcr_installed():
        """Check if fastBCR package is installed in R"""
        try:
            # Hide console window on Windows
            if sys.platform == "win32":
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                startupinfo.wShowWindow = subprocess.SW_HIDE
                creationflags = subprocess.CREATE_NO_WINDOW
            else:
                startupinfo = None
                creationflags = 0
            
            result = subprocess.run(
                ["Rscript", "-e", "library(fastBCR); cat('OK')"],
                capture_output=True,
                text=True,
                timeout=30,
                startupinfo=startupinfo,
                creationflags=creationflags
            )
            return "OK" in result.stdout
        except:
            return False


# =============================================================================
# R ANALYSIS THREAD (FastBCR Method)
# =============================================================================

class FastBCRAnalysisThread(QThread):
    """Thread to run FastBCR R analysis as subprocess"""
    progress_updated = pyqtSignal(str)
    analysis_completed = pyqtSignal(str)
    error_occurred = pyqtSignal(str)
    
    def __init__(self, iso1_file, iso2_file, output_dir, parameters):
        super().__init__()
        self.iso1_file = iso1_file
        self.iso2_file = iso2_file
        self.output_dir = output_dir
        self.parameters = parameters
        self.process = None
        self._is_running = True
    
    def run(self):
        try:
            # Create the R script with user parameters
            temp_script = self.create_r_script()
            
            self.progress_updated.emit("=" * 60)
            self.progress_updated.emit("Starting FastBCR Analysis...")
            self.progress_updated.emit("=" * 60)
            self.progress_updated.emit(f"Isotype 1 file: {self.iso1_file}")
            self.progress_updated.emit(f"Isotype 2 file: {self.iso2_file}")
            self.progress_updated.emit(f"Output directory: {self.output_dir}")
            self.progress_updated.emit(f"Parameters: {json.dumps(self.parameters, indent=2)}")
            self.progress_updated.emit("-" * 60)
            
            # Configure subprocess to hide R console window on Windows
            if sys.platform == "win32":
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                startupinfo.wShowWindow = subprocess.SW_HIDE
                creationflags = subprocess.CREATE_NO_WINDOW
            else:
                startupinfo = None
                creationflags = 0
            
            # Run R script
            self.process = subprocess.Popen(
                ["Rscript", "--vanilla", temp_script],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True,
                startupinfo=startupinfo,
                creationflags=creationflags
            )
                
            # Stream output
            for line in iter(self.process.stdout.readline, ''):
                if not self._is_running:
                    break
                if line:
                    self.progress_updated.emit(line.strip())
            
            self.process.wait()
            
            # Clean up temp script
            if os.path.exists(temp_script):
                os.remove(temp_script)
            
            if self.process.returncode == 0:
                self.analysis_completed.emit("FastBCR analysis completed successfully!")
            else:
                self.error_occurred.emit(f"R analysis failed with code {self.process.returncode}")
                
        except Exception as e:
            self.error_occurred.emit(f"Error running R analysis: {str(e)}\n{traceback.format_exc()}")
    
    def create_r_script(self):
        """Create a temporary R script with user parameters"""
        
        iso1_name = self.parameters.get('iso1_name', 'IgE')
        iso2_name = self.parameters.get('iso2_name', 'IgG1')
        junction_threshold = self.parameters.get('junction_threshold', 80)
        min_aa_diff = self.parameters.get('min_aa_diff', 2)
        min_depth = self.parameters.get('min_depth', 3)
        productive_only = self.parameters.get('productive_only', False)
        n_top_shared = self.parameters.get('n_top_shared', 10)
        n_top_divergent = self.parameters.get('n_top_divergent', 10)
        
        script_content = f'''
################################################################################
#                BCR REPERTOIRE ANALYSIS: {iso1_name} vs {iso2_name} COMPLETE PIPELINE
#                        ENHANCED VERSION WITH ADDITIONAL VISUALIZATIONS
################################################################################
#
# This script performs:
# 1. Data loading and preprocessing for {iso1_name} and {iso2_name}
# 2. BCR clustering analysis
# 3. SHM (Somatic Hypermutation) calculation
# 4. Shared cluster identification between {iso1_name} and {iso2_name}
# 5. Divergent {iso1_name} clone identification (ALL mutation levels)
# 6. Lineage tree visualization for shared and divergent clusters
# 7. MSA alignments showing divergent {iso1_name} sequences
# 8. Additional visualizations (V/J usage, junction similarity, etc.)
#
# Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
#######################################################################################

# =============================================================================------
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================-----

rm(list = ls())

# Enable warning capture for detailed reporting
options(warn = 1)  # Print warnings as they occur
warning_log <- c()  # Store warnings for summary

# Custom warning handler to capture details
withCallingHandlers_warnings <- function() {{
  last_warning <<- NULL
}}

# Track start time and initial memory for performance reporting
analysis_start_time <- Sys.time()
initial_memory <- gc(reset = TRUE)
cat("[MEMORY] Analysis started at:", format(analysis_start_time, "%Y-%m-%d %H:%M:%S"), "\\n")
cat("[MEMORY] Initial memory cleared\\n\\n")

# Track start time and initial memory for performance reporting
analysis_start_time <- Sys.time()
initial_memory <- gc(reset = TRUE)
cat("[MEMORY] Analysis started at:", format(analysis_start_time, "%Y-%m-%d %H:%M:%S"), "\\n")
cat("[MEMORY] Initial memory cleared\\n\\n")

# Load required libraries
required_packages <- c("fastBCR", "ggplot2", "dplyr", "stringdist", "gridExtra", 
                       "RColorBrewer", "VennDiagram", "ape", "ggtree", "Biostrings",
                       "msa", "grid", "ggpubr", "pheatmap", "seqinr", "tidyr",
                       "data.table", "parallel")

for(pkg in required_packages) {{
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {{
    message(paste("Installing package:", pkg))
    if(pkg %in% c("Biostrings", "msa", "ggtree")) {{
      if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    }} else {{
      install.packages(pkg)
    }}
    library(pkg, character.only = TRUE)
  }}
}}

# -----------------------------------------------------------------------------
# PERFORMANCE OPTIMIZATION SETUP (WINDOWS COMPATIBLE)
# -----------------------------------------------------------------------------
library(data.table)
library(parallel)

N_CORES <- max(1, detectCores() - 1)
setDTthreads(N_CORES)

# Detect OS and setup cluster for Windows
IS_WINDOWS <- .Platform$OS.type == "windows"

if(IS_WINDOWS) {{
  cat("[PERFORMANCE] Windows detected - using parLapply with", N_CORES, "cores\n")
  CL <- makeCluster(N_CORES)
}} else {{
  cat("[PERFORMANCE] Unix/Mac detected - using mclapply with", N_CORES, "cores\n")
  CL <- NULL
}}

# Universal parallel lapply function
par_lapply <- function(X, FUN, ...) {{
  if(IS_WINDOWS) {{
    parLapply(CL, X, FUN, ...)
  }} else {{
    mclapply(X, FUN, ..., mc.cores = N_CORES)
  }}
}}

cat("[PERFORMANCE] Parallel processing initialized\n\n")

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------

ISO1_FILE <- "{self.iso1_file.replace(os.sep, '/')}"
ISO2_FILE <- "{self.iso2_file.replace(os.sep, '/')}"

ISO1_SAMPLE_NAME <- "Sample_{iso1_name}"
ISO2_SAMPLE_NAME <- "Sample_{iso2_name}"

ISO1_NAME <- "{iso1_name}"
ISO2_NAME <- "{iso2_name}"

# Analysis parameters
JUNCTION_SIMILARITY_THRESHOLD <- {junction_threshold}
DIVERGENT_MIN_DIFFERENCES <- {min_aa_diff}
PRODUCTIVE_ONLY <- {str(productive_only).upper()}

# Clustering parameters
MIN_DEPTH_THRESHOLD <- {min_depth}
MAX_DEPTH_THRESHOLD <- 1000
OVERLAP_THRESHOLD <- 0.1
CONSENSUS_THRESHOLD <- 0.80

# Output directory
OUTPUT_DIR <- "{self.output_dir.replace(os.sep, '/')}"
setwd(OUTPUT_DIR)

# Number of top clusters to analyze in detail
N_TOP_SHARED_CLUSTERS <- {n_top_shared}
N_TOP_DIVERGENT_CLUSTERS <- {n_top_divergent}

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

calculate_shm <- function(seq, germline) {{
  if(is.na(seq) || is.na(germline)) return(NA)
  seq_chars <- strsplit(as.character(seq), "")[[1]]
  germ_chars <- strsplit(as.character(germline), "")[[1]]
  min_len <- min(length(seq_chars), length(germ_chars))
  if(min_len == 0) return(NA)
  seq_chars <- seq_chars[1:min_len]
  germ_chars <- germ_chars[1:min_len]
  valid_positions <- which(seq_chars != "." & germ_chars != "." & 
                             seq_chars != "N" & germ_chars != "N" &
                             seq_chars != "-" & germ_chars != "-")
  if(length(valid_positions) == 0) return(NA)
  mismatches <- sum(seq_chars[valid_positions] != germ_chars[valid_positions])
  return(mismatches / length(valid_positions) * 100)
}}

extract_gene <- function(gene_call) {{
  gene_call <- as.character(gene_call)
  return(gsub("\\\\*.*", "", gene_call))
}}

calculate_junction_similarity <- function(junction1, junction2) {{
  if(is.na(junction1) || is.na(junction2) || junction1 == "" || junction2 == "") return(0)
  dist <- stringdist(junction1, junction2, method = "lv")
  max_len <- max(nchar(junction1), nchar(junction2))
  return((1 - (dist / max_len)) * 100)
}}

calculate_aa_differences <- function(seq1, seq2) {{
  if(is.na(seq1) || is.na(seq2) || seq1 == "" || seq2 == "") return(NA)
  len1 <- nchar(seq1)
  len2 <- nchar(seq2)
  if(len1 == 0 || len2 == 0) return(NA)
  if(len1 == len2) {{
    seq1_chars <- strsplit(seq1, "")[[1]]
    seq2_chars <- strsplit(seq2, "")[[1]]
    return(sum(seq1_chars != seq2_chars))
  }} else {{
    return(stringdist(seq1, seq2, method = "lv"))
  }}
}}

get_cluster_consensus_genes <- function(cluster_data) {{
  v_genes <- sapply(cluster_data$v_call, extract_gene)
  j_genes <- sapply(cluster_data$j_call, extract_gene)
  return(list(
    v_gene = names(sort(table(v_genes), decreasing = TRUE))[1],
    j_gene = names(sort(table(j_genes), decreasing = TRUE))[1]
  ))
}}

is_divergent_clone <- function(junction, reference_junctions, min_differences) {{
  differences <- sapply(reference_junctions, function(ref_seq) {{
    calculate_aa_differences(junction, ref_seq)
  }})
  differences <- differences[!is.na(differences)]
  if(length(differences) == 0) return(FALSE)
  return(min(differences) >= min_differences)
}}


# Helper function to save plots as both PDF and PNG
save_plot <- function(plot_obj, filename_base, width = 10, height = 8, dpi = 150) {{
  pdf_file <- file.path(OUTPUT_DIR, paste0(filename_base, ".pdf"))
  png_file <- file.path(OUTPUT_DIR, paste0(filename_base, ".png"))
  
  ggsave(pdf_file, plot = plot_obj, width = width, height = height)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  
  cat("  SAVED:", filename_base, ".pdf and .png\\n")
}}

# Helper function to save plots as both PDF and PNG
save_plot_dual <- function(plot_obj, filename_base, width = 10, height = 8, dpi = 150) {{
  pdf_file <- file.path(OUTPUT_DIR, paste0(filename_base, ".pdf"))
  png_file <- file.path(OUTPUT_DIR, paste0(filename_base, ".png"))
  
  ggsave(pdf_file, plot = plot_obj, width = width, height = height)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  
  cat("  SAVED:", filename_base, ".pdf and .png\\n")
}}

get_clean_id <- function(seq_id) {{
  seq_id <- as.character(seq_id)
  nums <- regmatches(seq_id, gregexpr("[0-9]+", seq_id))[[1]]
  if(length(nums) > 0) {{
    return(nums[length(nums)])
  }}
  return(seq_id)
}}


# Helper function to run code and capture warnings with context
run_with_warning_capture <- function(expr, section_name) {{
  tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {{
        cat("[WARNING] In", section_name, ":", conditionMessage(w), "\\n")
        invokeRestart("muffleWarning")
      }}
    ),
    error = function(e) {{
      cat("[ERROR] In", section_name, ":", conditionMessage(e), "\\n")
    }}
  )
}}



cat("\\n")
cat("================================================================================\\n")
cat("         BCR REPERTOIRE ANALYSIS:", ISO1_NAME, "vs", ISO2_NAME, "- ENHANCED VERSION\\n")
cat("================================================================================\\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

# =============================================================================
# SECTION 2: DATA LOADING AND PREPROCESSING
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 2: DATA LOADING AND PREPROCESSING\\n")
cat("================================================================================\\n\\n")

raw_data_iso1 <- read.csv(ISO1_FILE, sep = "\\t", stringsAsFactors = FALSE)
cat("  -", ISO1_NAME, "sequences loaded:", nrow(raw_data_iso1), "\\n")

raw_data_iso2 <- read.csv(ISO2_FILE, sep = "\\t", stringsAsFactors = FALSE)
cat("  -", ISO2_NAME, "sequences loaded:", nrow(raw_data_iso2), "\\n\\n")

raw_sample_list_iso1 <- list(raw_data_iso1)
names(raw_sample_list_iso1) <- ISO1_SAMPLE_NAME

raw_sample_list_iso2 <- list(raw_data_iso2)
names(raw_sample_list_iso2) <- ISO2_SAMPLE_NAME


# Count non-productive sequences for exclusion reporting
cat("\\n================================================================================\\n")
cat("[EXCLUSION] NON-PRODUCTIVE SEQUENCES:\\n")
if("productive" %in% colnames(raw_data_iso1)) {{
  nonprod_iso1 <- sum(raw_data_iso1$productive == "F" | raw_data_iso1$productive == FALSE, na.rm = TRUE)
  prod_iso1 <- sum(raw_data_iso1$productive == "T" | raw_data_iso1$productive == TRUE, na.rm = TRUE)
  cat("[EXCLUSION]", ISO1_NAME, "non-productive sequences:", nonprod_iso1, "of", nrow(raw_data_iso1), "\\n")
  cat("[EXCLUSION]", ISO1_NAME, "productive sequences:", prod_iso1, "\\n")
}} else {{
  cat("[EXCLUSION]", ISO1_NAME, "productive column not found - assuming all productive\\n")
}}

if("productive" %in% colnames(raw_data_iso2)) {{
  nonprod_iso2 <- sum(raw_data_iso2$productive == "F" | raw_data_iso2$productive == FALSE, na.rm = TRUE)
  prod_iso2 <- sum(raw_data_iso2$productive == "T" | raw_data_iso2$productive == TRUE, na.rm = TRUE)
  cat("[EXCLUSION]", ISO2_NAME, "non-productive sequences:", nonprod_iso2, "of", nrow(raw_data_iso2), "\\n")
  cat("[EXCLUSION]", ISO2_NAME, "productive sequences:", prod_iso2, "\\n")
}} else {{
  cat("[EXCLUSION]", ISO2_NAME, "productive column not found - assuming all productive\\n")
}}
cat("================================================================================\\n\\n")




# =============================================================================
# SECTION 3: BCR CLUSTERING
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 3: BCR CLUSTERING\\n")
cat("================================================================================\\n\\n")

cat("Preprocessing and clustering", ISO1_NAME, "data...\\n")
pro_sample_list_iso1 <- data.preprocess(raw_data_list = raw_sample_list_iso1,
                                        productive_only = PRODUCTIVE_ONLY)

cluster_list_iso1 <- data.BCR.clusters(pro_sample_list_iso1,
                                       min_depth_thre = MIN_DEPTH_THRESHOLD,
                                       min_depth_thre_adjustment = TRUE,
                                       max_depth_thre = MAX_DEPTH_THRESHOLD,
                                       overlap_thre = OVERLAP_THRESHOLD,
                                       consensus_thre = CONSENSUS_THRESHOLD, 
                                       paired = FALSE,
                                       singletons_backtrack = TRUE)

clusters_iso1 <- cluster_list_iso1[[ISO1_SAMPLE_NAME]]
cat("  -", ISO1_NAME, "clusters found:", length(clusters_iso1), "\\n")

cat("Preprocessing and clustering", ISO2_NAME, "data...\\n")
pro_sample_list_iso2 <- data.preprocess(raw_data_list = raw_sample_list_iso2,
                                        productive_only = PRODUCTIVE_ONLY)

cluster_list_iso2 <- data.BCR.clusters(pro_sample_list_iso2,
                                       min_depth_thre = MIN_DEPTH_THRESHOLD,
                                       min_depth_thre_adjustment = TRUE,
                                       max_depth_thre = MAX_DEPTH_THRESHOLD,
                                       overlap_thre = OVERLAP_THRESHOLD,
                                       consensus_thre = CONSENSUS_THRESHOLD, 
                                       paired = FALSE,
                                       singletons_backtrack = TRUE)

clusters_iso2 <- cluster_list_iso2[[ISO2_SAMPLE_NAME]]
cat("  -", ISO2_NAME, "clusters found:", length(clusters_iso2), "\\n\\n")

# Get clustered/unclustered sequences
seqs_list_iso1 <- Clustered.seqs(pro_data_list = pro_sample_list_iso1, clusters_list = cluster_list_iso1)
clustered_seqs_iso1 <- seqs_list_iso1$clustered_seqs[[ISO1_SAMPLE_NAME]]
unclustered_seqs_iso1 <- seqs_list_iso1$unclustered_seqs[[ISO1_SAMPLE_NAME]]

seqs_list_iso2 <- Clustered.seqs(pro_data_list = pro_sample_list_iso2, clusters_list = cluster_list_iso2)
clustered_seqs_iso2 <- seqs_list_iso2$clustered_seqs[[ISO2_SAMPLE_NAME]]
unclustered_seqs_iso2 <- seqs_list_iso2$unclustered_seqs[[ISO2_SAMPLE_NAME]]

# Output exclusion statistics for log viewer
cat("\\n================================================================================\\n")
cat("[EXCLUSION] UNCLUSTERED/SINGLETON SEQUENCES:\\n")
cat("[EXCLUSION]", ISO1_NAME, "unclustered (singleton) sequences:", nrow(unclustered_seqs_iso1), "\\n")
cat("[EXCLUSION]", ISO2_NAME, "unclustered (singleton) sequences:", nrow(unclustered_seqs_iso2), "\\n")
cat("[EXCLUSION] Total singletons excluded from clustering:", nrow(unclustered_seqs_iso1) + nrow(unclustered_seqs_iso2), "\\n")
cat("================================================================================\\n\\n")

# =============================================================================
# SECTION 4: FASTBCR BUILT-IN VISUALIZATIONS
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 4: FASTBCR BUILT-IN VISUALIZATIONS\\n")
cat("================================================================================\\n\\n")

# Cluster visualization (bubble plot)
cat("Creating cluster bubble visualizations...\\n")
pdf(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_clusters_bubble.pdf")), width = 8, height = 8)
Clusters.visualization(pro_data_list = pro_sample_list_iso1, clusters_list = cluster_list_iso1, index = 1)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_clusters_bubble.png")), width = 800, height = 800, res = 100)
Clusters.visualization(pro_data_list = pro_sample_list_iso1, clusters_list = cluster_list_iso1, index = 1)
dev.off()

pdf(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_clusters_bubble.pdf")), width = 8, height = 8)
Clusters.visualization(pro_data_list = pro_sample_list_iso2, clusters_list = cluster_list_iso2, index = 1)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_clusters_bubble.png")), width = 800, height = 800, res = 100)
Clusters.visualization(pro_data_list = pro_sample_list_iso2, clusters_list = cluster_list_iso2, index = 1)
dev.off()

# V gene usage pie charts
cat("Creating V/J gene usage plots...\\n")
pdf(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_V_gene_usage.pdf")), width = 10, height = 8)
pie.freq.plot(clustered_seqs = clustered_seqs_iso1, colname = "v_call")
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_V_gene_usage.png")), width = 1000, height = 800, res = 100)
pie.freq.plot(clustered_seqs = clustered_seqs_iso1, colname = "v_call")
dev.off()

pdf(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_J_gene_usage.pdf")), width = 10, height = 8)
pie.freq.plot(clustered_seqs = clustered_seqs_iso1, colname = "j_call")
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_J_gene_usage.png")), width = 1000, height = 800, res = 100)
pie.freq.plot(clustered_seqs = clustered_seqs_iso1, colname = "j_call")
dev.off()

pdf(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_V_gene_usage.pdf")), width = 10, height = 8)
pie.freq.plot(clustered_seqs = clustered_seqs_iso2, colname = "v_call")
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_V_gene_usage.png")), width = 1000, height = 800, res = 100)
pie.freq.plot(clustered_seqs = clustered_seqs_iso2, colname = "v_call")
dev.off()

pdf(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_J_gene_usage.pdf")), width = 10, height = 8)
pie.freq.plot(clustered_seqs = clustered_seqs_iso2, colname = "j_call")
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_J_gene_usage.png")), width = 1000, height = 800, res = 100)
pie.freq.plot(clustered_seqs = clustered_seqs_iso2, colname = "j_call")
dev.off()

# V-J pairing heatmap
cat("Creating V-J pairing heatmaps...\\n")
pdf(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_VJ_pairing.pdf")), width = 12, height = 10)
vjpair.sample.plot(clustered_seqs = clustered_seqs_iso1)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_VJ_pairing.png")), width = 1200, height = 1000, res = 100)
vjpair.sample.plot(clustered_seqs = clustered_seqs_iso1)
dev.off()

pdf(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_VJ_pairing.pdf")), width = 12, height = 10)
vjpair.sample.plot(clustered_seqs = clustered_seqs_iso2)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_VJ_pairing.png")), width = 1200, height = 1000, res = 100)
vjpair.sample.plot(clustered_seqs = clustered_seqs_iso2)
dev.off()

# Junction length distribution
cat("Creating junction length plots...\\n")
pdf(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_junction_length.pdf")), width = 10, height = 6)
len.sample.plot(clustered_seqs = clustered_seqs_iso1)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_junction_length.png")), width = 1000, height = 600, res = 100)
len.sample.plot(clustered_seqs = clustered_seqs_iso1)
dev.off()

pdf(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_junction_length.pdf")), width = 10, height = 6)
len.sample.plot(clustered_seqs = clustered_seqs_iso2)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_junction_length.png")), width = 1000, height = 600, res = 100)
len.sample.plot(clustered_seqs = clustered_seqs_iso2)
dev.off()

# Clustered vs unclustered length comparison
pdf(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_clustered_vs_unclustered_length.pdf")), width = 10, height = 6)
len.clustered.plot(clustered_seqs = clustered_seqs_iso1, unclustered_seqs = unclustered_seqs_iso1)
dev.off()

png(file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_clustered_vs_unclustered_length.png")), width = 1000, height = 600, res = 100)
len.clustered.plot(clustered_seqs = clustered_seqs_iso1, unclustered_seqs = unclustered_seqs_iso1)
dev.off()

# =============================================================================
# SECTION 5: SHM ANALYSIS FOR ALL CLUSTERS (OPTIMIZED)
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 5: SOMATIC HYPERMUTATION (SHM) ANALYSIS\\n")
cat("================================================================================\\n\\n")

# Detect which SHM method to use based on available columns
detect_shm_method <- function(data) {{
  # Method 1: Calculate from alignments (IMGT style)
  has_alignments <- "germline_alignment" %in% colnames(data) && 
                    "sequence_alignment" %in% colnames(data) &&
                    sum(!is.na(data$germline_alignment) & data$germline_alignment != "") > 0
  
  # Method 2: Pre-calculated mutation rate (MiXCR style)
  has_mutation_rate <- "Nt mutations rate in V gene" %in% colnames(data) &&
                       sum(!is.na(data$`Nt mutations rate in V gene`)) > 0
  
  # Method 3: v_identity
  has_v_identity <- "v_identity" %in% colnames(data) &&
                    sum(!is.na(data$v_identity)) > 0
  
  if(has_alignments) {{
    return("alignment")
  }} else if(has_mutation_rate) {{
    return("mutation_rate")
  }} else if(has_v_identity) {{
    return("v_identity")
  }} else {{
    return("none")
  }}
}}

# Function to get SHM values based on method
get_shm_values <- function(data, idx, method) {{
  if(method == "alignment") {{
    # Calculate from sequence vs germline alignment
    shm_values <- mapply(calculate_shm, 
                         data$sequence_alignment[idx], 
                         data$germline_alignment[idx])
  }} else if(method == "mutation_rate") {{
    # Use pre-calculated mutation rate (convert to percentage)
    shm_values <- data$`Nt mutations rate in V gene`[idx] * 100
  }} else if(method == "v_identity") {{
    # Calculate from v_identity (1 - identity = mutation rate)
    shm_values <- (1 - data$v_identity[idx]) * 100
  }} else {{
    return(NULL)
  }}
  
  shm_values <- as.numeric(shm_values)
  shm_values <- shm_values[!is.na(shm_values)]
  return(shm_values)
}}

# Detect method for each isotype
shm_method_iso1 <- detect_shm_method(raw_data_iso1)
shm_method_iso2 <- detect_shm_method(raw_data_iso2)

cat("SHM calculation method:\\n")
cat("  ", ISO1_NAME, ":", shm_method_iso1, "\\n")
cat("  ", ISO2_NAME, ":", shm_method_iso2, "\\n\\n")

# Calculate SHM for ISO1 clusters
cat("Calculating SHM for", ISO1_NAME, "clusters (parallel processing)...\\n")

if(IS_WINDOWS) {{
  clusterExport(CL, c("clusters_iso1", "raw_data_iso1", "calculate_shm", 
                      "get_shm_values", "shm_method_iso1"), envir = environment())
}}

cluster_shm_iso1_list <- par_lapply(1:length(clusters_iso1), function(i) {{
  cluster_seqs <- clusters_iso1[[i]]
  idx <- which(raw_data_iso1$junction_aa %in% cluster_seqs$junction_aa)
  
  if(length(idx) > 0) {{
    shm_values <- get_shm_values(raw_data_iso1, idx, shm_method_iso1)
    
    if(!is.null(shm_values) && length(shm_values) > 0) {{
      return(data.frame(
        cluster_id = i, cluster_size = nrow(cluster_seqs),
        mean_shm = mean(shm_values),
        sd_shm = ifelse(length(shm_values) > 1, sd(shm_values), NA),
        min_shm = min(shm_values), max_shm = max(shm_values)
      ))
    }}
  }}
  return(NULL)
}})

cluster_shm_iso1 <- do.call(rbind, Filter(Negate(is.null), cluster_shm_iso1_list))
if(is.null(cluster_shm_iso1)) {{
  cluster_shm_iso1 <- data.frame(cluster_id=integer(), cluster_size=integer(), 
                                  mean_shm=numeric(), sd_shm=numeric(), 
                                  min_shm=numeric(), max_shm=numeric())
}}
cat("  -", ISO1_NAME, "clusters with SHM data:", nrow(cluster_shm_iso1), "\\n")

# Calculate SHM for ISO2 clusters
cat("Calculating SHM for", ISO2_NAME, "clusters (parallel processing)...\\n")

if(IS_WINDOWS) {{
  clusterExport(CL, c("clusters_iso2", "raw_data_iso2", "shm_method_iso2"), envir = environment())
}}

cluster_shm_iso2_list <- par_lapply(1:length(clusters_iso2), function(i) {{
  cluster_seqs <- clusters_iso2[[i]]
  idx <- which(raw_data_iso2$junction_aa %in% cluster_seqs$junction_aa)
  
  if(length(idx) > 0) {{
    shm_values <- get_shm_values(raw_data_iso2, idx, shm_method_iso2)
    
    if(!is.null(shm_values) && length(shm_values) > 0) {{
      return(data.frame(
        cluster_id = i, cluster_size = nrow(cluster_seqs),
        mean_shm = mean(shm_values),
        sd_shm = ifelse(length(shm_values) > 1, sd(shm_values), NA),
        min_shm = min(shm_values), max_shm = max(shm_values)
      ))
    }}
  }}
  return(NULL)
}})

cluster_shm_iso2 <- do.call(rbind, Filter(Negate(is.null), cluster_shm_iso2_list))
if(is.null(cluster_shm_iso2)) {{
  cluster_shm_iso2 <- data.frame(cluster_id=integer(), cluster_size=integer(), 
                                  mean_shm=numeric(), sd_shm=numeric(), 
                                  min_shm=numeric(), max_shm=numeric())
}}
cat("  -", ISO2_NAME, "clusters with SHM data:", nrow(cluster_shm_iso2), "\\n\\n")

write.csv(cluster_shm_iso1, file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_cluster_SHM.csv")), row.names = FALSE)
write.csv(cluster_shm_iso2, file.path(OUTPUT_DIR, paste0(ISO2_NAME, "_cluster_SHM.csv")), row.names = FALSE)

# =============================================================================
# SECTION 6: IDENTIFY SHARED CLUSTERS (OPTIMIZED)
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 6: SHARED CLUSTER IDENTIFICATION\\n")
cat("================================================================================\\n\\n")

cat("Finding shared clusters (V/J match + junction similarity >=", 
    JUNCTION_SIMILARITY_THRESHOLD, "%) - OPTIMIZED...\\n")

iso1_consensus <- lapply(clusters_iso1, get_cluster_consensus_genes)
iso2_consensus <- lapply(clusters_iso2, get_cluster_consensus_genes)

# Create data.tables for fast V/J matching
iso1_vj_dt <- data.table(
  cluster_id = 1:length(clusters_iso1),
  v_gene = sapply(iso1_consensus, `[[`, "v_gene"),
  j_gene = sapply(iso1_consensus, `[[`, "j_gene")
)

iso2_vj_dt <- data.table(
  cluster_id = 1:length(clusters_iso2),
  v_gene = sapply(iso2_consensus, `[[`, "v_gene"),
  j_gene = sapply(iso2_consensus, `[[`, "j_gene")
)

# Fast join on V/J genes - eliminates most comparisons immediately
setkey(iso1_vj_dt, v_gene, j_gene)
setkey(iso2_vj_dt, v_gene, j_gene)
candidate_pairs <- merge(iso1_vj_dt, iso2_vj_dt, by = c("v_gene", "j_gene"),
                         suffixes = c("_iso1", "_iso2"), allow.cartesian = TRUE)

cat("  - V/J matched candidate pairs:", nrow(candidate_pairs), "(reduced from", 
    length(clusters_iso1) * length(clusters_iso2), ")\\n")

# Pre-extract junctions for faster access
iso1_junctions_list <- lapply(clusters_iso1, function(x) x$junction_aa)
iso2_junctions_list <- lapply(clusters_iso2, function(x) x$junction_aa)

# Store all junction similarities for histogram
all_junction_similarities <- c()

# Process candidate pairs in parallel (Windows Compatible)
cat("  - Processing candidate pairs (parallel)...\\n")

# Export required objects to cluster for Windows
if(IS_WINDOWS) {{
  clusterExport(CL, c("candidate_pairs", "iso1_junctions_list", "iso2_junctions_list", 
                      "JUNCTION_SIMILARITY_THRESHOLD"), envir = environment())
  clusterEvalQ(CL, library(stringdist))
}}

n_candidates <- nrow(candidate_pairs)

shared_clusters_list <- par_lapply(1:n_candidates, function(row_idx) {{
  
  i <- candidate_pairs$cluster_id_iso1[row_idx]
  j <- candidate_pairs$cluster_id_iso2[row_idx]
  v_gene <- candidate_pairs$v_gene[row_idx]
  j_gene <- candidate_pairs$j_gene[row_idx]
  
  iso1_juncs <- iso1_junctions_list[[i]]
  iso2_juncs <- iso2_junctions_list[[j]]
  
  # Remove NAs
  iso1_juncs <- iso1_juncs[!is.na(iso1_juncs) & iso1_juncs != ""]
  iso2_juncs <- iso2_juncs[!is.na(iso2_juncs) & iso2_juncs != ""]
  
  if(length(iso1_juncs) == 0 || length(iso2_juncs) == 0) {{
    return(list(result = NULL, all_sims = NULL))
  }}
  
  # Compute distance matrix in one vectorized call
  dist_mat <- stringdist::stringdistmatrix(iso1_juncs, iso2_juncs, method = "lv")
  max_len_mat <- outer(nchar(iso1_juncs), nchar(iso2_juncs), pmax)
  sim_mat <- (1 - dist_mat / max_len_mat) * 100
  
  matching_mask <- sim_mat >= JUNCTION_SIMILARITY_THRESHOLD
  matching_pairs <- sum(matching_mask)
  
  if(matching_pairs > 0) {{
    similarities <- sim_mat[matching_mask]
    return(list(
      result = data.frame(
        iso1_cluster = i, iso2_cluster = j,
        iso1_size = length(iso1_juncs), iso2_size = length(iso2_juncs),
        iso1_v_gene = v_gene, iso1_j_gene = j_gene,
        matching_pairs = matching_pairs,
        avg_junction_similarity = mean(similarities),
        max_junction_similarity = max(similarities),
        stringsAsFactors = FALSE
      ),
      all_sims = as.vector(sim_mat)
    ))
  }}
  return(list(result = NULL, all_sims = as.vector(sim_mat)))
}})

cat("  - Combining results...\\n")

# Combine results
shared_clusters <- do.call(rbind, lapply(shared_clusters_list, function(x) x$result))
all_junction_similarities <- unlist(lapply(shared_clusters_list, function(x) x$all_sims))

if(is.null(shared_clusters) || nrow(shared_clusters) == 0) {{
  shared_clusters <- data.frame(
    iso1_cluster = integer(), iso2_cluster = integer(),
    iso1_size = integer(), iso2_size = integer(),
    iso1_v_gene = character(), iso1_j_gene = character(),
    matching_pairs = integer(), avg_junction_similarity = numeric(),
    max_junction_similarity = numeric(), stringsAsFactors = FALSE
  )
}}

shared_clusters <- shared_clusters[order(-shared_clusters$matching_pairs), ]

cat("\\n  - Shared cluster pairs found:", nrow(shared_clusters), "\\n")
cat("  - Total matching sequence pairs:", sum(shared_clusters$matching_pairs), "\\n\\n")

unique_iso1_in_shared <- unique(shared_clusters$iso1_cluster)
unique_iso2_in_shared <- unique(shared_clusters$iso2_cluster)

write.csv(shared_clusters, file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_", ISO2_NAME, "_shared_clusters.csv")), row.names = FALSE)

# =============================================================================
# SECTION 7: JUNCTION SIMILARITY DISTRIBUTION
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 7: JUNCTION SIMILARITY DISTRIBUTION\\n")
cat("================================================================================\\n\\n")

cat("Creating junction similarity histogram...\\n")

matching_similarities <- all_junction_similarities[all_junction_similarities >= JUNCTION_SIMILARITY_THRESHOLD]

if(length(matching_similarities) > 0) {{
  mean_sim <- mean(matching_similarities)
  
  p_junction_hist <- ggplot(data.frame(similarity = matching_similarities), aes(x = similarity)) +
    geom_histogram(binwidth = 5, fill = "coral", color = "white", alpha = 0.8) +
    geom_vline(xintercept = JUNCTION_SIMILARITY_THRESHOLD, linetype = "dashed", color = "darkblue", size = 1.2) +
    labs(title = "Distribution of Junction Similarity in Matching Pairs",
         subtitle = paste0("n = ", length(matching_similarities), " matching pairs (>=", 
                           JUNCTION_SIMILARITY_THRESHOLD, "% similarity)"),
         x = "Junction Similarity (%)",
         y = "Number of Matching Pairs") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11))
  
  ggsave(file.path(OUTPUT_DIR, "Junction_similarity_distribution.pdf"), 
         plot = p_junction_hist, width = 10, height = 6)
  ggsave(file.path(OUTPUT_DIR, "Junction_similarity_distribution.png"), 
         plot = p_junction_hist, width = 10, height = 6, dpi = 150)
}}

# =============================================================================
# SECTION 8: MULTI-PANEL SHARED CLUSTER ANALYSIS
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 8: MULTI-PANEL SHARED CLUSTER ANALYSIS\\n")
cat("================================================================================\\n\\n")

cat("Creating multi-panel shared cluster visualization...\\n")

if(nrow(shared_clusters) > 0) {{
  
  # Panel 1: ISO1 vs ISO2 cluster sizes
  p1 <- ggplot(shared_clusters, aes(x = iso1_size, y = iso2_size)) +
    geom_point(aes(size = matching_pairs, color = avg_junction_similarity), alpha = 0.7) +
    scale_color_gradient(low = "gold", high = "red", name = "Avg Junction\\nSimilarity (%)") +
    scale_size_continuous(range = c(2, 15), name = "Matching\\nPairs") +
    labs(title = paste0(ISO1_NAME, " vs ", ISO2_NAME, " Cluster Sizes (Shared Clusters)"),
         x = paste0(ISO1_NAME, " Cluster Size"), y = paste0(ISO2_NAME, " Cluster Size")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
  
  # Panel 2: Distribution of average junction similarity
  p2 <- ggplot(shared_clusters, aes(x = avg_junction_similarity)) +
    geom_histogram(binwidth = 2, fill = "steelblue", color = "white") +
    geom_vline(xintercept = JUNCTION_SIMILARITY_THRESHOLD, linetype = "dashed", 
               color = "red", size = 1) +
    labs(title = "Distribution of Average Junction Similarity in Shared Clusters",
         x = "Average Junction Similarity (%)", y = "Number of Cluster Pairs") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
  
  # Panel 3: Matching pairs vs junction similarity
  p3 <- ggplot(shared_clusters, aes(x = matching_pairs, y = avg_junction_similarity)) +
    geom_point(aes(color = iso1_size + iso2_size), alpha = 0.6, size = 3) +
    scale_color_gradient(low = "lightblue", high = "darkblue", name = "Combined\\nCluster Size") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
    labs(title = "Matching Pairs vs Junction Similarity",
         x = "Number of Matching Pairs", y = "Average Junction Similarity (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
  
  # Panel 4: Cluster size distribution violin
  size_data <- data.frame(
    size = c(shared_clusters$iso1_size, shared_clusters$iso2_size),
    isotype = rep(c(ISO1_NAME, ISO2_NAME), each = nrow(shared_clusters))
  )
  
  p4 <- ggplot(size_data, aes(x = isotype, y = size, fill = isotype)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", outlier.color = "red") +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1, color = "red") +
    scale_fill_manual(values = c("coral", "steelblue")) +
    labs(title = "Cluster Size Distribution in Shared Clusters",
         x = "Isotype", y = "Cluster Size") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "none")
  
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                                top = paste0(ISO1_NAME, "-", ISO2_NAME, " Shared Clusters Analysis (", 
                                             JUNCTION_SIMILARITY_THRESHOLD, "% Junction Similarity)"))
  
  ggsave(file.path(OUTPUT_DIR, "Shared_clusters_multi_panel.pdf"), 
         plot = combined_plot, width = 14, height = 12)
  ggsave(file.path(OUTPUT_DIR, "Shared_clusters_multi_panel.png"), 
         plot = combined_plot, width = 14, height = 12, dpi = 150)
}}

# =============================================================================
# SECTION 9: V/J GENE COMBINATIONS IN SHARED CLUSTERS
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 9: V/J GENE COMBINATIONS IN SHARED CLUSTERS\\n")
cat("================================================================================\\n\\n")

cat("Creating V/J gene combination plot...\\n")

if(nrow(shared_clusters) > 0) {{
  vj_usage <- as.data.frame(table(paste(shared_clusters$iso1_v_gene, 
                                        shared_clusters$iso1_j_gene, sep = " + ")))
  colnames(vj_usage) <- c("VJ_combination", "Frequency")
  vj_usage <- vj_usage[order(-vj_usage$Frequency), ]
  vj_usage_top <- head(vj_usage, 20)
  
  p_vj <- ggplot(vj_usage_top, aes(x = reorder(VJ_combination, Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Top 20 V/J Gene Combinations in Shared Clusters",
         x = "V + J Gene Combination",
         y = "Number of Shared Cluster Pairs") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y = element_text(size = 8))
  
  ggsave(file.path(OUTPUT_DIR, "VJ_combinations_shared_clusters.pdf"), 
         plot = p_vj, width = 12, height = 10)

  ggsave(file.path(OUTPUT_DIR, "VJ_combinations_shared_clusters.png"), 
             plot = p_vj, width = 12, height = 10, dpi = 150)
}}

# =============================================================================
# OVERLAP PIE CHARTS
# =============================================================================

cat("Creating overlap pie charts...\\n")

total_iso1_clusters <- length(clusters_iso1)
total_iso2_clusters <- length(clusters_iso2)

iso1_with_iso2 <- length(unique_iso1_in_shared)
iso1_not_with_iso2 <- total_iso1_clusters - iso1_with_iso2

iso2_with_iso1 <- length(unique_iso2_in_shared)
iso2_not_with_iso1 <- total_iso2_clusters - iso2_with_iso1

iso1_overlap_data <- data.frame(
  category = c(paste0(ISO1_NAME, " clusters shared with ", ISO2_NAME), paste0(ISO1_NAME, " clusters not in ", ISO2_NAME)),
  count = c(iso1_with_iso2, iso1_not_with_iso2)
)
iso1_overlap_data$percentage <- round(iso1_overlap_data$count / total_iso1_clusters * 100, 1)
iso1_overlap_data$label <- paste0(iso1_overlap_data$count, "\\n(", iso1_overlap_data$percentage, "%)")

iso2_overlap_data <- data.frame(
  category = c(paste0(ISO2_NAME, " clusters shared with ", ISO1_NAME), paste0(ISO2_NAME, " clusters not in ", ISO1_NAME)),
  count = c(iso2_with_iso1, iso2_not_with_iso1)
)
iso2_overlap_data$percentage <- round(iso2_overlap_data$count / total_iso2_clusters * 100, 1)
iso2_overlap_data$label <- paste0(iso2_overlap_data$count, "\\n(", iso2_overlap_data$percentage, "%)")

p_iso1_pie <- ggplot(iso1_overlap_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("coral", "gray70")) +
  labs(title = paste0(ISO1_NAME, " Clusters: Overlap with ", ISO2_NAME),
       subtitle = paste0("Total ", ISO1_NAME, " clusters: ", total_iso1_clusters)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

p_iso2_pie <- ggplot(iso2_overlap_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 5, fontface = "bold") +
  scale_fill_manual(values = c("steelblue", "gray70")) +
  labs(title = paste0(ISO2_NAME, " Clusters: Overlap with ", ISO1_NAME),
       subtitle = paste0("Total ", ISO2_NAME, " clusters: ", total_iso2_clusters)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

combined_pie <- grid.arrange(p_iso1_pie, p_iso2_pie, ncol = 2)

ggsave(file.path(OUTPUT_DIR, "Cluster_overlap_pie_charts.pdf"), 
       plot = combined_pie, width = 12, height = 6)

ggsave(file.path(OUTPUT_DIR, "Cluster_overlap_pie_charts.png"), 
       plot = combined_pie, width = 12, height = 6, dpi = 150)

cat("  SAVED: Cluster_overlap_pie_charts.pdf\\n")

cat("\\nOVERLAP SUMMARY:\\n")
cat("  ", ISO1_NAME, "clusters with", ISO2_NAME, ":", iso1_with_iso2, "/", total_iso1_clusters, 
    "(", round(iso1_with_iso2/total_iso1_clusters*100, 1), "%)\\n")
cat("  ", ISO2_NAME, "clusters with", ISO1_NAME, ":", iso2_with_iso1, "/", total_iso2_clusters, 
    "(", round(iso2_with_iso1/total_iso2_clusters*100, 1), "%)\\n")

# =============================================================================
# SECTION 10: SHM FOR SHARED CLUSTERS WITH MUTATION LEVELS
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 10: SHM ANALYSIS FOR SHARED CLUSTERS\\n")
cat("================================================================================\\n\\n")

shared_iso1_shm <- data.frame(
  cluster_id = integer(), cluster_size = integer(),
  mean_shm = numeric(), sd_shm = numeric(),
  isotype = character(), stringsAsFactors = FALSE
)

cat("Using SHM method:", shm_method_iso1, "for shared cluster analysis\\n\\n")

for(i in unique_iso1_in_shared) {{
  cluster_seqs <- clusters_iso1[[i]]
  idx <- which(raw_data_iso1$junction_aa %in% cluster_seqs$junction_aa)
  
  if(length(idx) > 0) {{
    # Use the unified get_shm_values function from Section 5
    # This automatically handles both IMGT (alignment) and MiXCR (mutation_rate) formats
    shm_values <- get_shm_values(raw_data_iso1, idx, shm_method_iso1)
    
    if(!is.null(shm_values) && length(shm_values) > 0) {{
      shared_iso1_shm <- rbind(shared_iso1_shm, data.frame(
        cluster_id = i, cluster_size = nrow(cluster_seqs),
        mean_shm = mean(shm_values),
        sd_shm = ifelse(length(shm_values) > 1, sd(shm_values), 0),
        isotype = ISO1_NAME, stringsAsFactors = FALSE
      ))
    }}
  }}
}}

# Define mutation level thresholds
if(nrow(shared_iso1_shm) > 0) {{
  low_threshold <- quantile(shared_iso1_shm$mean_shm, 0.33, na.rm = TRUE)
  high_threshold <- quantile(shared_iso1_shm$mean_shm, 0.67, na.rm = TRUE)
  
  cat("SHM Thresholds:\\n")
  cat("  Low threshold:", round(low_threshold, 2), "%\\n")
  cat("  High threshold:", round(high_threshold, 2), "%\\n\\n")
  
  shared_iso1_shm$mutation_level <- cut(
    shared_iso1_shm$mean_shm,
    breaks = c(-Inf, low_threshold, high_threshold, Inf),
    labels = c("Low", "Moderate", "High")
  )
  
  shared_iso1_shm <- shared_iso1_shm[order(shared_iso1_shm$mean_shm), ]
  
  cat("Shared", ISO1_NAME, "Clusters by Mutation Level:\\n")
  print(table(shared_iso1_shm$mutation_level))
  cat("\\n")
  
  write.csv(shared_iso1_shm, file.path(OUTPUT_DIR, paste0("Shared_", ISO1_NAME, "_clusters_SHM.csv")), row.names = FALSE)

  # =============================================================================
  # SHM PER SHARED CLUSTER VISUALIZATION (Ascending Order)
  # =============================================================================
  cat("Creating SHM per shared cluster visualization...\\n")
  
  # Order by mean_shm for x-axis
  shared_iso1_shm_ordered <- shared_iso1_shm[order(shared_iso1_shm$mean_shm), ]
  shared_iso1_shm_ordered$order_idx <- 1:nrow(shared_iso1_shm_ordered)
  
  p_shm_ordered <- ggplot(shared_iso1_shm_ordered, 
                          aes(x = factor(order_idx), y = mean_shm)) +
    geom_point(aes(size = cluster_size, color = mean_shm), alpha = 0.8) +
    scale_color_gradient(low = "steelblue", high = "orange", name = "Mean SHM (%)") +
    scale_size_continuous(range = c(2, 12), name = "Cluster Size") +
    geom_hline(yintercept = low_threshold, linetype = "dashed", color = "green4", size = 1) +
    geom_hline(yintercept = high_threshold, linetype = "dashed", color = "darkred", size = 1) +
    annotate("text", x = 1, y = low_threshold + 0.15, label = "Low", 
             color = "green4", hjust = 0, fontface = "bold", size = 4) +
    annotate("text", x = 1, y = high_threshold + 0.15, label = "High", 
             color = "darkred", hjust = 0, fontface = "bold", size = 4) +
    labs(title = paste0("Mean SHM per Shared ", ISO1_NAME, " Cluster (Ascending Order)"),
         subtitle = paste0("n = ", nrow(shared_iso1_shm_ordered), " clusters | Method: ", shm_method_iso1),
         x = "Cluster ID (Ordered by SHM)",
         y = "Mean SHM (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "right")
  
  ggsave(file.path(OUTPUT_DIR, paste0("Shared_", ISO1_NAME, "_SHM_ascending.pdf")), 
         plot = p_shm_ordered, width = 12, height = 6)
  ggsave(file.path(OUTPUT_DIR, paste0("Shared_", ISO1_NAME, "_SHM_ascending.png")), 
         plot = p_shm_ordered, width = 12, height = 6, dpi = 150)
  
  cat("  SAVED: Shared_", ISO1_NAME, "_SHM_ascending.pdf and .png\\n")

}} else {{
  low_threshold <- 0
  high_threshold <- 0
  cat("Warning: No SHM data available for shared clusters\\n")
  cat("  Possible reasons:\\n")
  cat("  - No matching junction_aa between clusters and raw data\\n")
  cat("  - Missing mutation data columns\\n")
  cat("  Available columns:", paste(head(colnames(raw_data_iso1), 20), collapse = ", "), "...\\n")
}}

# =============================================================================
# SECTION 11: DIVERGENT CLONE IDENTIFICATION (OPTIMIZED)
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 11: DIVERGENT", ISO1_NAME, "CLONE IDENTIFICATION\\n")
cat("================================================================================\\n\\n")

cat("Identifying divergent", ISO1_NAME, "clones (>=", DIVERGENT_MIN_DIFFERENCES, 
    " AA differences from all", ISO2_NAME, ") - OPTIMIZED...\\n")
cat("Analyzing ALL mutation levels (Low, Moderate, High)...\\n\\n")

all_iso2_junctions <- unique(raw_data_iso2$junction_aa)
all_iso2_junctions <- all_iso2_junctions[!is.na(all_iso2_junctions) & all_iso2_junctions != ""]
cat("  - Pre-indexed", length(all_iso2_junctions), ISO2_NAME, "junctions for comparison\\n")

# Initialize empty results dataframe
divergent_results <- data.frame(
  iso1_cluster_id = integer(), iso2_cluster_id = integer(),
  mutation_level = character(), iso1_cluster_size = integer(),
  num_divergent = integer(), pct_divergent = numeric(),
  mean_shm = numeric(), divergent_junctions = character(),
  v_gene = character(), j_gene = character(),
  stringsAsFactors = FALSE
)

# Optimized batch divergent checking function
check_divergent_batch <- function(query_junctions, ref_junctions, min_diff) {{
  if(length(query_junctions) == 0) return(character(0))
  
  # Use stringdistmatrix for vectorized comparison
  dist_mat <- stringdist::stringdistmatrix(query_junctions, ref_junctions, method = "lv")
  min_dists <- apply(dist_mat, 1, min, na.rm = TRUE)
  
  return(query_junctions[min_dists >= min_diff])
}}

if(nrow(shared_iso1_shm) > 0) {{
  cat("  - Processing", nrow(shared_iso1_shm), "shared clusters (parallel)...\\n")
  
  # Export required objects to cluster for Windows
  if(IS_WINDOWS) {{
    clusterExport(CL, c("shared_iso1_shm", "clusters_iso1", "shared_clusters", 
                        "all_iso2_junctions", "DIVERGENT_MIN_DIFFERENCES",
                        "extract_gene", "check_divergent_batch"), envir = environment())
    clusterEvalQ(CL, library(stringdist))
  }}
  
  n_shared <- nrow(shared_iso1_shm)
  
  divergent_results_list <- par_lapply(1:n_shared, function(row_i) {{
    
    iso1_cluster_id <- shared_iso1_shm$cluster_id[row_i]
    mutation_level <- as.character(shared_iso1_shm$mutation_level[row_i])
    mean_shm <- shared_iso1_shm$mean_shm[row_i]
    
    iso1_cluster <- clusters_iso1[[iso1_cluster_id]]
    iso1_junctions <- unique(iso1_cluster$junction_aa)
    iso1_junctions <- iso1_junctions[!is.na(iso1_junctions) & iso1_junctions != ""]
    
    v_gene <- extract_gene(iso1_cluster$v_call[1])
    j_gene <- extract_gene(iso1_cluster$j_call[1])
    
    matching_rows <- which(shared_clusters$iso1_cluster == iso1_cluster_id)
    if(length(matching_rows) > 0) {{
      iso2_cluster_id <- shared_clusters$iso2_cluster[matching_rows[1]]
    }} else {{
      iso2_cluster_id <- NA
    }}
    
    # Batch check for divergent sequences
    divergent_seqs <- check_divergent_batch(iso1_junctions, all_iso2_junctions, DIVERGENT_MIN_DIFFERENCES)
    
    if(length(divergent_seqs) > 0) {{
      return(data.frame(
        iso1_cluster_id = iso1_cluster_id,
        iso2_cluster_id = iso2_cluster_id,
        mutation_level = mutation_level,
        iso1_cluster_size = nrow(iso1_cluster),
        num_divergent = length(divergent_seqs),
        pct_divergent = (length(divergent_seqs) / nrow(iso1_cluster)) * 100,
        mean_shm = mean_shm,
        divergent_junctions = paste(divergent_seqs, collapse = "; "),
        v_gene = v_gene,
        j_gene = j_gene,
        stringsAsFactors = FALSE
      ))
    }}
    
    return(NULL)
  }})
  
  cat("  - Combining divergent results...\\n")
  
  divergent_results <- do.call(rbind, Filter(Negate(is.null), divergent_results_list))
  
  if(is.null(divergent_results) || nrow(divergent_results) == 0) {{
    divergent_results <- data.frame(
      iso1_cluster_id = integer(), iso2_cluster_id = integer(),
      mutation_level = character(), iso1_cluster_size = integer(),
      num_divergent = integer(), pct_divergent = numeric(),
      mean_shm = numeric(), divergent_junctions = character(),
      v_gene = character(), j_gene = character(),
      stringsAsFactors = FALSE
    )
  }}
}}

cat("\\n")
cat("DIVERGENT", ISO1_NAME, "SUMMARY BY MUTATION LEVEL:\\n")
cat("=========================================\\n")

if(nrow(divergent_results) > 0) {{
  divergent_summary <- divergent_results %>%
    group_by(mutation_level) %>%
    summarise(
      n_clusters = n(),
      total_divergent = sum(num_divergent),
      mean_pct_divergent = round(mean(pct_divergent), 1),
      mean_shm = round(mean(mean_shm), 2),
      .groups = "drop"
    )
  
  print(divergent_summary)
  cat("\\n")
  
  cat("Divergent", ISO1_NAME, "in LOW mutation clusters:", 
      sum(divergent_results$num_divergent[divergent_results$mutation_level == "Low"]), "\\n")
  cat("Divergent", ISO1_NAME, "in MODERATE mutation clusters:", 
      sum(divergent_results$num_divergent[divergent_results$mutation_level == "Moderate"]), "\\n")
  cat("Divergent", ISO1_NAME, "in HIGH mutation clusters:", 
      sum(divergent_results$num_divergent[divergent_results$mutation_level == "High"]), "\\n\\n")
}} else {{
  cat("No divergent", ISO1_NAME, "clones found.\\n\\n")
}}
# =============================================================================
# NON-DIVERGENT CLONE ANALYSIS BY MUTATION LEVEL
# =============================================================================

cat("\\n")
cat("NON-DIVERGENT", ISO1_NAME, "SUMMARY BY MUTATION LEVEL:\\n")
cat("=========================================\\n")

# Initialize non-divergent results dataframe
nondivergent_results <- data.frame(
  mutation_level = character(),
  total_clusters = integer(),
  total_sequences = integer(),
  divergent_clusters = integer(),
  divergent_sequences = integer(),
  nondivergent_clusters = integer(),
  nondivergent_sequences = integer(),
  pct_nondivergent_clusters = numeric(),
  pct_nondivergent_sequences = numeric(),
  stringsAsFactors = FALSE
)

# Calculate for each mutation level
if(nrow(shared_iso1_shm) > 0) {{
  for(mut_level in c("Low", "Moderate", "High")) {{
    
    # Get clusters in this mutation level
    level_clusters <- shared_iso1_shm$cluster_id[shared_iso1_shm$mutation_level == mut_level]
    n_level_clusters <- length(level_clusters)
    
    # Total sequences in this mutation level
    total_seqs_level <- sum(sapply(level_clusters, function(id) nrow(clusters_iso1[[id]])))
    
    # Divergent clusters and sequences in this level
    div_in_level <- divergent_results[divergent_results$mutation_level == mut_level, ]
    n_div_clusters <- nrow(div_in_level)
    n_div_seqs <- if(n_div_clusters > 0) sum(div_in_level$num_divergent) else 0
    
    # Non-divergent clusters (clusters with 0 divergent sequences)
    div_cluster_ids <- if(n_div_clusters > 0) div_in_level$iso1_cluster_id else c()
    nondiv_cluster_ids <- setdiff(level_clusters, div_cluster_ids)
    n_nondiv_clusters <- length(nondiv_cluster_ids)
    
    # Non-divergent sequences = total - divergent
    n_nondiv_seqs <- total_seqs_level - n_div_seqs
    
    # Percentages
    pct_nondiv_clusters <- if(n_level_clusters > 0) round(n_nondiv_clusters/n_level_clusters*100, 1) else 0
    pct_nondiv_seqs <- if(total_seqs_level > 0) round(n_nondiv_seqs/total_seqs_level*100, 1) else 0
    
    nondivergent_results <- rbind(nondivergent_results, data.frame(
      mutation_level = mut_level,
      total_clusters = n_level_clusters,
      total_sequences = total_seqs_level,
      divergent_clusters = n_div_clusters,
      divergent_sequences = n_div_seqs,
      nondivergent_clusters = n_nondiv_clusters,
      nondivergent_sequences = n_nondiv_seqs,
      pct_nondivergent_clusters = pct_nondiv_clusters,
      pct_nondivergent_sequences = pct_nondiv_seqs,
      stringsAsFactors = FALSE
    ))
    
    cat(mut_level, "mutation level:\\n")
    cat("  Total clusters:", n_level_clusters, "| Total sequences:", total_seqs_level, "\\n")
    cat("  Divergent clusters:", n_div_clusters, "| Divergent sequences:", n_div_seqs, "\\n")
    cat("  NON-divergent clusters:", n_nondiv_clusters, "(", pct_nondiv_clusters, "%)\\n")
    cat("  NON-divergent sequences:", n_nondiv_seqs, "(", pct_nondiv_seqs, "%)\\n\\n")
  }}
  
  # Calculate totals
  total_iso1_in_shared <- sum(nondivergent_results$total_sequences)
  total_divergent_seqs <- sum(nondivergent_results$divergent_sequences)
  total_nondivergent_seqs <- sum(nondivergent_results$nondivergent_sequences)
  total_nondiv_clusters <- sum(nondivergent_results$nondivergent_clusters)
  total_shared_clusters <- sum(nondivergent_results$total_clusters)
  
  cat("OVERALL TOTALS:\\n")
  cat("  Total", ISO1_NAME, "sequences in shared clusters:", total_iso1_in_shared, "\\n")
  cat("  Total divergent", ISO1_NAME, "sequences:", total_divergent_seqs, "\\n")
  cat("  Total NON-divergent", ISO1_NAME, "sequences:", total_nondivergent_seqs, 
      "(", round(total_nondivergent_seqs/total_iso1_in_shared*100, 1), "%)\\n")
  cat("  Total NON-divergent", ISO1_NAME, "clusters:", total_nondiv_clusters,
      "(", round(total_nondiv_clusters/total_shared_clusters*100, 1), "%)\\n\\n")
  
}} else {{
  total_iso1_in_shared <- 0
  total_divergent_seqs <- 0
  total_nondivergent_seqs <- 0
  total_nondiv_clusters <- 0
  total_shared_clusters <- 0
  cat("No shared clusters with SHM data available.\\n\\n")
}}

# Save non-divergent results
write.csv(nondivergent_results, file.path(OUTPUT_DIR, paste0("NonDivergent_", ISO1_NAME, "_by_mutation_level.csv")), row.names = FALSE)
cat("  SAVED: NonDivergent_", ISO1_NAME, "_by_mutation_level.csv\\n\\n")

write.csv(divergent_results, file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_clones_all_levels.csv")), row.names = FALSE)

# Filter: Remove clusters with only 1 divergent sequence
cat("\\n================================================================================\\n")
cat("[EXCLUSION] SINGLE-COPY DIVERGENT", ISO1_NAME, "FILTERING:\\n")
cat("[EXCLUSION] Clusters before filtering:", nrow(divergent_results), "\\n")

# Count and store excluded single-copy divergent clusters
excluded_single_copy <- divergent_results[divergent_results$num_divergent == 1, ]
n_excluded_clusters <- nrow(excluded_single_copy)
n_excluded_sequences <- sum(excluded_single_copy$num_divergent)

divergent_results <- divergent_results[divergent_results$num_divergent >= 2, ]

cat("[EXCLUSION] Clusters after filtering:", nrow(divergent_results), "\\n")
cat("[EXCLUSION] Single-copy divergent", ISO1_NAME, "clusters EXCLUDED:", n_excluded_clusters, "\\n")
cat("[EXCLUSION] Single-copy divergent", ISO1_NAME, "sequences EXCLUDED:", n_excluded_sequences, "\\n")

if(n_excluded_clusters > 0) {{
  cat("[EXCLUSION] Excluded cluster IDs:", paste(excluded_single_copy$iso1_cluster_id, collapse = ", "), "\\n")
}}
cat("================================================================================\\n\\n")

write.csv(divergent_results, file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_clones_filtered.csv")), row.names = FALSE)

# =============================================================================
# SECTION 12: DIVERGENT CLONE VISUALIZATIONS
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 12: DIVERGENT", ISO1_NAME, "VISUALIZATIONS\\n")
cat("================================================================================\\n\\n")

if(nrow(divergent_results) > 0) {{
  
  # Bar plot by mutation level
  p_div_bar <- ggplot(divergent_results, aes(x = mutation_level, fill = mutation_level)) +
    geom_bar() +
    scale_fill_manual(values = c("Low" = "green3", "Moderate" = "gold", "High" = "red2")) +
    geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 6, fontface = "bold") +
    labs(title = paste0("Clusters with Divergent ", ISO1_NAME, " by Mutation Level"),
         subtitle = paste0("Divergent = >=", DIVERGENT_MIN_DIFFERENCES, " AA differences from all ", ISO2_NAME),
         x = "Mutation Level", y = "Number of Clusters") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none")
  
  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_by_mutation_level.pdf")), plot = p_div_bar, width = 10, height = 6)

  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_by_mutation_level.png")), plot = p_div_bar, width = 10, height = 6, dpi = 150)
  
  
  # Total divergent sequences by mutation level
  total_by_level <- divergent_results %>%
    group_by(mutation_level) %>%
    summarise(total = sum(num_divergent), .groups = "drop")
  
  total_by_level$mutation_level <- factor(total_by_level$mutation_level, levels = c("Low", "Moderate", "High"))
  
  p_total_div <- ggplot(total_by_level, aes(x = mutation_level, y = total, fill = mutation_level)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Low" = "green3", "Moderate" = "gold", "High" = "red2")) +
    geom_text(aes(label = total), vjust = -0.5, size = 6, fontface = "bold") +
    labs(title = paste0("Total Divergent ", ISO1_NAME, " Sequences by Mutation Level"),
         x = "Mutation Level", y = "Total Divergent Sequences") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none")
  
  ggsave(file.path(OUTPUT_DIR, paste0("Total_divergent_", ISO1_NAME, "_by_mutation_level.pdf")), plot = p_total_div, width = 10, height = 6)

  ggsave(file.path(OUTPUT_DIR, paste0("Total_divergent_", ISO1_NAME, "_by_mutation_level.png")), plot = p_total_div, width = 10, height = 6, dpi = 150)
  
  # SHM vs % Divergent
  p_shm_div <- ggplot(divergent_results, aes(x = mean_shm, y = pct_divergent, 
                                             color = mutation_level, size = num_divergent)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("Low" = "green3", "Moderate" = "gold2", "High" = "red2"),
                       name = "Mutation Level") +
    scale_size_continuous(range = c(3, 15), name = "# Divergent") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    labs(title = paste0("Mean SHM vs % Divergent ", ISO1_NAME, " Sequences"),
         x = "Mean SHM (%)", y = paste0("% Divergent ", ISO1_NAME)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_SHM_vs_pct.pdf")), plot = p_shm_div, width = 10, height = 7)
  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_SHM_vs_pct.png")), plot = p_shm_div, width = 10, height = 7, dpi = 150)
}}

# =============================================================================
# SECTION 12B: SHARED CLUSTERS HEATMAP
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 12B: SHARED CLUSTERS HEATMAP\\n")
cat("================================================================================\\n\\n")

cat("Creating shared clusters heatmap...\\n")

if(nrow(shared_clusters) > 0) {{
  
  p_heatmap <- ggplot(shared_clusters, aes(x = factor(iso2_cluster), 
                                           y = factor(iso1_cluster), 
                                           fill = matching_pairs)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient(low = "lightblue", high = "darkred", 
                        name = "Matching\\nPairs") +
    labs(title = paste0("Shared Clusters: ", ISO1_NAME, " vs ", ISO2_NAME, "\\n(V/J match + Junction >=", 
                        JUNCTION_SIMILARITY_THRESHOLD, "% similar)"),
         x = paste0(ISO2_NAME, " Cluster ID"),
         y = paste0(ISO1_NAME, " Cluster ID")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "right")
  
  ggsave(file.path(OUTPUT_DIR, "Shared_clusters_heatmap.pdf"), 
         plot = p_heatmap, width = 14, height = 10)
         
  ggsave(file.path(OUTPUT_DIR, "Shared_clusters_heatmap.png"), 
         plot = p_heatmap, width = 14, height = 10, dpi = 150)
  
  cat("  SAVED: Shared_clusters_heatmap.pdf\\n")
}}

# =============================================================================
# SECTION 12C: DIVERGENT CLONE COPY SIZE VISUALIZATION
# =============================================================================

cat("\\n================================================================================\\n")
cat("SECTION 12C: DIVERGENT", ISO1_NAME, "COPY SIZE ANALYSIS\\n")
cat("================================================================================\\n\\n")

if(nrow(divergent_results) > 0) {{
  
  # Add ISO2 cluster size to divergent_results
  divergent_results$iso2_cluster_size <- sapply(divergent_results$iso2_cluster_id, function(id) {{
    if(is.na(id)) return(NA)
    nrow(clusters_iso2[[id]])
  }})
  
  # Create long format for grouped bar plot
  divergent_copy_data <- divergent_results %>%
    select(iso1_cluster_id, mutation_level, iso1_cluster_size, 
           num_divergent, iso2_cluster_size) %>%
    tidyr::pivot_longer(
      cols = c(iso1_cluster_size, num_divergent, iso2_cluster_size),
      names_to = "Size_Type",
      values_to = "Count"
    )
  
  divergent_copy_data$Size_Type <- recode(divergent_copy_data$Size_Type,
                                          "iso1_cluster_size" = paste0(ISO1_NAME, " Cluster Size"),
                                          "num_divergent" = paste0("Divergent ", ISO1_NAME),
                                          "iso2_cluster_size" = paste0(ISO2_NAME, " Cluster Size"))
  
  divergent_copy_data$iso1_cluster_id <- factor(divergent_copy_data$iso1_cluster_id)
  divergent_copy_data$Size_Type <- factor(divergent_copy_data$Size_Type, 
                                          levels = c(paste0(ISO1_NAME, " Cluster Size"), 
                                                     paste0("Divergent ", ISO1_NAME), 
                                                     paste0(ISO2_NAME, " Cluster Size")))
  divergent_copy_data$mutation_level <- factor(divergent_copy_data$mutation_level,
                                               levels = c("Low", "Moderate", "High"))
  
  p_copy_sizes <- ggplot(divergent_copy_data, 
                         aes(x = iso1_cluster_id, y = Count, fill = Size_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = Count), 
              position = position_dodge(width = 0.8), 
              vjust = -0.3, size = 3) +
    facet_wrap(~mutation_level, scales = "free_x", ncol = 3) +
    scale_fill_manual(values = c("coral", "red3", "steelblue"),
                      name = "Copy Size Type") +
    labs(title = paste0("Divergent ", ISO1_NAME, ": Cluster Copy Sizes by Mutation Level"),
         subtitle = paste0(ISO1_NAME, " cluster size, number of divergent ", ISO1_NAME, " sequences, and clonally related ", ISO2_NAME, " cluster size"),
         x = paste0(ISO1_NAME, " Cluster ID"),
         y = "Number of Sequences") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 12))
  
  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_copy_sizes.pdf")), 
         plot = p_copy_sizes, width = 14, height = 8)

  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_copy_sizes.png")), 
         plot = p_copy_sizes, width = 14, height = 8, dpi = 150)
  
  cat("  SAVED: Divergent_", ISO1_NAME, "_copy_sizes.pdf\\n")
  
  # Scatter plot comparing ISO1 vs ISO2 cluster sizes
  p_scatter <- ggplot(divergent_results, aes(x = iso1_cluster_size, y = iso2_cluster_size)) +
    geom_point(aes(size = num_divergent, color = mutation_level), alpha = 0.7) +
    geom_text(aes(label = iso1_cluster_id), hjust = -0.3, vjust = -0.3, size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Low" = "green3", "Moderate" = "gold2", "High" = "red2"),
                       name = "Mutation Level") +
    scale_size_continuous(range = c(3, 15), name = paste0("# Divergent ", ISO1_NAME)) +
    labs(title = paste0(ISO1_NAME, " vs ", ISO2_NAME, " Cluster Sizes (Clusters with Divergent ", ISO1_NAME, ")"),
         subtitle = paste0("Point size = number of divergent ", ISO1_NAME, "; dashed line = equal sizes"),
         x = paste0(ISO1_NAME, " Cluster Size (copy number)"),
         y = paste0(ISO2_NAME, " Cluster Size (copy number)")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11))
  
  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_vs_", ISO2_NAME, "_sizes.pdf")), 
         plot = p_scatter, width = 10, height = 8)

  ggsave(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_vs_", ISO2_NAME, "_sizes.png")), 
         plot = p_scatter, width = 10, height = 8, dpi = 150)
  
  cat("  SAVED: Divergent_", ISO1_NAME, "_vs_", ISO2_NAME, "_sizes.pdf\\n")
  
  # Summary table
  summary_table <- divergent_results %>%
    select(iso1_cluster_id, mutation_level, iso1_cluster_size, 
           num_divergent, pct_divergent, iso2_cluster_size, v_gene, j_gene) %>%
    arrange(mutation_level, desc(num_divergent)) %>%
    mutate(pct_divergent = round(pct_divergent, 1))
  
  colnames(summary_table) <- c(paste0(ISO1_NAME, " ID"), "Mutation", paste0(ISO1_NAME, " Size"), "# Divergent", 
                               "% Divergent", paste0(ISO2_NAME, " Size"), "V Gene", "J Gene")
  
  table_plot <- tableGrob(
    summary_table,
    rows = NULL,
    theme = ttheme_minimal(
      core = list(fg_params = list(fontsize = 9)),
      colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
    )
  )
  
  pdf(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_summary_table.pdf")), 
      width = 12, height = max(4, nrow(summary_table) * 0.4))
  grid.draw(table_plot)
  dev.off()

  png(file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_summary_table.png")), 
      width = 1200, height = max(400, nrow(summary_table) * 40), res = 100)
  grid.draw(table_plot)
  dev.off()
  
  cat("  SAVED: Divergent_", ISO1_NAME, "_summary_table.pdf\\n")
  
  write.csv(divergent_results, file.path(OUTPUT_DIR, paste0("Divergent_", ISO1_NAME, "_clones_with_", ISO2_NAME, "_size.csv")), row.names = FALSE)
  cat("  SAVED: Divergent_", ISO1_NAME, "_clones_with_", ISO2_NAME, "_size.csv\\n")
  
  cat("\\nDIVERGENT", ISO1_NAME, "COPY SIZE SUMMARY:\\n")
  cat("================================\\n")
  print(summary_table)
}}

cat("\\n================================================================================\\n")
cat("SECTION 12C COMPLETE\\n")
cat("================================================================================\\n\\n")

# =============================================================================
# SECTION 13: LINEAGE TREES AND MSA ALIGNMENTS
# =============================================================================

cat("================================================================================\\n")
cat("SECTION 13: LINEAGE TREES AND MSA ALIGNMENTS\\n")
cat("================================================================================\\n\\n")

# Try to load additional packages for trees
tree_packages_available <- TRUE
tryCatch({{
  library(ape)
  library(ggtree)
  library(Biostrings)
  library(pwalign)
  library(msa)
}}, error = function(e) {{
  cat("Warning: Some tree/MSA packages not available. Skipping lineage analysis.\\n")
  tree_packages_available <<- FALSE
}})


# Try to load ggmsa
  if(!require(ggmsa, quietly = TRUE)) {{
    if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("ggmsa")
    library(ggmsa)
  }}


if(tree_packages_available && nrow(shared_clusters) > 0) {{




  # Function to create lineage tree
 # Function to create lineage tree - supports both IMGT and MiXCR formats
  create_lineage_tree <- function(iso1_cluster_id, iso2_cluster_id,
                                  clusters_iso1, clusters_iso2,
                                  raw_iso1, raw_iso2,
                                  divergent_junctions = NULL,
                                  output_prefix, output_dir) {{
    
    cat("  Creating lineage tree for", ISO1_NAME, "cluster", iso1_cluster_id, "...\\n")
    
    iso1_cluster <- clusters_iso1[[iso1_cluster_id]]
    iso2_cluster <- clusters_iso2[[iso2_cluster_id]]
    
    iso1_indices <- which(raw_iso1$junction_aa %in% iso1_cluster$junction_aa)
    iso2_indices <- which(raw_iso2$junction_aa %in% iso2_cluster$junction_aa)
    
    if(length(iso1_indices) == 0) {{
      cat("    Warning: No", ISO1_NAME, "sequences found\\n")
      return(NULL)
    }}
    
    # =================================================================
    # AUTO-DETECT SEQUENCE COLUMN (supports both IMGT and MiXCR)
    # =================================================================
    detect_seq_column <- function(data) {{
      # Priority 1: IMGT style sequence_alignment
      if("sequence_alignment" %in% colnames(data)) {{
        n_valid <- sum(!is.na(data$sequence_alignment) & 
                       data$sequence_alignment != "" & 
                       nchar(data$sequence_alignment) > 10, na.rm = TRUE)
        if(n_valid > 0) return("sequence_alignment")
      }}
      # Priority 2: MiXCR style sequence
      if("sequence" %in% colnames(data)) {{
        n_valid <- sum(!is.na(data$sequence) & 
                       data$sequence != "" & 
                       nchar(data$sequence) > 10, na.rm = TRUE)
        if(n_valid > 0) return("sequence")
      }}
      return(NULL)
    }}
    
    seq_col_iso1 <- detect_seq_column(raw_iso1)
    seq_col_iso2 <- detect_seq_column(raw_iso2)
    
    if(is.null(seq_col_iso1)) {{
      cat("    Warning: No valid sequence column found in", ISO1_NAME, "data\\n")
      cat("    Available columns:", paste(colnames(raw_iso1)[1:min(10, ncol(raw_iso1))], collapse = ", "), "...\\n")
      return(NULL)
    }}
    
    if(is.null(seq_col_iso2)) {{
      cat("    Warning: No valid sequence column found in", ISO2_NAME, "data\\n")
      return(NULL)
    }}
    
    cat("    Using sequence columns:", seq_col_iso1, "(", ISO1_NAME, "),", seq_col_iso2, "(", ISO2_NAME, ")\\n")
    
    # =================================================================
    # BUILD SEQUENCE DATA FRAME
    # =================================================================
    seq_data <- data.frame(
      label = character(),
      sequence = character(),
      junction_aa = character(),
      isotype = character(),
      is_divergent = logical(),
      stringsAsFactors = FALSE
    )
    
    # Check for germline (IMGT style)
    germline_seq <- NA
    if("germline_alignment" %in% colnames(raw_iso1)) {{
      germ_candidates <- raw_iso1$germline_alignment[iso1_indices]
      germ_candidates <- germ_candidates[!is.na(germ_candidates) & germ_candidates != ""]
      if(length(germ_candidates) > 0) {{
        germline_seq <- germ_candidates[1]
      }}
    }}
    
    if(!is.na(germline_seq) && nchar(germline_seq) > 10) {{
      seq_data <- rbind(seq_data, data.frame(
        label = "naive",
        sequence = germline_seq,
        junction_aa = "germline",
        isotype = "Germline",
        is_divergent = FALSE,
        stringsAsFactors = FALSE
      ))
      cat("    Added germline sequence (IMGT style)\\n")
    }} else {{
      cat("    No germline sequence available\\n")
    }}
    
    # =================================================================
    # ADD ISO1 SEQUENCES
    # =================================================================
    for(idx in iso1_indices) {{
      junction <- raw_iso1$junction_aa[idx]
      seq <- raw_iso1[[seq_col_iso1]][idx]
      
      # Robust NA/empty checks
      if(is.null(seq) || length(seq) == 0) next
      if(length(seq) > 1) seq <- seq[1]  # Handle vector case
      if(is.na(seq) || seq == "") next
      if(nchar(seq) < 10) next
      
      if(is.null(junction) || length(junction) == 0) next
      if(length(junction) > 1) junction <- junction[1]
      if(is.na(junction) || junction == "") next
      
      is_div <- FALSE
      if(!is.null(divergent_junctions) && length(divergent_junctions) > 0) {{
        is_div <- junction %in% divergent_junctions
      }}
      
      clean_id <- get_clean_id(raw_iso1$sequence_id[idx])
      label <- paste0(substr(junction, 1, 15), "-", ISO1_NAME, "-ID", clean_id)
      
      seq_data <- rbind(seq_data, data.frame(
        label = label,
        sequence = seq,
        junction_aa = junction,
        isotype = ISO1_NAME,
        is_divergent = is_div,
        stringsAsFactors = FALSE
      ))
    }}
    
    # =================================================================
    # ADD ISO2 SEQUENCES
    # =================================================================
    for(idx in iso2_indices) {{
      junction <- raw_iso2$junction_aa[idx]
      seq <- raw_iso2[[seq_col_iso2]][idx]
      
      # Robust NA/empty checks
      if(is.null(seq) || length(seq) == 0) next
      if(length(seq) > 1) seq <- seq[1]
      if(is.na(seq) || seq == "") next
      if(nchar(seq) < 10) next
      
      if(is.null(junction) || length(junction) == 0) next
      if(length(junction) > 1) junction <- junction[1]
      if(is.na(junction) || junction == "") next
      
      clean_id <- get_clean_id(raw_iso2$sequence_id[idx])
      label <- paste0(substr(junction, 1, 15), "-", ISO2_NAME, "-ID", clean_id)
      
      seq_data <- rbind(seq_data, data.frame(
        label = label,
        sequence = seq,
        junction_aa = junction,
        isotype = ISO2_NAME,
        is_divergent = FALSE,
        stringsAsFactors = FALSE
      ))
    }}
    
    # Remove duplicates
    seq_data <- seq_data[!duplicated(seq_data$label), ]
    
    cat("    Collected", nrow(seq_data), "sequences\\n")
    
    # =================================================================
    # LIMIT SEQUENCES IF TOO MANY
    # =================================================================
    if(nrow(seq_data) > 100) {{
      cat("    Limiting to 100 sequences\\n")
      germline_rows <- seq_data[seq_data$isotype == "Germline", ]
      div_rows <- seq_data[seq_data$is_divergent == TRUE, ]
      other_rows <- seq_data[seq_data$isotype != "Germline" & !seq_data$is_divergent, ]
      
      n_other <- 100 - nrow(germline_rows) - nrow(div_rows)
      if(nrow(other_rows) > n_other && n_other > 0) {{
        set.seed(123)
        other_rows <- other_rows[sample(nrow(other_rows), n_other), ]
      }}
      seq_data <- rbind(germline_rows, div_rows, other_rows)
    }}
    
    if(nrow(seq_data) < 3) {{
      cat("    Warning: Too few sequences (n=", nrow(seq_data), ") - skipping tree\\n")
      return(NULL)
    }}
    
    # =================================================================
    # BUILD PHYLOGENETIC TREE
    # =================================================================
    cat("    Building tree with", nrow(seq_data), "sequences...\\n")
    
    # Trim sequences to same length
    seq_lengths <- nchar(seq_data$sequence)
    min_len <- min(seq_lengths)
    if(min_len < 50) {{
      cat("    Warning: Sequences too short (min=", min_len, ")\\n")
      return(NULL)
    }}
    seq_data$sequence <- substr(seq_data$sequence, 1, min_len)
    
    tryCatch({{
      dna_seqs <- DNAStringSet(seq_data$sequence)
      names(dna_seqs) <- seq_data$label
      
      dist_matrix <- pwalign::stringDist(dna_seqs, method = "hamming")
      tree <- nj(as.dist(dist_matrix))
      
      if("naive" %in% tree$tip.label) {{
        tree <- root(tree, outgroup = "naive", resolve.root = TRUE)
        cat("    Tree rooted at germline\\n")
      }}
      
      tip_df <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
      tip_df$isotype <- seq_data$isotype[match(tip_df$label, seq_data$label)]
      tip_df$is_divergent <- seq_data$is_divergent[match(tip_df$label, seq_data$label)]
      tip_df$isotype[is.na(tip_df$isotype)] <- "Unknown"
      tip_df$is_divergent[is.na(tip_df$is_divergent)] <- FALSE
      
      n_div <- sum(tip_df$is_divergent, na.rm = TRUE)
      n_iso1 <- sum(tip_df$isotype == ISO1_NAME, na.rm = TRUE)
      n_iso2 <- sum(tip_df$isotype == ISO2_NAME, na.rm = TRUE)
      
      p <- ggtree(tree, layout = "rectangular") %<+% tip_df +
        geom_tiplab(aes(color = isotype), size = 2.5, align = FALSE, hjust = -0.02) +
        geom_tippoint(aes(color = isotype, 
                          shape = ifelse(is_divergent, "Divergent", "Normal")), 
                      size = 3) +
        geom_nodepoint(color = "orange", size = 2, alpha = 0.8) +
        scale_color_manual(values = c("red", "blue", "darkgreen", "gray"),
                           name = "Isotype") +
        scale_shape_manual(values = c("Divergent" = 17, "Normal" = 16),
                           name = "Status") +
        theme_tree2() +
        ggtitle(paste0("Lineage Tree - Cluster ", iso1_cluster_id, 
                       " (n=", nrow(seq_data), 
                       ", ", ISO1_NAME, "=", n_iso1, ", ", ISO2_NAME, "=", n_iso2,
                       if(n_div > 0) paste0(", ", n_div, " divergent") else "",
                       ")")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
              legend.position = "right") +
        xlim(0, max(node.depth.edgelength(tree)) * 1.5)
      
      pdf_file <- file.path(output_dir, paste0(output_prefix, "_lineage_tree.pdf"))
      png_file <- file.path(output_dir, paste0(output_prefix, "_lineage_tree.png"))
      tree_height <- max(10, nrow(seq_data) * 0.3)
      ggsave(pdf_file, plot = p, width = 14, height = tree_height)
      ggsave(png_file, plot = p, width = 14, height = tree_height, dpi = 150)
      
      nwk_file <- file.path(output_dir, paste0(output_prefix, "_tree.nwk"))
      write.tree(tree, file = nwk_file)
      
      cat("    SAVED:", basename(pdf_file), "\\n")
      
      return(list(tree = tree, plot = p, data = seq_data))
      
    }}, error = function(e) {{
      cat("    ERROR in tree building:", e$message, "\\n")
      return(NULL)
    }})
  }}
# =============================================================================
  # FUNCTION 2: Create MSA Alignment
  # =============================================================================

  create_msa_alignment <- function(iso1_cluster_id, iso2_cluster_id,
                                   clusters_iso1, clusters_iso2,
                                   raw_iso1, raw_iso2,
                                   divergent_junctions = NULL,
                                   output_prefix, output_dir) {{
    
    cat("  Creating MSA for cluster", iso1_cluster_id, "...\\n")
    
    iso1_cluster <- clusters_iso1[[iso1_cluster_id]]
    iso2_cluster <- clusters_iso2[[iso2_cluster_id]]
    
    iso1_indices <- which(raw_iso1$junction_aa %in% iso1_cluster$junction_aa)
    iso2_indices <- which(raw_iso2$junction_aa %in% iso2_cluster$junction_aa)
    
    all_seqs <- c()
    all_labels <- c()
    
    for(idx in iso1_indices) {{
      junction <- raw_iso1$junction_aa[idx]
      if(is.na(junction) || junction == "") next
      
      clean_id <- get_clean_id(raw_iso1$sequence_id[idx])
      is_div <- !is.null(divergent_junctions) && junction %in% divergent_junctions
      
      label <- paste0(junction, "-", ISO1_NAME, "-ID", clean_id, if(is_div) "_DIV" else "")
      
      all_seqs <- c(all_seqs, junction)
      all_labels <- c(all_labels, label)
    }}
    
    for(idx in iso2_indices) {{
      junction <- raw_iso2$junction_aa[idx]
      if(is.na(junction) || junction == "") next
      
      clean_id <- get_clean_id(raw_iso2$sequence_id[idx])
      label <- paste0(junction, "-", ISO2_NAME, "-ID", clean_id)
      
      all_seqs <- c(all_seqs, junction)
      all_labels <- c(all_labels, label)
    }}
    
    dup_idx <- !duplicated(all_labels)
    all_seqs <- all_seqs[dup_idx]
    all_labels <- all_labels[dup_idx]
    
    if(length(all_seqs) > 30) {{
      div_idx <- grep("_DIV", all_labels)
      other_idx <- setdiff(seq_along(all_seqs), div_idx)
      
      n_keep <- 30 - length(div_idx)
      if(length(other_idx) > n_keep && n_keep > 0) {{
        set.seed(123)
        other_idx <- sample(other_idx, n_keep)
      }}
      
      keep_idx <- c(div_idx, other_idx)
      all_seqs <- all_seqs[keep_idx]
      all_labels <- all_labels[keep_idx]
    }}
    
    if(length(all_seqs) < 2) {{
      cat("    Warning: Too few sequences\\n")
      return(NULL)
    }}
    
    cat("    Aligning", length(all_seqs), "sequences...\\n")
    
    tryCatch({{
      aa_seqs <- AAStringSet(all_seqs)
      names(aa_seqs) <- all_labels
      
      alignment <- msa::msa(aa_seqs, method = "ClustalOmega")
      
      fasta_file <- file.path(output_dir, paste0(output_prefix, "_alignment.fasta"))
      writeXStringSet(AAStringSet(alignment), file = fasta_file)
      cat("    SAVED:", basename(fasta_file), "\\n")
      
      aligned_seqs <- AAStringSet(alignment)
      
      p_msa <- ggmsa(aligned_seqs, 
                     start = 1, 
                     end = width(aligned_seqs)[1],
                     char_width = 0.5,
                     seq_name = TRUE,
                     color = "Chemistry_AA") +
        ggtitle(paste0("Junction MSA - Cluster ", iso1_cluster_id)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
      
      msa_pdf <- file.path(output_dir, paste0(output_prefix, "_MSA.pdf"))
      msa_png <- file.path(output_dir, paste0(output_prefix, "_MSA.png"))
      msa_height <- max(6, length(all_seqs) * 0.3)
      ggsave(msa_pdf, plot = p_msa, width = 14, height = msa_height)
      ggsave(msa_png, plot = p_msa, width = 14, height = msa_height, dpi = 150)
      cat("    SAVED:", basename(msa_pdf), "and .png\\n")
      
      tryCatch({{
        pretty_pdf <- file.path(output_dir, paste0(output_prefix, "_MSA_pretty.pdf"))
        msa::msaPrettyPrint(alignment, 
                            file = pretty_pdf,
                            output = "pdf",
                            showConsensus = "bottom",
                            showLogo = "top",
                            logoColors = "chemistry",
                            askForOverwrite = FALSE,
                            verbose = FALSE)
        cat("    SAVED:", basename(pretty_pdf), "\\n")
      }}, error = function(e) {{
        cat("    Note: msaPrettyPrint skipped\\n")
      }})
      
      return(alignment)
      
    }}, error = function(e) {{
      cat("    ERROR:", e$message, "\\n")
      return(NULL)
    }})
  }}
 
  # Generate trees for top shared clusters
  cat("\\n--- TOP SHARED CLUSTERS ---\\n\\n")
  
  n_shared <- min(N_TOP_SHARED_CLUSTERS, nrow(shared_clusters))
  
  for(i in 1:n_shared) {{
    iso1_id <- shared_clusters$iso1_cluster[i]
    iso2_id <- shared_clusters$iso2_cluster[i]
    output_prefix <- paste0("Shared_Rank", i, "_", ISO1_NAME, iso1_id, "_", ISO2_NAME, iso2_id)
    
    cat("\\n[", i, "/", n_shared, "]", ISO1_NAME, iso1_id, "-", ISO2_NAME, iso2_id, "\\n")
    
    create_lineage_tree(iso1_id, iso2_id, clusters_iso1, clusters_iso2,
                        raw_data_iso1, raw_data_iso2, NULL, output_prefix, OUTPUT_DIR)
    
    create_msa_alignment(iso1_id, iso2_id, clusters_iso1, clusters_iso2,
                         raw_data_iso1, raw_data_iso2, NULL, output_prefix, OUTPUT_DIR)
  }}

  # Generate trees for divergent clusters
  cat("\\n\\n--- DIVERGENT", ISO1_NAME, "CLUSTERS ---\\n\\n")
  
  if(nrow(divergent_results) > 0) {{
    
    for(mut_level in c("Low", "Moderate", "High")) {{
      level_data <- divergent_results[divergent_results$mutation_level == mut_level, ]
      
      if(nrow(level_data) > 0) {{
        cat("\\n===", toupper(mut_level), "MUTATION (n=", nrow(level_data), ") ===\\n")
        
        level_data <- level_data[order(-level_data$num_divergent), ]
        n_trees <- min(N_TOP_DIVERGENT_CLUSTERS, nrow(level_data))
        
        for(i in 1:n_trees) {{
          iso1_id <- level_data$iso1_cluster_id[i]
          iso2_id <- level_data$iso2_cluster_id[i]
          div_junctions <- strsplit(level_data$divergent_junctions[i], "; ")[[1]]
          output_prefix <- paste0("Divergent_", mut_level, "_Rank", i, "_", ISO1_NAME, iso1_id, "_", ISO2_NAME, iso2_id)
          
          cat("\\n[", i, "]", level_data$num_divergent[i], "divergent in", ISO1_NAME, iso1_id, "\\n")
          
          create_lineage_tree(iso1_id, iso2_id, clusters_iso1, clusters_iso2,
                              raw_data_iso1, raw_data_iso2, div_junctions, output_prefix, OUTPUT_DIR)
          
          create_msa_alignment(iso1_id, iso2_id, clusters_iso1, clusters_iso2,
                               raw_data_iso1, raw_data_iso2, div_junctions, output_prefix, OUTPUT_DIR)
        }}
      }}
    }}
  }}
}} else {{
  cat("Skipping lineage tree generation (packages not available or no shared clusters)\\n")
}}

cat("\\n\\n================================================================================\\n")
cat("SECTION 13 COMPLETE\\n")
cat("================================================================================\\n\\n")


# =============================================================================
# SECTION 14: CLUSTER SUMMARY STATISTICS
# =============================================================================

cat("\\n")
cat("================================================================================\\n")
cat("SECTION 14: CLUSTER SUMMARY STATISTICS\\n")
cat("================================================================================\\n\\n")

cat("Generating cluster summary statistics...\\n")

clusters_summary_iso1 <- Clusters.summary(pro_data_list = pro_sample_list_iso1, 
                                          clusters_list = cluster_list_iso1)
print(clusters_summary_iso1[[ISO1_SAMPLE_NAME]])

clusters_summary_iso2 <- Clusters.summary(pro_data_list = pro_sample_list_iso2, 
                                          clusters_list = cluster_list_iso2)
print(clusters_summary_iso2[[ISO2_SAMPLE_NAME]])


# =============================================================================
# SECTION 14B: EXPORT STATISTICS TO CSV FILES
# =============================================================================

cat("\\n")
cat("================================================================================\\n")
cat("SECTION 14B: EXPORTING STATISTICS TO CSV FILES\\n")
cat("================================================================================\\n\\n")

# Create summary statistics dataframe
summary_stats <- data.frame(
  Statistic = c(
    paste0("Total ", ISO1_NAME, " sequences"),
    paste0("Total ", ISO2_NAME, " sequences"),
    paste0(ISO1_NAME, " clusters"),
    paste0(ISO2_NAME, " clusters"),
    "Junction similarity threshold (%)",
    "Shared cluster pairs",
    paste0(ISO1_NAME, " clusters overlapping with ", ISO2_NAME),
    paste0(ISO1_NAME, " overlap percentage (%)"),
    paste0(ISO2_NAME, " clusters overlapping with ", ISO1_NAME),
    paste0(ISO2_NAME, " overlap percentage (%)"),
    "---DIVERGENT ANALYSIS---",
    "Minimum AA differences for divergence",
    paste0("Total divergent ", ISO1_NAME, " clusters"),
    paste0("Total divergent ", ISO1_NAME, " sequences"),
    paste0("Divergent ", ISO1_NAME, " - Low mutation clusters"),
    paste0("Divergent ", ISO1_NAME, " - Low mutation sequences"),
    paste0("Divergent ", ISO1_NAME, " - Moderate mutation clusters"),
    paste0("Divergent ", ISO1_NAME, " - Moderate mutation sequences"),
    paste0("Divergent ", ISO1_NAME, " - High mutation clusters"),
    paste0("Divergent ", ISO1_NAME, " - High mutation sequences"),
    "---NON-DIVERGENT ANALYSIS---",
    paste0("Total ", ISO1_NAME, " sequences in shared clusters"),
    paste0("Total NON-divergent ", ISO1_NAME, " clusters"),
    paste0("Total NON-divergent ", ISO1_NAME, " sequences"),
    paste0("NON-divergent ", ISO1_NAME, " percentage (%)"),
    paste0("NON-divergent ", ISO1_NAME, " - Low mutation clusters"),
    paste0("NON-divergent ", ISO1_NAME, " - Low mutation sequences"),
    paste0("NON-divergent ", ISO1_NAME, " - Moderate mutation clusters"),
    paste0("NON-divergent ", ISO1_NAME, " - Moderate mutation sequences"),
    paste0("NON-divergent ", ISO1_NAME, " - High mutation clusters"),
    paste0("NON-divergent ", ISO1_NAME, " - High mutation sequences")
  ),
  Value = c(
    nrow(raw_data_iso1),
    nrow(raw_data_iso2),
    length(clusters_iso1),
    length(clusters_iso2),
    JUNCTION_SIMILARITY_THRESHOLD,
    nrow(shared_clusters),
    length(unique_iso1_in_shared),
    round(length(unique_iso1_in_shared)/length(clusters_iso1)*100, 1),
    length(unique_iso2_in_shared),
    round(length(unique_iso2_in_shared)/length(clusters_iso2)*100, 1),
    "",
    DIVERGENT_MIN_DIFFERENCES,
    nrow(divergent_results),
    if(nrow(divergent_results) > 0) sum(divergent_results$num_divergent) else 0,
    if(nrow(divergent_results) > 0) sum(divergent_results$mutation_level == "Low") else 0,
    if(nrow(divergent_results) > 0) sum(divergent_results$num_divergent[divergent_results$mutation_level == "Low"]) else 0,
    if(nrow(divergent_results) > 0) sum(divergent_results$mutation_level == "Moderate") else 0,
    if(nrow(divergent_results) > 0) sum(divergent_results$num_divergent[divergent_results$mutation_level == "Moderate"]) else 0,
    if(nrow(divergent_results) > 0) sum(divergent_results$mutation_level == "High") else 0,
    if(nrow(divergent_results) > 0) sum(divergent_results$num_divergent[divergent_results$mutation_level == "High"]) else 0,
    "",
    total_iso1_in_shared,
    total_nondiv_clusters,
    total_nondivergent_seqs,
    if(total_iso1_in_shared > 0) round(total_nondivergent_seqs/total_iso1_in_shared*100, 1) else 0,
    if(nrow(nondivergent_results) > 0) nondivergent_results$nondivergent_clusters[nondivergent_results$mutation_level == "Low"] else 0,
    if(nrow(nondivergent_results) > 0) nondivergent_results$nondivergent_sequences[nondivergent_results$mutation_level == "Low"] else 0,
    if(nrow(nondivergent_results) > 0) nondivergent_results$nondivergent_clusters[nondivergent_results$mutation_level == "Moderate"] else 0,
    if(nrow(nondivergent_results) > 0) nondivergent_results$nondivergent_sequences[nondivergent_results$mutation_level == "Moderate"] else 0,
    if(nrow(nondivergent_results) > 0) nondivergent_results$nondivergent_clusters[nondivergent_results$mutation_level == "High"] else 0,
    if(nrow(nondivergent_results) > 0) nondivergent_results$nondivergent_sequences[nondivergent_results$mutation_level == "High"] else 0
  ),
  stringsAsFactors = FALSE
)

# Export summary statistics as CSV
write.csv(summary_stats, file.path(OUTPUT_DIR, paste0(ISO1_NAME, "_", ISO2_NAME, "_summary_statistics.csv")), row.names = FALSE)
cat("  SAVED:", paste0(ISO1_NAME, "_", ISO2_NAME, "_summary_statistics.csv"), "\\n")

cat("\\n  CSV files exported successfully!\\n\\n")



# =============================================================================
# SECTION 15: FINAL SUMMARY REPORT
# =============================================================================

cat("\\n")
cat("================================================================================\\n")
cat("SECTION 15: ANALYSIS COMPLETE - FINAL SUMMARY\\n")
cat("================================================================================\\n\\n")

cat("INPUT DATA:\\n")
cat("  ", ISO1_NAME, "sequences:", nrow(raw_data_iso1), "\\n")
cat("  ", ISO2_NAME, "sequences:", nrow(raw_data_iso2), "\\n\\n")

cat("CLUSTERING RESULTS:\\n")
cat("  ", ISO1_NAME, "clusters:", length(clusters_iso1), "\\n")
cat("  ", ISO2_NAME, "clusters:", length(clusters_iso2), "\\n\\n")

cat("SHARED CLUSTER ANALYSIS:\\n")
cat("  Junction similarity threshold:", JUNCTION_SIMILARITY_THRESHOLD, "%\\n")
cat("  Shared cluster pairs found:", nrow(shared_clusters), "\\n")
cat("  ", ISO1_NAME, "clusters in", ISO2_NAME, ":", length(unique_iso1_in_shared), 
    "(", round(length(unique_iso1_in_shared)/length(clusters_iso1)*100, 1), "%)\\n\\n")

if(nrow(shared_iso1_shm) > 0) {{
  cat("SHM THRESHOLDS:\\n")
  cat("  Low: <", round(low_threshold, 2), "%\\n")
  cat("  Moderate:", round(low_threshold, 2), "-", round(high_threshold, 2), "%\\n")
  cat("  High: >", round(high_threshold, 2), "%\\n\\n")
}}

cat("DIVERGENT", ISO1_NAME, "ANALYSIS:\\n")
cat("  Minimum AA differences:", DIVERGENT_MIN_DIFFERENCES, "\\n")
cat("  Total clusters with divergent", ISO1_NAME, ":", nrow(divergent_results), "\\n")

if(nrow(divergent_results) > 0) {{
  cat("    - Low mutation:", sum(divergent_results$mutation_level == "Low"), "clusters,",
      sum(divergent_results$num_divergent[divergent_results$mutation_level == "Low"]), "sequences\\n")
  cat("    - Moderate mutation:", sum(divergent_results$mutation_level == "Moderate"), "clusters,",
      sum(divergent_results$num_divergent[divergent_results$mutation_level == "Moderate"]), "sequences\\n")
  cat("    - High mutation:", sum(divergent_results$mutation_level == "High"), "clusters,",
      sum(divergent_results$num_divergent[divergent_results$mutation_level == "High"]), "sequences\\n\\n")
}}

cat("NON-DIVERGENT", ISO1_NAME, "ANALYSIS:\\n")
cat("  Total non-divergent", ISO1_NAME, "sequences:", total_nondivergent_seqs, 
    "(", round(total_nondivergent_seqs/total_iso1_in_shared*100, 1), "%)\\n")
cat("  Total non-divergent", ISO1_NAME, "clusters:", total_nondiv_clusters, "\\n")
if(nrow(nondivergent_results) > 0) {{
  cat("    - Low mutation:", nondivergent_results$nondivergent_clusters[nondivergent_results$mutation_level == "Low"], "clusters,",
      nondivergent_results$nondivergent_sequences[nondivergent_results$mutation_level == "Low"], "sequences\\n")
  cat("    - Moderate mutation:", nondivergent_results$nondivergent_clusters[nondivergent_results$mutation_level == "Moderate"], "clusters,",
      nondivergent_results$nondivergent_sequences[nondivergent_results$mutation_level == "Moderate"], "sequences\\n")
  cat("    - High mutation:", nondivergent_results$nondivergent_clusters[nondivergent_results$mutation_level == "High"], "clusters,",
      nondivergent_results$nondivergent_sequences[nondivergent_results$mutation_level == "High"], "sequences\\n\\n")
}}

cat("OUTPUT FILES CREATED:\\n")
cat("  CSV files:\\n")
cat("    -", paste0(ISO1_NAME, "_cluster_SHM.csv"), "\\n")
cat("    -", paste0(ISO2_NAME, "_cluster_SHM.csv"), "\\n")
cat("    -", paste0(ISO1_NAME, "_", ISO2_NAME, "_shared_clusters.csv"), "\\n")
cat("    -", paste0("Shared_", ISO1_NAME, "_clusters_SHM.csv"), "\\n")
cat("    -", paste0("Divergent_", ISO1_NAME, "_clones_all_levels.csv"), "\\n")
cat("  Visualizations:\\n")
cat("    - Cluster bubble plots, V/J gene usage, VJ pairing heatmaps\\n")
cat("    - Junction similarity distribution\\n")
cat("    - Multi-panel shared cluster analysis\\n")
cat("    - Divergent", ISO1_NAME, "by mutation level\\n")
cat("  Lineage trees and MSA:\\n")
cat("    - Top", N_TOP_SHARED_CLUSTERS, "shared cluster trees\\n")
cat("    - Divergent cluster trees for all mutation levels\\n\\n")

# Calculate and report runtime and memory usage
analysis_end_time <- Sys.time()
total_runtime_secs <- as.numeric(difftime(analysis_end_time, analysis_start_time, units = "secs"))
total_runtime_mins <- as.numeric(difftime(analysis_end_time, analysis_start_time, units = "mins"))

cat("\\n================================================================================\\n")
cat("[MEMORY] PERFORMANCE STATISTICS\\n")
cat("================================================================================\\n")
cat("[MEMORY] Analysis start time:", format(analysis_start_time, "%Y-%m-%d %H:%M:%S"), "\\n")
cat("[MEMORY] Analysis end time:", format(analysis_end_time, "%Y-%m-%d %H:%M:%S"), "\\n")
cat("[MEMORY] Total runtime:", round(total_runtime_mins, 2), "minutes (", round(total_runtime_secs, 1), "seconds)\\n")

# Memory usage statistics
final_memory <- gc()
cat("[MEMORY] Memory used (Ncells):", final_memory[1,2], "Mb\\n")
cat("[MEMORY] Memory used (Vcells):", final_memory[2,2], "Mb\\n")
cat("[MEMORY] Total memory consumed:", sum(final_memory[,2]), "Mb\\n")
cat("[MEMORY] Max memory used (Ncells):", final_memory[1,6], "Mb\\n")
cat("[MEMORY] Max memory used (Vcells):", final_memory[2,6], "Mb\\n")
cat("[MEMORY] Peak memory consumption:", sum(final_memory[,6]), "Mb\\n")

# System info
cat("[MEMORY] R version:", R.version$version.string, "\\n")
cat("[MEMORY] Platform:", R.version$platform, "\\n")
cat("[MEMORY] OS:", Sys.info()["sysname"], Sys.info()["release"], "\\n")
cat("================================================================================\\n")

# Capture and display all warnings in detail
cat("\\n================================================================================\\n")
cat("[WARNING] DETAILED WARNING REPORT\\n")
cat("================================================================================\\n")

all_warnings <- warnings()
if(length(all_warnings) > 0) {{
  cat("[WARNING] Total warnings encountered:", length(all_warnings), "\\n\\n")
  
  for(i in seq_along(all_warnings)) {{
    warning_call <- names(all_warnings)[i]
    warning_msg <- as.character(all_warnings[[i]])
    
    cat("[WARNING]", i, ":\\n")
    cat("[WARNING]   Function:", warning_call, "\\n")
    cat("[WARNING]   Message:", warning_msg, "\\n")
    cat("[WARNING]   ---\\n")
  }}
}} else {{
  cat("[WARNING] No warnings were generated during analysis.\\n")
}}

cat("================================================================================\\n")

# Clean up parallel cluster (Windows)
if(IS_WINDOWS && exists("CL") && !is.null(CL)) {{
  stopCluster(CL)
  cat("[PERFORMANCE] Parallel cluster stopped\\n")
}}

cat("\\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
cat("================================================================================\\n")

sessionInfo_file <- file.path(OUTPUT_DIR, "session_info.txt")
writeLines(capture.output(sessionInfo()), sessionInfo_file)
'''
        
        # Write to temp file
        temp_script_path = os.path.join(self.output_dir, "_temp_fastbcr_analysis.R")
        
        with open(temp_script_path, 'w') as f:
            f.write(script_content)
        
        return temp_script_path
    
    def stop(self):
        self._is_running = False
        if self.process:
            self.process.terminate()


# =============================================================================
# FASTBCR ANALYSIS TAB
# =============================================================================

class FastBCRTab(QWidget):
    """Tab for FastBCR (VJ + Junction Similarity) analysis"""
    
    def __init__(self):
        super().__init__()
        self.r_thread = None
        self.current_log_content = ""  # ADD THIS LINE
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Title and description
        title_layout = QVBoxLayout()
        title = QLabel("FastBCR Analysis (VJ + Kmer)")
        title.setFont(QFont("Arial", 16, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("color: #cccccc;")
        title_layout.addWidget(title)
        
        description = QLabel(
            "This method uses V+J genes and k-mer seeding for clustering within each isotype. "
            "Shared clusters between isotypes are identified using V+J gene matching and "
            "junction amino acid similarity. Requires AIRR-formatted TSV files (db-pass.tsv)."
        )
        description.setAlignment(Qt.AlignCenter)
        description.setStyleSheet("color: #888888; font-size: 11px;")
        title_layout.addWidget(description)
        
        layout.addLayout(title_layout)
        
        # R Status
        self.r_status_label = QLabel()
        self.check_r_status()
        layout.addWidget(self.r_status_label)
        
        # Scroll area for inputs
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Isotype Configuration
        isotype_group = QGroupBox("Isotype Configuration")
        isotype_group.setStyleSheet(self.get_group_style())
        isotype_layout = QGridLayout()
        
        isotype_layout.addWidget(QLabel("Isotype 1 Name:"), 0, 0)
        self.iso1_name = QLineEdit("IgE")
        self.iso1_name.setStyleSheet(self.get_input_style())
        isotype_layout.addWidget(self.iso1_name, 0, 1)
        
        isotype_layout.addWidget(QLabel("Isotype 2 Name:"), 0, 2)
        self.iso2_name = QLineEdit("IgG1")
        self.iso2_name.setStyleSheet(self.get_input_style())
        isotype_layout.addWidget(self.iso2_name, 0, 3)
        
        isotype_group.setLayout(isotype_layout)
        scroll_layout.addWidget(isotype_group)
        
        # File inputs
        file_group = QGroupBox("Input Files (AIRR-formatted TSV)")
        file_group.setStyleSheet(self.get_group_style())
        file_layout = QVBoxLayout()
        
        # Iso1 file
        iso1_layout = QHBoxLayout()
        self.iso1_label = QLabel("Isotype 1 TSV File:")
        iso1_layout.addWidget(self.iso1_label)
        self.iso1_input = QLineEdit()
        self.iso1_input.setStyleSheet(self.get_input_style())
        iso1_layout.addWidget(self.iso1_input)
        iso1_browse = QPushButton("Browse")
        iso1_browse.setStyleSheet(self.get_button_style())
        iso1_browse.clicked.connect(lambda: self.browse_file(self.iso1_input))
        iso1_layout.addWidget(iso1_browse)
        file_layout.addLayout(iso1_layout)
        
        # Iso2 file
        iso2_layout = QHBoxLayout()
        self.iso2_label = QLabel("Isotype 2 TSV File:")
        iso2_layout.addWidget(self.iso2_label)
        self.iso2_input = QLineEdit()
        self.iso2_input.setStyleSheet(self.get_input_style())
        iso2_layout.addWidget(self.iso2_input)
        iso2_browse = QPushButton("Browse")
        iso2_browse.setStyleSheet(self.get_button_style())
        iso2_browse.clicked.connect(lambda: self.browse_file(self.iso2_input))
        iso2_layout.addWidget(iso2_browse)
        file_layout.addLayout(iso2_layout)
        
        # Output directory
        output_layout = QHBoxLayout()
        output_layout.addWidget(QLabel("Output Directory:"))
        self.output_input = QLineEdit()
        self.output_input.setText(str(Path.home() / "Downloads" / "FastBCR_Analysis"))
        self.output_input.setStyleSheet(self.get_input_style())
        output_layout.addWidget(self.output_input)
        output_browse = QPushButton("Browse")
        output_browse.setStyleSheet(self.get_button_style())
        output_browse.clicked.connect(self.browse_output_dir)
        output_layout.addWidget(output_browse)
        file_layout.addLayout(output_layout)
        
        file_group.setLayout(file_layout)
        scroll_layout.addWidget(file_group)
        
        # Logging Options Section (INSERT THIS ENTIRE BLOCK)
        logging_group = QGroupBox("Logging Options")
        logging_group.setStyleSheet(self.get_group_style())
        logging_layout = QHBoxLayout()
        
        self.enable_logging_checkbox = QCheckBox("Enable Detailed Logging")
        self.enable_logging_checkbox.setChecked(True)
        self.enable_logging_checkbox.setStyleSheet("color: #cccccc;")
        
        logging_info = QLabel("(Tracks: clusters, exclusions, singletons, non-productive, runtime)")
        logging_info.setStyleSheet("color: #666666; font-size: 10px;")
        
        self.view_log_button = QPushButton("View Log")
        self.view_log_button.setEnabled(False)
        self.view_log_button.setStyleSheet(self.get_button_style())
        self.view_log_button.clicked.connect(self.show_log_viewer)
        
        logging_layout.addWidget(self.enable_logging_checkbox)
        logging_layout.addWidget(logging_info)
        logging_layout.addStretch()
        logging_layout.addWidget(self.view_log_button)
        logging_group.setLayout(logging_layout)
        scroll_layout.addWidget(logging_group)
        
        # Analysis Parameters
        param_group = QGroupBox("Analysis Parameters")
        
        # Analysis Parameters
        param_group = QGroupBox("Analysis Parameters")
        param_group.setStyleSheet(self.get_group_style())
        param_layout = QGridLayout()
        
        # Junction similarity threshold
        param_layout.addWidget(QLabel("Junction Similarity Threshold (%):"), 0, 0)
        self.junction_threshold = QSpinBox()
        self.junction_threshold.setRange(50, 100)
        self.junction_threshold.setValue(80)
        self.junction_threshold.setStyleSheet(self.get_input_style())
        param_layout.addWidget(self.junction_threshold, 0, 1)
        
        help_label1 = QLabel("(Min % similarity for shared clusters)")
        help_label1.setStyleSheet("color: #666666; font-size: 10px;")
        param_layout.addWidget(help_label1, 0, 2)
        
        # Divergent AA difference threshold
        param_layout.addWidget(QLabel("Divergent AA Differences (≥):"), 1, 0)
        self.min_aa_diff = QSpinBox()
        self.min_aa_diff.setRange(1, 10)
        self.min_aa_diff.setValue(2)
        self.min_aa_diff.setStyleSheet(self.get_input_style())
        param_layout.addWidget(self.min_aa_diff, 1, 1)
        
        help_label2 = QLabel("(Min AA differences for divergent clones)")
        help_label2.setStyleSheet("color: #666666; font-size: 10px;")
        param_layout.addWidget(help_label2, 1, 2)
        
        # Minimum cluster depth
        param_layout.addWidget(QLabel("Minimum Cluster Depth:"), 2, 0)
        self.min_depth = QSpinBox()
        self.min_depth.setRange(1, 20)
        self.min_depth.setValue(3)
        self.min_depth.setStyleSheet(self.get_input_style())
        param_layout.addWidget(self.min_depth, 2, 1)
        
        help_label3 = QLabel("(Min sequences per cluster)")
        help_label3.setStyleSheet("color: #666666; font-size: 10px;")
        param_layout.addWidget(help_label3, 2, 2)
        
        # Productive only checkbox
        # Productive only checkbox
        self.productive_only = QCheckBox("Productive sequences only")
        self.productive_only.setChecked(False)
        self.productive_only.setStyleSheet("color: #cccccc;")
        param_layout.addWidget(self.productive_only, 3, 0, 1, 2)
        
        param_group.setLayout(param_layout)
        scroll_layout.addWidget(param_group)

        # Divergent Clone Information Note
        divergent_info_frame = QFrame()
        divergent_info_frame.setStyleSheet("""
            QFrame {
                background-color: rgba(255, 152, 0, 0.15);
                border: 1px solid #FF9800;
                border-radius: 5px;
                padding: 5px;
                margin: 5px 0px;
            }
        """)
        divergent_info_layout = QHBoxLayout()
        divergent_info_layout.setContentsMargins(10, 8, 10, 8)
        
        warning_icon = QLabel("⚠️")
        warning_icon.setFont(QFont("Arial", 14))
        
        divergent_note = QLabel(
            "<b>Note:</b> Divergent clones are identified for <b style='color: #000000;'>Isotype 1</b> only. "
            "These are sequences in shared clusters that differ by ≥ the specified AA differences from all Isotype 2 sequences."
        )
        divergent_note.setWordWrap(True)
        divergent_note.setStyleSheet("color: #000000; font-size: 11px;")
        
        divergent_info_layout.addWidget(warning_icon)
        divergent_info_layout.addWidget(divergent_note, 1)
        divergent_info_frame.setLayout(divergent_info_layout)
        scroll_layout.addWidget(divergent_info_frame)
        
        scroll.setWidget(scroll_widget)
        layout.addWidget(scroll)
        
        # Connect isotype name changes to update labels
        self.iso1_name.textChanged.connect(self.update_labels)
        self.iso2_name.textChanged.connect(self.update_labels)
        
        # Progress and log
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)  # Indeterminate
        self.progress_bar.setVisible(False)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #666666;
                border-radius: 6px;
                text-align: center;
                background-color: #2a2a2a;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
            }
        """)
        layout.addWidget(self.progress_bar)
        
        self.log_output = QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setFont(QFont("Courier", 9))
        self.log_output.setMaximumHeight(250)
        self.log_output.setStyleSheet("""
            QPlainTextEdit {
                background-color: #1a1a1a;
                color: #00ff00;
                border: 1px solid #444444;
                border-radius: 5px;
            }
        """)
        layout.addWidget(self.log_output)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.start_button = QPushButton("Start FastBCR Analysis")
        self.start_button.clicked.connect(self.start_analysis)
        self.start_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                padding: 12px 24px;
                font-size: 14px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover { background-color: #45a049; }
            QPushButton:disabled { background-color: #666666; }
        """)
        button_layout.addWidget(self.start_button)
        
        # ADD THIS - Missing Stop Button
        self.stop_button = QPushButton("Stop Analysis")
        self.stop_button.clicked.connect(self.stop_analysis)
        self.stop_button.setEnabled(False)
        self.stop_button.setStyleSheet("""
            QPushButton {
                background-color: #f44336;
                color: white;
                padding: 12px 24px;
                font-size: 14px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover { background-color: #d32f2f; }
            QPushButton:disabled { background-color: #666666; }
        """)
        button_layout.addWidget(self.stop_button)
        
        self.view_plots_button = QPushButton("View Plots")
        self.view_plots_button.clicked.connect(self.view_plots)
        self.view_plots_button.setEnabled(False)
        self.view_plots_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                padding: 12px 24px;
                font-size: 14px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover { background-color: #1976D2; }
            QPushButton:disabled { background-color: #666666; }
        """)
        button_layout.addWidget(self.view_plots_button)
        
        self.clear_button = QPushButton("Clear Log")
        self.clear_button.clicked.connect(self.clear_log)
        self.clear_button.setStyleSheet(self.get_button_style())
        button_layout.addWidget(self.clear_button)
        
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def get_group_style(self):
        return """
            QGroupBox {
                color: #cccccc;
                border: 1px solid #666666;
                border-radius: 5px;
                margin-top: 10px;
                padding: 15px;
                font-weight: bold;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
        """
    
    def get_input_style(self):
        return """
            QLineEdit, QSpinBox, QComboBox {
                background-color: #2a2a2a;
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px;
                border-radius: 4px;
            }
            QLineEdit:hover, QSpinBox:hover, QComboBox:hover {
                border: 1px solid #999999;
            }
        """
    
    def get_button_style(self):
        return """
            QPushButton {
                background-color: #555555;
                color: white;
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #777777;
            }
        """
    
    def check_r_status(self):
        if Config.check_r_installation():
            if Config.check_fastbcr_installed():
                self.r_status_label.setText("✓ R and fastBCR are installed and available")
                self.r_status_label.setStyleSheet("color: green; font-weight: bold; padding: 5px;")
            else:
                self.r_status_label.setText("⚠ R is installed but fastBCR package is missing. It will be installed on first run.")
                self.r_status_label.setStyleSheet("color: orange; font-weight: bold; padding: 5px;")
        else:
            self.r_status_label.setText("✗ R not found. Please install R from https://cran.r-project.org/")
            self.r_status_label.setStyleSheet("color: red; font-weight: bold; padding: 5px;")
    
    def update_labels(self):
        """Update file input labels when isotype names change"""
        iso1 = self.iso1_name.text() or "Isotype 1"
        iso2 = self.iso2_name.text() or "Isotype 2"
        self.iso1_label.setText(f"{iso1} TSV File:")
        self.iso2_label.setText(f"{iso2} TSV File:")
    
    def browse_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select TSV File", "", "TSV Files (*.tsv);;All Files (*)"
        )
        if file_path:
            line_edit.setText(file_path)
    
    def browse_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_input.setText(dir_path)
    
    def start_analysis(self):
        # Validate inputs
        if not self.iso1_input.text() or not self.iso2_input.text():
            QMessageBox.warning(self, "Missing Input", "Please select both isotype files.")
            return
        
        if not Config.check_r_installation():
            QMessageBox.critical(self, "R Not Found", 
                "R is not installed or not in PATH. Please install R first.")
            return
        
        # Clear previous log content (ADD THIS LINE)
        self.current_log_content = ""
        
        # Create output directory
        output_dir = self.output_input.text()
        os.makedirs(output_dir, exist_ok=True)
        
        # Get parameters
        parameters = {
            'iso1_name': self.iso1_name.text() or "IgE",
            'iso2_name': self.iso2_name.text() or "IgG1",
            'junction_threshold': self.junction_threshold.value(),
            'min_aa_diff': self.min_aa_diff.value(),
            'min_depth': self.min_depth.value(),
            'productive_only': self.productive_only.isChecked()
        }
        
        # Start analysis thread
        self.r_thread = FastBCRAnalysisThread(
            iso1_file=self.iso1_input.text(),
            iso2_file=self.iso2_input.text(),
            output_dir=output_dir,
            parameters=parameters
        )
        
        self.r_thread.progress_updated.connect(self.update_log)
        self.r_thread.analysis_completed.connect(self.on_complete)
        self.r_thread.error_occurred.connect(self.on_error)
        
        self.r_thread.start()
        
        self.start_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.log_output.appendPlainText(f"\n{'='*60}\nStarting FastBCR Analysis at {datetime.now()}\n{'='*60}\n")
    
    def stop_analysis(self):
        if self.r_thread:
            self.r_thread.stop()
            self.log_output.appendPlainText("\n[Analysis stopped by user]\n")
        self.reset_ui()
    
    def update_log(self, message):
        self.log_output.appendPlainText(message)
        self.current_log_content += message + "\n"  # ADD THIS LINE
        scrollbar = self.log_output.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())
    
    def on_complete(self, message):
        self.log_output.appendPlainText(f"\n{'='*60}\n{message}\n{'='*60}\n")
        self.current_log_content += f"\n{'='*60}\n{message}\n{'='*60}\n"
        QMessageBox.information(self, "Complete", message)
        self.reset_ui()
        self.view_log_button.setEnabled(True)
        self.view_plots_button.setEnabled(True)  # Also enable View Plots
    
    def on_error(self, message):
        self.log_output.appendPlainText(f"\n[ERROR] {message}\n")
        QMessageBox.critical(self, "Error", message)
        self.reset_ui()
    
    def reset_ui(self):
        self.start_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    def clear_log(self):
        """Clear the log output and stored content."""
        self.log_output.clear()
        self.current_log_content = ""
    
    def show_log_viewer(self):
        """Open the log viewer dialog for FastBCR analysis."""
        if self.current_log_content:
            dialog = FastBCRLogViewer(self.current_log_content, self)
            dialog.exec_()
        else:
            QMessageBox.information(self, "No Log", "No log data available. Run an analysis first.")
    
    def update_log_content(self, log_line):
        """Append log line to stored content."""
        self.current_log_content += log_line + "\n"

    def view_plots(self):
        """Open the output directory or display plots in a viewer"""
        output_dir = self.output_input.text()
        
        if not output_dir or not os.path.exists(output_dir):
            QMessageBox.warning(self, "No Output", "No output directory found. Run analysis first.")
            return
        
        # Check for plot files
        plot_files = []
        for ext in ['*.pdf', '*.png', '*.jpg']:
            import glob
            plot_files.extend(glob.glob(os.path.join(output_dir, ext)))
        
        if not plot_files:
            QMessageBox.information(self, "No Plots", "No plot files found in output directory.")
            # Still open the folder
            self.open_output_folder()
            return
        
        # Show plot viewer dialog
        dialog = FastBCRPlotViewer(output_dir, self)
        dialog.exec_()
    
    def open_output_folder(self):
        """Open the output folder in file explorer"""
        output_dir = self.output_input.text()
        if output_dir and os.path.exists(output_dir):
            import subprocess
            import platform
            
            if platform.system() == "Windows":
                os.startfile(output_dir)
            elif platform.system() == "Darwin":  # macOS
                subprocess.run(["open", output_dir])
            else:  # Linux
                subprocess.run(["xdg-open", output_dir])
        else:
            QMessageBox.warning(self, "Error", "Output directory not found.")


# Add this new class before class AnalysisThread (around line 20):

class LogHandler(logging.Handler):
    """Custom log handler that emits signals for GUI updates."""
    def __init__(self, signal):
        super().__init__()
        self.signal = signal
    
    def emit(self, record):
        msg = self.format(record)
        self.signal.emit(msg)


class AnalysisLogger:
    """Logger for tracking analysis metrics and exclusions."""
    def __init__(self):
        self.logs = []
        self.start_time = None
        self.exclusion_counts = {
            'non_productive_iso1': 0,
            'non_productive_iso2': 0,
            'no_vdj_iso1': 0,
            'no_vdj_iso2': 0,
            'missing_vgene': 0,
            'missing_dgene': 0,
            'missing_jgene': 0,
            'low_copy_number': 0,
        }
        self.validation_results = []
        self.memory_snapshots = []
    
    def start_timer(self):
        self.start_time = time.time()
        self.log_info("Analysis started")
        self.log_memory("Initial memory state")
    
    def get_elapsed_time(self):
        if self.start_time:
            return time.time() - self.start_time
        return 0
    
    def log_memory(self, label=""):
        if not PSUTIL_AVAILABLE:
            self.logs.append(f"[MEMORY] {label}: psutil not available")
            return
        try:
            process = psutil.Process()
            mem_info = process.memory_info()
            mem_mb = mem_info.rss / (1024 * 1024)
            snapshot = f"[MEMORY] {label}: {mem_mb:.2f} MB"
            self.memory_snapshots.append(snapshot)
            self.logs.append(snapshot)
        except Exception as e:
            self.logs.append(f"[MEMORY] Could not get memory info: {e}")
    
    def log_info(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.logs.append(f"[{timestamp}] [INFO] {message}")
    
    def log_warning(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.logs.append(f"[{timestamp}] [WARNING] {message}")
    
    def log_error(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.logs.append(f"[{timestamp}] [ERROR] {message}")
    
    def log_exclusion(self, category, count, reason=""):
        self.exclusion_counts[category] = count
        timestamp = datetime.now().strftime("%H:%M:%S")
        msg = f"[{timestamp}] [EXCLUSION] {category}: {count} sequences excluded"
        if reason:
            msg += f" - Reason: {reason}"
        self.logs.append(msg)
    
    def log_validation(self, file_type, file_path, status, details=""):
        timestamp = datetime.now().strftime("%H:%M:%S")
        result = {
            'file_type': file_type,
            'file_path': file_path,
            'status': status,
            'details': details
        }
        self.validation_results.append(result)
        self.logs.append(f"[{timestamp}] [VALIDATION] {file_type}: {status} - {details}")
    
    def log_dataframe_info(self, df_name, df, key_columns=None):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.logs.append(f"[{timestamp}] [DATAFRAME] {df_name}:")
        self.logs.append(f"    - Shape: {df.shape[0]} rows x {df.shape[1]} columns")
        if key_columns:
            for col in key_columns:
                if col in df.columns:
                    non_null = df[col].notna().sum()
                    null_count = df[col].isna().sum()
                    self.logs.append(f"    - {col}: {non_null} valid, {null_count} missing")
    
    def get_summary(self):
        elapsed = self.get_elapsed_time()
        summary = [
            "=" * 60,
            "ANALYSIS SUMMARY",
            "=" * 60,
            f"Total runtime: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)",
            "",
            "EXCLUSION SUMMARY:",
            "-" * 40,
        ]
        for category, count in self.exclusion_counts.items():
            summary.append(f"  {category}: {count}")
        
        summary.extend([
            "",
            "VALIDATION RESULTS:",
            "-" * 40,
        ])
        for result in self.validation_results:
            summary.append(f"  {result['file_type']}: {result['status']}")
        
        summary.extend([
            "",
            "MEMORY USAGE:",
            "-" * 40,
        ])
        summary.extend([f"  {snap}" for snap in self.memory_snapshots])
        
        summary.append("=" * 60)
        return "\n".join(summary)
    
    def get_full_log(self):
        return "\n".join(self.logs) + "\n\n" + self.get_summary()



class AnalysisThread(QThread):
    progress_updated = pyqtSignal(int)
    analysis_completed = pyqtSignal(str)
    error_occurred = pyqtSignal(str)
    log_updated = pyqtSignal(str)  # New signal for log updates

    def __init__(self, file_paths, isotype1_name, isotype2_name, cdr3_threshold=2, cdr3cdr2_threshold=2, enable_logging=True, output_dir=None):
        super().__init__()
        self.file_paths = file_paths
        # Store user-defined isotype names
        self.iso1 = isotype1_name  # e.g., "IgE", "IgA", "IgM"
        self.iso2 = isotype2_name  # e.g., "IgG1", "IgG2", "IgD"
        # Store mutation thresholds for divergent clone filtering
        self.cdr3_threshold = cdr3_threshold
        self.cdr3cdr2_threshold = cdr3cdr2_threshold
        # Logging
        self.enable_logging = enable_logging
        self.logger = AnalysisLogger()
        # Output directory (ADD THIS)
        if output_dir:
            self.output_dir = output_dir
        else:
            self.output_dir = os.path.join(os.path.expanduser("~"), "Downloads", "CompIgS_Analysis")

    def run(self):
        try:
            import pandas as pd
            import os
            from Bio import SeqIO
            from collections import Counter
            
            # Start logging
            self.logger.start_timer()
            self.logger.log_info(f"Starting BCR analysis: {self.iso1} vs {self.iso2}")
            self.logger.log_info(f"CDR3 threshold: >= {self.cdr3_threshold}, CDR3CDR2 threshold: >= {self.cdr3cdr2_threshold}")
            
            # Use dynamic isotype names throughout
            iso1 = self.iso1
            iso2 = self.iso2
            
            # Get mutation thresholds
            cdr3_thresh = self.cdr3_threshold
            cdr3cdr2_thresh = self.cdr3cdr2_threshold
            
            # Simulate file processing with progress updates
            for i in range(1, 101):
                self.msleep(20)  # Simulate processing delay
                self.progress_updated.emit(i)

            # Load files into DataFrames with dynamic naming
            self.logger.log_info("Loading input files...")
            
            stats_iso1 = pd.read_csv(self.file_paths[0], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso1} StatClonotype", self.file_paths[0], "LOADED", f"{len(stats_iso1)} rows")
            
            stats_iso2 = pd.read_csv(self.file_paths[1], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso2} StatClonotype", self.file_paths[1], "LOADED", f"{len(stats_iso2)} rows")
            
            stats_iso1_mut = pd.read_csv(self.file_paths[2], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso1} Mutation Stats", self.file_paths[2], "LOADED", f"{len(stats_iso1_mut)} rows")
            
            stats_iso2_mut = pd.read_csv(self.file_paths[3], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso2} Mutation Stats", self.file_paths[3], "LOADED", f"{len(stats_iso2_mut)} rows")
            
            stats_iso1_Aminotab = pd.read_csv(self.file_paths[4], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso1} AA Sequences", self.file_paths[4], "LOADED", f"{len(stats_iso1_Aminotab)} rows")
            
            stats_iso2_Aminotab = pd.read_csv(self.file_paths[5], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso2} AA Sequences", self.file_paths[5], "LOADED", f"{len(stats_iso2_Aminotab)} rows")
            
            stats_iso1_Aminotab_change = pd.read_csv(self.file_paths[6], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso1} AA Change", self.file_paths[6], "LOADED", f"{len(stats_iso1_Aminotab_change)} rows")
            
            stats_iso2_Aminotab_change = pd.read_csv(self.file_paths[7], sep='\t', low_memory=False)
            self.logger.log_validation(f"{iso2} AA Change", self.file_paths[7], "LOADED", f"{len(stats_iso2_Aminotab_change)} rows")
            
            fasta_file_path = self.file_paths[8]
            sequences = list(SeqIO.parse(fasta_file_path, "fasta"))
            self.logger.log_validation("IMGT Mouse VH Reference FASTA", fasta_file_path, "LOADED", f"{len(sequences)} sequences")
            
            # Validate germline FASTA file
            self.logger.log_info("Validating germline reference file...")
            valid_sequences = 0
            invalid_sequences = 0
            for seq in sequences:
                if len(str(seq.seq)) > 0 and all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in str(seq.seq).upper()):
                    valid_sequences += 1
                else:
                    invalid_sequences += 1
                    self.logger.log_warning(f"Invalid sequence in germline file: {seq.id}")
            self.logger.log_validation("Germline Validation", fasta_file_path, "COMPLETE", 
                                       f"{valid_sequences} valid, {invalid_sequences} invalid sequences")
            
            self.logger.log_memory("After loading all files")

            # Analysis logic starts here
            self.logger.log_info(f"Filtering productive sequences for {iso1}...")
            
            # Count exclusions before filtering
            total_iso1 = len(stats_iso1)
            non_productive_iso1 = len(stats_iso1[stats_iso1['functionality'] != "productive"])
            low_copy_iso1 = len(stats_iso1[(stats_iso1['functionality'] == "productive") & (stats_iso1['onecopy'] < 2)])
            missing_vgene_iso1 = len(stats_iso1[(stats_iso1['vgene'] == '') | (stats_iso1['vgene'].isna())])
            missing_dgene_iso1 = len(stats_iso1[(stats_iso1['dgene'] == '') | (stats_iso1['dgene'].isna())])
            missing_jgene_iso1 = len(stats_iso1[(stats_iso1['jgene'] == '') | (stats_iso1['jgene'].isna())])
            
            self.logger.log_info(f"{iso1} Initial sequences: {total_iso1}")
            self.logger.log_exclusion('non_productive_iso1', non_productive_iso1, 
                                      f"Sequences with functionality != 'productive'")
            self.logger.log_info(f"{iso1} Low copy number (< 2): {low_copy_iso1} excluded")
            self.logger.log_info(f"{iso1} Missing V-gene: {missing_vgene_iso1} excluded")
            self.logger.log_info(f"{iso1} Missing D-gene: {missing_dgene_iso1} excluded")
            self.logger.log_info(f"{iso1} Missing J-gene: {missing_jgene_iso1} excluded")
            
            # Filter productive sequences for Isotype 1
            stats_iso1_productive = stats_iso1[(stats_iso1['functionality'] == "productive") &
                                               (stats_iso1['onecopy'] >= 2) &
                                               ~(stats_iso1['vgene'] == '') & ~(stats_iso1['vgene'].isna()) &
                                               ~(stats_iso1['dgene'] == '') & ~(stats_iso1['dgene'].isna()) &
                                               ~(stats_iso1['jgene'] == '') & ~(stats_iso1['jgene'].isna())]
            
            self.logger.log_info(f"{iso1} Productive sequences after filtering: {len(stats_iso1_productive)}")

            # Create a copy of the productive Isotype 1 DataFrame
            stats_iso1_productive_VDJ = stats_iso1_productive.copy()

           # Check if the input DataFrame is empty before proceeding
            if not stats_iso1_productive.empty:
                # Concatenate vgene, dgene, and jgene to create the VDJ column
                stats_iso1_productive_VDJ = stats_iso1_productive.copy()
                stats_iso1_productive_VDJ[f'{iso1}_VDJ'] = stats_iso1_productive['vgene'].astype(str) + \
                                                          stats_iso1_productive['dgene'].astype(str) + \
                                                          stats_iso1_productive['jgene'].astype(str)

                # Remove spaces in the VDJ column
                stats_iso1_productive_VDJ[f'{iso1}_VDJ'] = stats_iso1_productive_VDJ[f'{iso1}_VDJ'].str.replace(' ', '', regex=False)

                # Calculate the number of unique clones
                unique_strings_iso1_productive_VDJ = stats_iso1_productive_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
                
                # Log clones with incomplete VDJ
                incomplete_vdj = len(stats_iso1_productive_VDJ[stats_iso1_productive_VDJ[f'{iso1}_VDJ'].str.contains('nan', case=False, na=True)])
                if incomplete_vdj > 0:
                    self.logger.log_warning(f"{iso1}: {incomplete_vdj} sequences have incomplete VDJ rearrangements")
                    self.logger.log_exclusion('no_vdj_iso1', incomplete_vdj, "Incomplete VDJ rearrangement (contains 'nan')")

                # Check if unique clones were found
                if unique_strings_iso1_productive_VDJ > 0:
                    print(f'Number of {iso1} clones: {unique_strings_iso1_productive_VDJ}')
                    self.logger.log_info(f"Number of {iso1} unique VDJ clones: {unique_strings_iso1_productive_VDJ}")
                else:
                    print(f"No {iso1}_VDJ clones found.")
                    self.logger.log_warning(f"No {iso1}_VDJ clones found after filtering")
            else:
                print(f"No entries found in the input DataFrame 'stats_{iso1}_productive'. Continuing analysis with zero clones.")
                self.logger.log_error(f"No productive {iso1} sequences found - DataFrame is empty")
                unique_strings_iso1_productive_VDJ = 0

            # Count exclusions for Isotype 2
            self.logger.log_info(f"Filtering productive sequences for {iso2}...")
            
            total_iso2 = len(stats_iso2)
            non_productive_iso2 = len(stats_iso2[stats_iso2['functionality'] != "productive"])
            low_copy_iso2 = len(stats_iso2[(stats_iso2['functionality'] == "productive") & (stats_iso2['onecopy'] < 2)])
            missing_vgene_iso2 = len(stats_iso2[(stats_iso2['vgene'] == '') | (stats_iso2['vgene'].isna())])
            missing_dgene_iso2 = len(stats_iso2[(stats_iso2['dgene'] == '') | (stats_iso2['dgene'].isna())])
            missing_jgene_iso2 = len(stats_iso2[(stats_iso2['jgene'] == '') | (stats_iso2['jgene'].isna())])
            
            self.logger.log_info(f"{iso2} Initial sequences: {total_iso2}")
            self.logger.log_exclusion('non_productive_iso2', non_productive_iso2, 
                                      f"Sequences with functionality != 'productive'")
            self.logger.log_info(f"{iso2} Low copy number (< 2): {low_copy_iso2} excluded")
            self.logger.log_info(f"{iso2} Missing V-gene: {missing_vgene_iso2} excluded")
            self.logger.log_info(f"{iso2} Missing D-gene: {missing_dgene_iso2} excluded")
            self.logger.log_info(f"{iso2} Missing J-gene: {missing_jgene_iso2} excluded")
            
            # Filter productive Isotype 2 sequences with additional criteria
            stats_iso2_productive = stats_iso2[(stats_iso2['functionality'] == "productive") &
                                               (stats_iso2['onecopy'] >= 2) &
                                               ~(stats_iso2['vgene'] == '') & ~(stats_iso2['vgene'].isna()) &
                                               ~(stats_iso2['dgene'] == '') & ~(stats_iso2['dgene'].isna()) &
                                               ~(stats_iso2['jgene'] == '') & ~(stats_iso2['jgene'].isna())]
            
            self.logger.log_info(f"{iso2} Productive sequences after filtering: {len(stats_iso2_productive)}")
            

            # Check if the filtered DataFrame is empty
            if not stats_iso2_productive.empty:
                # Create a copy and generate the VDJ column
                stats_iso2_productive_VDJ = stats_iso2_productive.copy()
                stats_iso2_productive_VDJ[f'{iso2}_VDJ'] = stats_iso2_productive['vgene'].astype(str) + \
                                                          stats_iso2_productive['dgene'].astype(str) + \
                                                          stats_iso2_productive['jgene'].astype(str)

                # Remove spaces from the VDJ column
                stats_iso2_productive_VDJ[f'{iso2}_VDJ'] = stats_iso2_productive_VDJ[f'{iso2}_VDJ'].str.replace(' ', '', regex=False)

                # Calculate the number of unique clones
                unique_strings_iso2_productive_VDJ = stats_iso2_productive_VDJ[f'{iso2}_VDJ'].nunique(dropna=True)
                print(f'Number of {iso2} clones: {unique_strings_iso2_productive_VDJ}')
            else:
                # Handle the case where no productive sequences are found
                print(f"No entries found in the input DataFrame 'stats_{iso2}_productive'. Continuing analysis with zero clones.")
                unique_strings_iso2_productive_VDJ = 0
                stats_iso2_productive_VDJ = stats_iso2_productive.copy()

            # Emit signal after successful analysis
            self.analysis_completed.emit("Initial processing completed!")

            # Make a copy for comparison
            stats_iso1_productive_VDJ_compare = stats_iso1_productive_VDJ.copy()
            stats_iso2_productive_VDJ_compare = stats_iso2_productive_VDJ.copy()

            # Compare the VDJ values between the two DataFrames
            stats_iso1_productive_VDJ_compare[f'Match_in_{iso2}_VDJ'] = stats_iso1_productive_VDJ_compare[f'{iso1}_VDJ'].apply(
                lambda x: x in set(stats_iso2_productive_VDJ_compare[f'{iso2}_VDJ']))
            stats_iso2_productive_VDJ_compare[f'Match_in_{iso1}_VDJ'] = stats_iso2_productive_VDJ_compare[f'{iso2}_VDJ'].apply(
                lambda x: x in set(stats_iso1_productive_VDJ_compare[f'{iso1}_VDJ']))

            # Filter out shared Isotype 1 clones
            shared_iso1 = stats_iso1_productive_VDJ_compare[stats_iso1_productive_VDJ_compare[f'Match_in_{iso2}_VDJ']]
            unique_strings_shared_iso1 = shared_iso1[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of shared {iso1} clones: {unique_strings_shared_iso1}')
            
            if unique_strings_iso1_productive_VDJ > 0:
                percentage_shared_iso1_clones = (unique_strings_shared_iso1 / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_shared_iso1_clones = 0
            print(f'Percentage of shared {iso1} clones: {percentage_shared_iso1_clones}')

            # Filter out Unique Isotype 1 clones
            Unique_iso1 = stats_iso1_productive_VDJ_compare[~stats_iso1_productive_VDJ_compare[f'Match_in_{iso2}_VDJ']]
            unique_strings_unique_iso1 = Unique_iso1[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of unique {iso1} clones: {unique_strings_unique_iso1}')
            
            if unique_strings_iso1_productive_VDJ > 0:
                percentage_unique_iso1_clones = (unique_strings_unique_iso1 / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_unique_iso1_clones = 0
            print(f'Percentage of unique {iso1} clones: {percentage_unique_iso1_clones}')

            # Filter out shared Isotype 2 clones
            shared_iso2 = stats_iso2_productive_VDJ_compare[stats_iso2_productive_VDJ_compare[f'Match_in_{iso1}_VDJ']]
            unique_strings_shared_iso2 = shared_iso2[f'{iso2}_VDJ'].nunique(dropna=True)
            print(f'Number of shared {iso2} clones: {unique_strings_shared_iso2}')
            
            if unique_strings_iso2_productive_VDJ > 0:
                percentage_shared_iso2_clones = (unique_strings_shared_iso2 / unique_strings_iso2_productive_VDJ) * 100
            else:
                percentage_shared_iso2_clones = 0
            print(f'Percentage of shared {iso2} clones: {percentage_shared_iso2_clones}')

            # Filter out Unique Isotype 2 clones
            Unique_iso2 = stats_iso2_productive_VDJ_compare[~stats_iso2_productive_VDJ_compare[f'Match_in_{iso1}_VDJ']]
            unique_strings_unique_iso2 = Unique_iso2[f'{iso2}_VDJ'].nunique(dropna=True)
            print(f'Number of unique {iso2} clones: {unique_strings_unique_iso2}')
            
            if unique_strings_iso2_productive_VDJ > 0:
                percentage_unique_iso2_clones = (unique_strings_unique_iso2 / unique_strings_iso2_productive_VDJ) * 100
            else:
                percentage_unique_iso2_clones = 0
            print(f'Percentage of unique {iso2} clones: {percentage_unique_iso2_clones}')

            # Plot a graph of number of shared and unique clones
            types = [iso1, iso2, iso1, iso1, iso2, iso2]
            categories = [f'Total {iso1}', f'Total {iso2}', f'Shared {iso1}', f'Unique {iso1}', f'Shared {iso2}', f'Unique {iso2}']
            counts = [unique_strings_iso1_productive_VDJ, unique_strings_iso2_productive_VDJ, 
                      unique_strings_shared_iso1, unique_strings_unique_iso1,
                      unique_strings_shared_iso2, unique_strings_unique_iso2]

            # Create a DataFrame
            Graph1 = pd.DataFrame({
                'Type': types,
                'Category': categories,
                'Count': counts
            })

            # Set Seaborn style for scientific plots
            sns.set(style="whitegrid")

            # Filter data for each isotype
            df_iso1 = Graph1[Graph1['Type'] == iso1]
            df_iso2 = Graph1[Graph1['Type'] == iso2]

            # Calculate dynamic y-axis limits
            max_iso1 = df_iso1['Count'].max() if not df_iso1.empty else 0
            max_iso2 = df_iso2['Count'].max() if not df_iso2.empty else 0

            ylim_iso1 = max_iso1 + 10
            ylim_iso2 = max_iso2 + 100

            # Generate dynamic color palette
            unique_categories_iso1 = df_iso1['Category'].nunique()
            unique_categories_iso2 = df_iso2['Category'].nunique()

            palette_iso1 = sns.color_palette("Set2", unique_categories_iso1)
            palette_iso2 = sns.color_palette("Set2", unique_categories_iso2)

            # Create two subplots
            fig1, axes = plt.subplots(1, 2, figsize=(12, 6))

            # Plot for Isotype 1
            sns.barplot(
                data=df_iso1,
                x='Category',
                y='Count',
                hue='Category',
                palette=palette_iso1,
                dodge=False,
                ax=axes[0]
            )
            axes[0].set_title(f'Graph A: {iso1}', fontsize=14)
            axes[0].set_xlabel('', fontsize=12)
            axes[0].set_ylabel('Number of Clones', fontsize=12)
            axes[0].set_ylim(0, ylim_iso1)
            axes[0].legend([], [], frameon=False)
            for container in axes[0].containers:
                axes[0].bar_label(container, fmt='%d', label_type='edge', fontsize=10)

            # Plot for Isotype 2
            sns.barplot(
                data=df_iso2,
                x='Category',
                y='Count',
                hue='Category',
                palette=palette_iso2,
                dodge=False,
                ax=axes[1]
            )
            axes[1].set_title(f'Graph B: {iso2}', fontsize=14)
            axes[1].set_xlabel('', fontsize=12)
            axes[1].set_ylabel('')
            axes[1].set_ylim(0, ylim_iso2)
            axes[1].legend([], [], frameon=False)
            for container in axes[1].containers:
                axes[1].bar_label(container, fmt='%d', label_type='edge', fontsize=10)

            plt.tight_layout()

           # Create output folder with dynamic naming (USE self.output_dir)
            plots_folder = os.path.join(self.output_dir, f"Plots", f"Shared_{iso1}")
            os.makedirs(plots_folder, exist_ok=True)

            plot_path = os.path.join(plots_folder, f"Number_of_shared_and_unique_{iso1}_{iso2}_clones.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            self.logger.log_memory("After generating clone distribution plot")

            # Sort out shared Isotype 1 for VDJ rearrangement with two V gene call
            shared_iso1_two_vgene = shared_iso1.loc[
                shared_iso1[f'{iso1}_VDJ'].str.contains('or', case=False, na=False)
            ].copy()

            shared_iso1_two_vgene[f'{iso1}_VDJ'] = shared_iso1_two_vgene[f'{iso1}_VDJ'].apply(
                lambda x: x.split('or', 1)[-1].strip() if isinstance(x, str) and 'or' in x else x
            )

            shared_iso1_two_vgene_1 = shared_iso1_two_vgene.copy()
            shared_iso1_two_vgene_1[f'Match_in_Shared_{iso1}_VDJ'] = shared_iso1_two_vgene_1[f'{iso1}_VDJ'].apply(
                lambda x: x in set(shared_iso1[f'{iso1}_VDJ']))

            if shared_iso1_two_vgene_1[f'Match_in_Shared_{iso1}_VDJ'].any():
                shared_iso1_two_vgene_2 = shared_iso1_two_vgene_1[shared_iso1_two_vgene_1[f'Match_in_Shared_{iso1}_VDJ']]
                shared_iso1_two_vgene_3 = shared_iso1_two_vgene_2.drop(columns=[f'Match_in_Shared_{iso1}_VDJ'])
                shared_iso1_main = pd.concat([shared_iso1, shared_iso1_two_vgene_3], ignore_index=True)
            else:
                print(f"No matches found in shared {iso1} VDJ.")
                shared_iso1_main = shared_iso1.copy()

            # Sort out shared Isotype 2 for VDJ rearrangement with two V gene call
            shared_iso2_two_vgene = shared_iso2.loc[
                shared_iso2[f'{iso2}_VDJ'].str.contains('or', case=False, na=False)
            ].copy()

            shared_iso2_two_vgene[f'{iso2}_VDJ'] = shared_iso2_two_vgene[f'{iso2}_VDJ'].apply(
                lambda x: x.split('or', 1)[-1].strip() if isinstance(x, str) and 'or' in x else x
            )

            shared_iso2_two_vgene_1 = shared_iso2_two_vgene.copy()
            shared_iso2_two_vgene_1[f'Match_in_Shared_{iso2}_VDJ'] = shared_iso2_two_vgene_1[f'{iso2}_VDJ'].apply(
                lambda x: x in set(shared_iso2[f'{iso2}_VDJ']))

            if shared_iso2_two_vgene_1[f'Match_in_Shared_{iso2}_VDJ'].any():
                shared_iso2_two_vgene_2 = shared_iso2_two_vgene_1[shared_iso2_two_vgene_1[f'Match_in_Shared_{iso2}_VDJ']]
                shared_iso2_two_vgene_3 = shared_iso2_two_vgene_2.drop(columns=[f'Match_in_Shared_{iso2}_VDJ'])
                shared_iso2_main = pd.concat([shared_iso2, shared_iso2_two_vgene_3], ignore_index=True)
            else:
                print(f"No matches found in shared {iso2} VDJ.")
                shared_iso2_main = shared_iso2.copy()

            # Extract the cdr3aa, onecopy, VDJ columns
            shared_iso1_main_1 = shared_iso1_main[['cdr3aa', 'onecopy', f'{iso1}_VDJ']]

            # Group by VDJ and sum the onecopy values
            subtotal_shared_iso1 = shared_iso1_main_1.groupby(f'{iso1}_VDJ', as_index=False)['onecopy'].sum()
            subtotal_shared_iso1.rename(columns={'onecopy': 'subtotal_onecopy'}, inplace=True)

            shared_iso2_main_1 = shared_iso2_main[['cdr3aa', 'onecopy', f'{iso2}_VDJ']]
            subtotal_shared_iso2 = shared_iso2_main_1.groupby(f'{iso2}_VDJ', as_index=False)['onecopy'].sum()
            subtotal_shared_iso2.rename(columns={'onecopy': 'subtotal_onecopy'}, inplace=True)

            # Calculate the ratio
            if len(subtotal_shared_iso1) == len(subtotal_shared_iso2):
                ratio1 = subtotal_shared_iso1['subtotal_onecopy'] / subtotal_shared_iso2['subtotal_onecopy']
                ratio2 = subtotal_shared_iso2['subtotal_onecopy'] / subtotal_shared_iso1['subtotal_onecopy']
                subtotal_shared_iso1[f'{iso1}_{iso2}_ratio'] = ratio1
                subtotal_shared_iso2[f'{iso2}_{iso1}_ratio'] = ratio2

            # Log transform
            subtotal_shared_iso1_log = subtotal_shared_iso1.copy()
            subtotal_shared_iso2_log = subtotal_shared_iso2.copy()
            subtotal_shared_iso1_log[f'log10_{iso1}_{iso2}_onecopy'] = np.log10(subtotal_shared_iso1_log[f'{iso1}_{iso2}_ratio'])
            subtotal_shared_iso2_log[f'log10_{iso2}_{iso1}_onecopy'] = np.log10(subtotal_shared_iso2_log[f'{iso2}_{iso1}_ratio'])

            # Sort the DataFrame
            subtotal_shared_iso1_log_sorted = subtotal_shared_iso1_log.sort_values(by=f'log10_{iso1}_{iso2}_onecopy', ascending=False)

            sns.set(style="whitegrid")

            # Create a color map for positive and negative values
            colors = subtotal_shared_iso1_log_sorted[f'log10_{iso1}_{iso2}_onecopy'].apply(lambda y: 'green' if y >= 0 else 'red')

            # Plot the bar graph
            fig2 = plt.figure(figsize=(12, 6))
            bars = plt.bar(subtotal_shared_iso1_log_sorted[f'{iso1}_VDJ'],
                           subtotal_shared_iso1_log_sorted[f'log10_{iso1}_{iso2}_onecopy'],
                           color=colors)

            plt.xlabel('VH Gene', fontsize=12)
            plt.ylabel(f'Log10({iso1}/{iso2})', fontsize=12)
            plt.title(f'Ratio of {iso1}/{iso2} Copy Number', fontsize=14)
            plt.xticks(rotation=90, ha='center', fontsize=4)

            legend_labels = [f'{iso1}-biased clones', f'{iso2}-biased clones']
            handles = [plt.Rectangle((0, 0), 1, 1, color='green'), plt.Rectangle((0, 0), 1, 1, color='red')]
            plt.legend(handles, legend_labels, title="", fontsize=10)

            plt.tight_layout()

            plots_folder_base = os.path.join(self.output_dir, f"Plots")
            os.makedirs(plots_folder_base, exist_ok=True)

            plot_path = os.path.join(plots_folder_base, f"Shared_{iso1}", f"{iso1}_{iso2}_copy_number_ratio.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')

            # Plot clonal distribution of shared clones
            subtotal_shared_iso1_log_copy_sort = shared_iso1_main_1.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            subtotal_shared_iso2_log_copy_sort = shared_iso2_main_1.sort_values(by='onecopy', ascending=True).reset_index(drop=True)

            mean1 = subtotal_shared_iso1_log_copy_sort['onecopy'].mean()
            mean2 = subtotal_shared_iso2_log_copy_sort['onecopy'].mean()

            sns.set(style="whitegrid")

            fig3, axes = plt.subplots(2, 2, figsize=(14, 12))

            # Plot for Isotype 1 - Original
            sns.scatterplot(data=subtotal_shared_iso1_log_copy_sort, x=subtotal_shared_iso1_log_copy_sort.index,
                            y='onecopy', color='blue', s=100, ax=axes[0, 0], label=f'Shared {iso1}', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[0, 0].set_title(f'Clonal Size Distribution - Shared {iso1}', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()

            # Plot for Isotype 2 - Original
            sns.scatterplot(data=subtotal_shared_iso2_log_copy_sort, x=subtotal_shared_iso2_log_copy_sort.index,
                            y='onecopy', color='red', s=100, ax=axes[0, 1], label=f'Shared {iso2}', edgecolor='w')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[0, 1].set_title(f'Clonal Size Distribution - Shared {iso2}', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()

            # Define zoom range
            zoom_range_1 = (mean2 - 20, mean2 + 20)
            zoom_range_2 = (mean2 - 20, mean2 + 20)

            shared_iso1_zoom = subtotal_shared_iso1_log_copy_sort[
                (subtotal_shared_iso1_log_copy_sort['onecopy'] >= zoom_range_1[0]) &
                (subtotal_shared_iso1_log_copy_sort['onecopy'] <= zoom_range_1[1])]
            shared_iso2_zoom = subtotal_shared_iso2_log_copy_sort[
                (subtotal_shared_iso2_log_copy_sort['onecopy'] >= zoom_range_2[0]) &
                (subtotal_shared_iso2_log_copy_sort['onecopy'] <= zoom_range_2[1])]

            # Zoomed plots
            sns.scatterplot(data=shared_iso1_zoom, x=shared_iso1_zoom.index, y='onecopy', color='blue',
                            s=100, ax=axes[1, 0], label=f'Shared {iso1}', edgecolor='w')
            axes[1, 0].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[1, 0].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[1, 0].set_title(f'Zoomed Clonal Size Distribution - Shared {iso1}', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()

            sns.scatterplot(data=shared_iso2_zoom, x=shared_iso2_zoom.index, y='onecopy', color='red',
                            s=100, ax=axes[1, 1], label=f'Shared {iso2}', edgecolor='w')
            axes[1, 1].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[1, 1].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[1, 1].set_title(f'Zoomed Clonal Size Distribution - Shared {iso2}', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()

            plt.tight_layout()

            plot_path = os.path.join(plots_folder_base, f"Clonal_size_distribution_of_Shared_{iso1}_and_{iso2}.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')

            # Select biased clones
            iso1_biased_clones = subtotal_shared_iso1_log.copy()
            iso1_biased_clones = subtotal_shared_iso1_log[subtotal_shared_iso1_log[f'log10_{iso1}_{iso2}_onecopy'] >= 0.301]

            iso1_iso2_biased_clones_VDJ_compare = shared_iso1_main.copy()
            iso1_iso2_biased_clones_VDJ_compare[f'Match_in_{iso1}_biased_VDJ'] = iso1_iso2_biased_clones_VDJ_compare[f'{iso1}_VDJ'].apply(
                lambda x: x in set(iso1_biased_clones[f'{iso1}_VDJ']))

            iso1_biased_clones_VDJ_compare = iso1_iso2_biased_clones_VDJ_compare[iso1_iso2_biased_clones_VDJ_compare[f'Match_in_{iso1}_biased_VDJ']]

            iso2_biased_clones = subtotal_shared_iso2_log.copy()
            iso2_biased_clones = subtotal_shared_iso2_log[subtotal_shared_iso2_log[f'log10_{iso2}_{iso1}_onecopy'] >= 0.301]

            iso1_iso2_biased_clones_VDJ_compare[f'Match_in_{iso2}_biased_VDJ'] = iso1_iso2_biased_clones_VDJ_compare[f'{iso1}_VDJ'].apply(
                lambda x: x in set(iso2_biased_clones[f'{iso2}_VDJ']))

            iso2_biased_clones_VDJ_compare = iso1_iso2_biased_clones_VDJ_compare[iso1_iso2_biased_clones_VDJ_compare[f'Match_in_{iso2}_biased_VDJ']]

            # Cleanup Isotype 1 mutation table
            stats_iso1_mut_1 = stats_iso1_mut[stats_iso1_mut['V-DOMAIN Functionality'] == "productive"].copy()

            stats_iso1_mut_1.loc[:, 'V-REGION Nb of nucleotides-ex'] = stats_iso1_mut_1['V-REGION Nb of nucleotides'].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True)
            stats_iso1_mut_1.loc[:, 'V-REGION Nb of mutations-ex'] = stats_iso1_mut_1['V-REGION Nb of mutations'].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True)
            stats_iso1_mut_1.loc[:, 'CDR3-IMGT Nb of nonsilent mutations-ex'] = stats_iso1_mut_1['CDR3-IMGT Nb of nonsilent mutations'].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True)

            stats_iso1_mut_2 = stats_iso1_mut_1.copy()

            stats_iso1_mut_2['V-REGION Nb of nucleotides-ex'] = pd.to_numeric(stats_iso1_mut_2['V-REGION Nb of nucleotides-ex'], errors='coerce')
            stats_iso1_mut_2['CDR2-IMGT Nb of nonsilent mutations'] = pd.to_numeric(stats_iso1_mut_2['CDR2-IMGT Nb of nonsilent mutations'], errors='coerce')
            stats_iso1_mut_2['CDR3-IMGT Nb of nonsilent mutations-ex'] = pd.to_numeric(stats_iso1_mut_2['CDR3-IMGT Nb of nonsilent mutations-ex'], errors='coerce')

            stats_iso1_mut_2['percentage of CDR2-IMGT Nb of nonsilent mutations'] = (
                stats_iso1_mut_2['CDR2-IMGT Nb of nonsilent mutations'] /
                stats_iso1_mut_2['V-REGION Nb of nucleotides-ex']
            ) * 100

            stats_iso1_mut_2['percentage of CDR3-IMGT Nb of nonsilent mutations'] = (
                stats_iso1_mut_2['CDR3-IMGT Nb of nonsilent mutations-ex'] /
                stats_iso1_mut_2['V-REGION Nb of nucleotides-ex']
            ) * 100

            # Compare seqid of the mutation table to the seqid of biased clones
            stats_iso1_mut_3 = stats_iso1_mut_2.copy()
            stats_iso1_mut_3[f'Match_in_{iso1}_biased_seqid'] = stats_iso1_mut_3['Sequence ID'].apply(
                lambda x: x in set(iso1_biased_clones_VDJ_compare['seqid']))
            stats_iso1_mut_3[f'Match_in_{iso2}_biased_seqid'] = stats_iso1_mut_3['Sequence ID'].apply(
                lambda x: x in set(iso2_biased_clones_VDJ_compare['seqid']))
            stats_iso1_mut_3[f'Match_in_unique_{iso1}_seqid'] = stats_iso1_mut_3['Sequence ID'].apply(
                lambda x: x in set(Unique_iso1['seqid']))
            stats_iso1_mut_3[f'Match_in_Total_{iso1}_seqid'] = stats_iso1_mut_3['Sequence ID'].apply(
                lambda x: x in set(stats_iso1_productive_VDJ['seqid']))

            # Convert to numeric
            stats_iso1_mut_3['V-REGION Nb of mutations-ex'] = pd.to_numeric(stats_iso1_mut_3['V-REGION Nb of mutations-ex'], errors='coerce')

            # Convert boolean columns to string
            stats_iso1_mut_3[f'Match_in_{iso1}_biased_seqid'] = stats_iso1_mut_3[f'Match_in_{iso1}_biased_seqid'].astype(str).str.strip()
            stats_iso1_mut_3[f'Match_in_{iso2}_biased_seqid'] = stats_iso1_mut_3[f'Match_in_{iso2}_biased_seqid'].astype(str).str.strip()
            stats_iso1_mut_3[f'Match_in_unique_{iso1}_seqid'] = stats_iso1_mut_3[f'Match_in_unique_{iso1}_seqid'].astype(str).str.strip()
            stats_iso1_mut_3[f'Match_in_Total_{iso1}_seqid'] = stats_iso1_mut_3[f'Match_in_Total_{iso1}_seqid'].astype(str).str.strip()

            # Filter mutated and unmutated subsets
            stats_tab_mutated_iso1_bias = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_{iso1}_biased_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_mutated_iso2_bias = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_{iso2}_biased_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_mutated_Unique_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_unique_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_mutated_Total_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_Total_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] > 0)]

            stats_tab_unmutated_iso1_bias = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_{iso1}_biased_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            stats_tab_unmutated_iso2_bias = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_{iso2}_biased_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            stats_tab_unmutated_Unique_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_unique_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            stats_tab_unmutated_Total_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_Total_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] == 0)]

            # Calculate means for CDR2 and CDR3 mutations
            mean_CDR2_stats_tab_mutated_iso1_bias = stats_tab_mutated_iso1_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            if pd.isna(mean_CDR2_stats_tab_mutated_iso1_bias):
                mean_CDR2_stats_tab_mutated_iso1_bias = 0
            print(f'Mutated {iso1}_biased average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_iso1_bias}')

            mean_CDR3_stats_tab_mutated_iso1_bias = stats_tab_mutated_iso1_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            if pd.isna(mean_CDR3_stats_tab_mutated_iso1_bias):
                mean_CDR3_stats_tab_mutated_iso1_bias = 0
            print(f'Mutated {iso1}_biased average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_iso1_bias}')

            mean_CDR2_stats_tab_mutated_iso2_bias = stats_tab_mutated_iso2_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            if pd.isna(mean_CDR2_stats_tab_mutated_iso2_bias):
                mean_CDR2_stats_tab_mutated_iso2_bias = 0
            print(f'Mutated {iso2}_biased average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_iso2_bias}')

            mean_CDR3_stats_tab_mutated_iso2_bias = stats_tab_mutated_iso2_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            if pd.isna(mean_CDR3_stats_tab_mutated_iso2_bias):
                mean_CDR3_stats_tab_mutated_iso2_bias = 0
            print(f'Mutated {iso2}_biased average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_iso2_bias}')

            mean_CDR2_stats_tab_mutated_Unique_iso1 = stats_tab_mutated_Unique_iso1['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            if pd.isna(mean_CDR2_stats_tab_mutated_Unique_iso1):
                mean_CDR2_stats_tab_mutated_Unique_iso1 = 0
            print(f'Mutated unique_{iso1} average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_Unique_iso1}')

            mean_CDR3_stats_tab_mutated_Unique_iso1 = stats_tab_mutated_Unique_iso1['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            if pd.isna(mean_CDR3_stats_tab_mutated_Unique_iso1):
                mean_CDR3_stats_tab_mutated_Unique_iso1 = 0
            print(f'Mutated unique_{iso1} average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_Unique_iso1}')

            # Handle unmutated means with empty checks
            if stats_tab_unmutated_iso1_bias.empty or stats_tab_unmutated_iso1_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR2_stats_tab_unmutated_iso1_bias = 0
            else:
                mean_CDR2_stats_tab_unmutated_iso1_bias = stats_tab_unmutated_iso1_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated {iso1}_biased average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_iso1_bias}')

            if stats_tab_unmutated_iso1_bias.empty or stats_tab_unmutated_iso1_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR3_stats_tab_unmutated_iso1_bias = 0
            else:
                mean_CDR3_stats_tab_unmutated_iso1_bias = stats_tab_unmutated_iso1_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated {iso1}_biased average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_iso1_bias}')

            if stats_tab_unmutated_iso2_bias.empty or stats_tab_unmutated_iso2_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR2_stats_tab_unmutated_iso2_bias = 0
            else:
                mean_CDR2_stats_tab_unmutated_iso2_bias = stats_tab_unmutated_iso2_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated {iso2}_biased average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_iso2_bias}')

            if stats_tab_unmutated_iso2_bias.empty or stats_tab_unmutated_iso2_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR3_stats_tab_unmutated_iso2_bias = 0
            else:
                mean_CDR3_stats_tab_unmutated_iso2_bias = stats_tab_unmutated_iso2_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated {iso2}_biased average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_iso2_bias}')

            if stats_tab_unmutated_Unique_iso1.empty or stats_tab_unmutated_Unique_iso1['percentage of CDR2-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR2_stats_tab_unmutated_Unique_iso1 = 0
            else:
                mean_CDR2_stats_tab_unmutated_Unique_iso1 = stats_tab_unmutated_Unique_iso1['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated unique_{iso1} average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_Unique_iso1}')

            if stats_tab_unmutated_Unique_iso1.empty or stats_tab_unmutated_Unique_iso1['percentage of CDR3-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR3_stats_tab_unmutated_Unique_iso1 = 0
            else:
                mean_CDR3_stats_tab_unmutated_Unique_iso1 = stats_tab_unmutated_Unique_iso1['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated unique_{iso1} average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_Unique_iso1}')

            # Barplot of CDR2 and CDR3 mutations
            data = {
                'Subset': [
                    f'Mutated {iso1}-bias', f'Mutated {iso1}-bias',
                    f'Mutated {iso2}-bias', f'Mutated {iso2}-bias',
                    f'Mutated Unique {iso1}', f'Mutated Unique {iso1}',
                    f'Unmutated {iso1}-bias', f'Unmutated {iso1}-bias',
                    f'Unmutated {iso2}-bias', f'Unmutated {iso2}-bias',
                    f'Unmutated Unique {iso1}', f'Unmutated Unique {iso1}'
                ],
                'Region': ['CDR2', 'CDR3'] * 6,
                'Percentage': [
                    mean_CDR2_stats_tab_mutated_iso1_bias, mean_CDR3_stats_tab_mutated_iso1_bias,
                    mean_CDR2_stats_tab_mutated_iso2_bias, mean_CDR3_stats_tab_mutated_iso2_bias,
                    mean_CDR2_stats_tab_mutated_Unique_iso1, mean_CDR3_stats_tab_mutated_Unique_iso1,
                    mean_CDR2_stats_tab_unmutated_iso1_bias, mean_CDR3_stats_tab_unmutated_iso1_bias,
                    mean_CDR2_stats_tab_unmutated_iso2_bias, mean_CDR3_stats_tab_unmutated_iso2_bias,
                    mean_CDR2_stats_tab_unmutated_Unique_iso1, mean_CDR3_stats_tab_unmutated_Unique_iso1
                ]
            }

            plot2mut = pd.DataFrame(data)

            fig4 = plt.figure(figsize=(10, 6))
            sns.barplot(data=plot2mut, x='Region', y='Percentage', hue='Subset', palette='deep')

            plt.title('Non-Silent Mutations in the CDR Regions', fontsize=14)
            plt.xlabel('VH Region', fontsize=12)
            plt.ylabel('Mean percentage non silent mutation(%)', fontsize=12)
            plt.legend(title='Subset', fontsize=10)
            plt.ylim(0, 0.8)
            plt.tight_layout()

            plot_path = os.path.join(plots_folder_base, f"Mean_percentage_non-silent_mutation_in_CDR2_and_CDR3_{iso1}_{iso2}.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            self.logger.log_memory("After mutation analysis")


            # Create a copy of the dataframe for seqid comparison
            stats_iso1_productive_VDJ_seqid_compare = stats_iso1_productive_VDJ.copy()

            # Compare the seqid from Total Isotype 1 to the seqid of mutated subsets
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_{iso1}_biased_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_mutated_iso1_bias['Sequence ID']))
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_{iso2}_biased_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_mutated_iso2_bias['Sequence ID']))
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_unique_{iso1}_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_mutated_Unique_iso1['Sequence ID']))
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_Total_{iso1}_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_mutated_Total_iso1['Sequence ID']))

            # Compare the seqid to unmutated subsets
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_{iso1}_biased_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_unmutated_iso1_bias['Sequence ID']))
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_{iso2}_biased_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_unmutated_iso2_bias['Sequence ID']))
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_unique_{iso1}_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_unmutated_Unique_iso1['Sequence ID']))
            stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_Total_{iso1}_seqid'] = stats_iso1_productive_VDJ_seqid_compare['seqid'].apply(
                lambda x: x in set(stats_tab_unmutated_Total_iso1['Sequence ID']))

            # Convert boolean columns to string
            for col in [f'Match_in_mutated_{iso1}_biased_seqid', f'Match_in_mutated_{iso2}_biased_seqid',
                        f'Match_in_mutated_unique_{iso1}_seqid', f'Match_in_mutated_Total_{iso1}_seqid',
                        f'Match_in_unmutated_{iso1}_biased_seqid', f'Match_in_unmutated_{iso2}_biased_seqid',
                        f'Match_in_unmutated_unique_{iso1}_seqid', f'Match_in_unmutated_Total_{iso1}_seqid']:
                stats_iso1_productive_VDJ_seqid_compare[col] = stats_iso1_productive_VDJ_seqid_compare[col].astype(str).str.strip()

            # Filter mutated subsets with VDJ
            stats_tab_mutated_iso1_bias_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_{iso1}_biased_seqid'] == 'True']
            stats_tab_mutated_iso2_bias_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_{iso2}_biased_seqid'] == 'True']
            stats_tab_mutated_unique_iso1_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_unique_{iso1}_seqid'] == 'True']
            stats_tab_mutated_Total_iso1_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_mutated_Total_{iso1}_seqid'] == 'True']

            # Filter unmutated subsets with VDJ
            stats_tab_unmutated_iso1_bias_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_{iso1}_biased_seqid'] == 'True']
            stats_tab_unmutated_iso2_bias_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_{iso2}_biased_seqid'] == 'True']
            stats_tab_unmutated_unique_iso1_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_unique_{iso1}_seqid'] == 'True']
            stats_tab_unmutated_Total_iso1_VDJ = stats_iso1_productive_VDJ_seqid_compare[
                stats_iso1_productive_VDJ_seqid_compare[f'Match_in_unmutated_Total_{iso1}_seqid'] == 'True']

            # Number and percentage of mutated clones
            unique_strings_mutated_iso1_bias_VDJ = stats_tab_mutated_iso1_bias_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of mutated {iso1}-biased clones: {unique_strings_mutated_iso1_bias_VDJ}')

            if unique_strings_iso1_productive_VDJ > 0:
                percentage_mutated_iso1_bias_clones = (unique_strings_mutated_iso1_bias_VDJ / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_mutated_iso1_bias_clones = 0
            print(f'Percentage of mutated {iso1}-biased clones: {percentage_mutated_iso1_bias_clones}')

            mean_copynumber_mutated_iso1_bias_VDJ = stats_tab_mutated_iso1_bias_VDJ['onecopy'].mean()
            mean_copynumber_mutated_iso1_bias_VDJ = 0 if pd.isna(mean_copynumber_mutated_iso1_bias_VDJ) else mean_copynumber_mutated_iso1_bias_VDJ
            print(f'Mean copy number of mutated {iso1}-biased clones: {mean_copynumber_mutated_iso1_bias_VDJ}')

            unique_strings_mutated_iso2_bias_VDJ = stats_tab_mutated_iso2_bias_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of mutated {iso2}-biased clones: {unique_strings_mutated_iso2_bias_VDJ}')

            if unique_strings_iso1_productive_VDJ > 0:
                percentage_mutated_iso2_bias_clones = (unique_strings_mutated_iso2_bias_VDJ / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_mutated_iso2_bias_clones = 0
            print(f'Percentage of mutated {iso2}-biased clones: {percentage_mutated_iso2_bias_clones}')

            mean_copynumber_mutated_iso2_bias_VDJ = stats_tab_mutated_iso2_bias_VDJ['onecopy'].mean()
            mean_copynumber_mutated_iso2_bias_VDJ = 0 if pd.isna(mean_copynumber_mutated_iso2_bias_VDJ) else mean_copynumber_mutated_iso2_bias_VDJ
            print(f'Mean copy number of mutated {iso2}-biased clones: {mean_copynumber_mutated_iso2_bias_VDJ}')

            unique_strings_mutated_unique_iso1_VDJ = stats_tab_mutated_unique_iso1_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of mutated unique {iso1} clones: {unique_strings_mutated_unique_iso1_VDJ}')

            if unique_strings_iso1_productive_VDJ > 0:
                percentage_mutated_unique_iso1_clones = (unique_strings_mutated_unique_iso1_VDJ / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_mutated_unique_iso1_clones = 0
            print(f'Percentage of mutated unique {iso1} clones: {percentage_mutated_unique_iso1_clones}')

            mean_copynumber_mutated_unique_iso1_VDJ = stats_tab_mutated_unique_iso1_VDJ['onecopy'].mean()
            mean_copynumber_mutated_unique_iso1_VDJ = 0 if pd.isna(mean_copynumber_mutated_unique_iso1_VDJ) else mean_copynumber_mutated_unique_iso1_VDJ
            print(f'Mean copy number of mutated unique {iso1} clones: {mean_copynumber_mutated_unique_iso1_VDJ}')

            unique_strings_unmutated_Total_iso1_VDJ = stats_tab_unmutated_Total_iso1_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            unique_strings_mutated_Total_iso1_VDJ = stats_tab_mutated_Total_iso1_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of mutated Total {iso1} clones: {unique_strings_mutated_Total_iso1_VDJ}')

            try:
                percentage_mutated_Total_iso1_clones = (unique_strings_mutated_Total_iso1_VDJ / 
                    (unique_strings_mutated_Total_iso1_VDJ + unique_strings_unmutated_Total_iso1_VDJ)) * 100
                print(f'Percentage of mutated Total {iso1} clones: {percentage_mutated_Total_iso1_clones}')
            except ZeroDivisionError:
                percentage_mutated_Total_iso1_clones = 0
                print(f'Percentage of mutated Total {iso1} clones: 0% (No clones present)')

            mean_copynumber_mutated_Total_iso1_VDJ = stats_tab_mutated_Total_iso1_VDJ['onecopy'].mean()
            mean_copynumber_mutated_Total_iso1_VDJ = 0 if pd.isna(mean_copynumber_mutated_Total_iso1_VDJ) else mean_copynumber_mutated_Total_iso1_VDJ
            print(f'Mean copy number of mutated Total {iso1} clones: {mean_copynumber_mutated_Total_iso1_VDJ}')

            # Number and percentage of unmutated clones
            unique_strings_unmutated_iso1_bias_VDJ = stats_tab_unmutated_iso1_bias_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated {iso1}-biased clones: {unique_strings_unmutated_iso1_bias_VDJ}')

            if unique_strings_iso1_productive_VDJ > 0:
                percentage_unmutated_iso1_bias_clones = (unique_strings_unmutated_iso1_bias_VDJ / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_unmutated_iso1_bias_clones = 0
            print(f'Percentage of unmutated {iso1}-biased clones: {percentage_unmutated_iso1_bias_clones}')

            mean_copynumber_unmutated_iso1_bias_VDJ = stats_tab_unmutated_iso1_bias_VDJ['onecopy'].mean()
            mean_copynumber_unmutated_iso1_bias_VDJ = 0 if pd.isna(mean_copynumber_unmutated_iso1_bias_VDJ) else mean_copynumber_unmutated_iso1_bias_VDJ
            print(f'Mean copy number of unmutated {iso1}-biased clones: {mean_copynumber_unmutated_iso1_bias_VDJ}')

            unique_strings_unmutated_iso2_bias_VDJ = stats_tab_unmutated_iso2_bias_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated {iso2}-biased clones: {unique_strings_unmutated_iso2_bias_VDJ}')

            if unique_strings_iso1_productive_VDJ > 0:
                percentage_unmutated_iso2_bias_clones = (unique_strings_unmutated_iso2_bias_VDJ / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_unmutated_iso2_bias_clones = 0
            print(f'Percentage of unmutated {iso2}-biased clones: {percentage_unmutated_iso2_bias_clones}')

            mean_copynumber_unmutated_iso2_bias_VDJ = stats_tab_unmutated_iso2_bias_VDJ['onecopy'].mean()
            mean_copynumber_unmutated_iso2_bias_VDJ = 0 if pd.isna(mean_copynumber_unmutated_iso2_bias_VDJ) else mean_copynumber_unmutated_iso2_bias_VDJ
            print(f'Mean copy number of unmutated {iso2}-biased clones: {mean_copynumber_unmutated_iso2_bias_VDJ}')

            unique_strings_unmutated_unique_iso1_VDJ = stats_tab_unmutated_unique_iso1_VDJ[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated unique {iso1} clones: {unique_strings_unmutated_unique_iso1_VDJ}')

            if unique_strings_iso1_productive_VDJ > 0:
                percentage_unmutated_unique_iso1_clones = (unique_strings_unmutated_unique_iso1_VDJ / unique_strings_iso1_productive_VDJ) * 100
            else:
                percentage_unmutated_unique_iso1_clones = 0
            print(f'Percentage of unmutated unique {iso1} clones: {percentage_unmutated_unique_iso1_clones}')

            mean_copynumber_unmutated_unique_iso1_VDJ = stats_tab_unmutated_unique_iso1_VDJ['onecopy'].mean()
            mean_copynumber_unmutated_unique_iso1_VDJ = 0 if pd.isna(mean_copynumber_unmutated_unique_iso1_VDJ) else mean_copynumber_unmutated_unique_iso1_VDJ
            print(f'Mean copy number of unmutated unique {iso1} clones: {mean_copynumber_unmutated_unique_iso1_VDJ}')

            print(f'Number of unmutated Total {iso1} clones: {unique_strings_unmutated_Total_iso1_VDJ}')

            total_unique_strings = unique_strings_mutated_Total_iso1_VDJ + unique_strings_unmutated_Total_iso1_VDJ
            if total_unique_strings > 0:
                percentage_unmutated_Total_iso1_clones = (unique_strings_unmutated_Total_iso1_VDJ / total_unique_strings) * 100
            else:
                percentage_unmutated_Total_iso1_clones = 0
            print(f'Percentage of unmutated Total {iso1} clones: {percentage_unmutated_Total_iso1_clones}')

            mean_copynumber_unmutated_Total_iso1_VDJ = stats_tab_unmutated_Total_iso1_VDJ['onecopy'].mean()
            mean_copynumber_unmutated_Total_iso1_VDJ = 0 if pd.isna(mean_copynumber_unmutated_Total_iso1_VDJ) else mean_copynumber_unmutated_Total_iso1_VDJ
            print(f'Mean copy number of unmutated Total {iso1} clones: {mean_copynumber_unmutated_Total_iso1_VDJ}')

            # Cleanup Isotype 2 mutation table
            # Check for column name variations
            if 'V.DOMAIN.Functionality' in stats_iso2_mut.columns:
                column_name = 'V.DOMAIN.Functionality'
            elif 'V-DOMAIN Functionality' in stats_iso2_mut.columns:
                column_name = 'V-DOMAIN Functionality'
            else:
                raise KeyError("Neither 'V.DOMAIN.Functionality' nor 'V-DOMAIN Functionality' column is present.")

            stats_iso2_mut_1 = stats_iso2_mut[stats_iso2_mut[column_name] == "productive"].copy()

            # Handle column name variations for nucleotides
            if 'V.REGION.Nb.of.nucleotides' in stats_iso2_mut_1.columns:
                col_nucleotides = 'V.REGION.Nb.of.nucleotides'
            elif 'V-REGION Nb of nucleotides' in stats_iso2_mut_1.columns:
                col_nucleotides = 'V-REGION Nb of nucleotides'
            else:
                raise KeyError("Neither 'V.REGION.Nb.of.nucleotides' nor 'V-REGION Nb of nucleotides' column is present.")

            stats_iso2_mut_1['V-REGION Nb of nucleotides-ex'] = stats_iso2_mut_1[col_nucleotides].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True)

            # Handle mutations column
            if 'V.REGION.Nb.of.mutations' in stats_iso2_mut_1.columns:
                col_mutations = 'V.REGION.Nb.of.mutations'
            elif 'V-REGION Nb of mutations' in stats_iso2_mut_1.columns:
                col_mutations = 'V-REGION Nb of mutations'
            else:
                raise KeyError("Neither 'V.REGION.Nb.of.mutations' nor 'V-REGION Nb of mutations' column is present.")

            stats_iso2_mut_1['V-REGION Nb of mutations-ex'] = stats_iso2_mut_1[col_mutations].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True)

            # Handle CDR3 column
            if 'CDR3.IMGT.Nb.of.nonsilent.mutations' in stats_iso2_mut_1.columns:
                col_cdr3 = 'CDR3.IMGT.Nb.of.nonsilent.mutations'
            elif 'CDR3-IMGT Nb of nonsilent mutations' in stats_iso2_mut_1.columns:
                col_cdr3 = 'CDR3-IMGT Nb of nonsilent mutations'
            else:
                raise KeyError("Neither 'CDR3.IMGT.Nb.of.nonsilent.mutations' nor 'CDR3-IMGT Nb of nonsilent mutations' column is present.")

            stats_iso2_mut_1['CDR3-IMGT Nb of nonsilent mutations-ex'] = stats_iso2_mut_1[col_cdr3].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True)

            stats_iso2_mut_3 = stats_iso2_mut_1.copy()

            # Handle Sequence ID column
            if 'Sequence.ID' in stats_iso2_mut_3.columns:
                col_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_iso2_mut_3.columns:
                col_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present.")

            stats_iso2_mut_3[f'Match_in_Total_{iso2}_seqid'] = stats_iso2_mut_3[col_seqid].isin(stats_iso2_productive_VDJ['seqid'])

            stats_iso2_mut_3['V-REGION Nb of mutations-ex'] = pd.to_numeric(stats_iso2_mut_3['V-REGION Nb of mutations-ex'], errors='coerce')

            # Filter mutated and unmutated Total Isotype 2
            stats_tab_mutated_Total_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_Total_{iso2}_seqid']) & (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] > 0)]

            stats_tab_unmutated_Total_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_Total_{iso2}_seqid']) & (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] == 0)]

            # Create seqid sets
            mutated_seqids = set(stats_tab_mutated_Total_iso2[col_seqid])
            unmutated_seqids = set(stats_tab_unmutated_Total_iso2[col_seqid])

            stats_iso2_productive_VDJ_seqid_compare = stats_iso2_productive_VDJ.copy()
            stats_iso2_productive_VDJ_seqid_compare[f'Match_in_mutated_Total_{iso2}_seqid'] = stats_iso2_productive_VDJ_seqid_compare['seqid'].isin(mutated_seqids)
            stats_iso2_productive_VDJ_seqid_compare[f'Match_in_unmutated_Total_{iso2}_seqid'] = stats_iso2_productive_VDJ_seqid_compare['seqid'].isin(unmutated_seqids)

            stats_iso2_productive_VDJ_seqid_compare[f'Match_in_mutated_Total_{iso2}_seqid'] = stats_iso2_productive_VDJ_seqid_compare[f'Match_in_mutated_Total_{iso2}_seqid'].astype(str).str.strip()
            stats_iso2_productive_VDJ_seqid_compare[f'Match_in_unmutated_Total_{iso2}_seqid'] = stats_iso2_productive_VDJ_seqid_compare[f'Match_in_unmutated_Total_{iso2}_seqid'].astype(str).str.strip()

            stats_tab_mutated_Total_iso2_VDJ = stats_iso2_productive_VDJ_seqid_compare[
                stats_iso2_productive_VDJ_seqid_compare[f'Match_in_mutated_Total_{iso2}_seqid'] == 'True']
            stats_tab_unmutated_Total_iso2_VDJ = stats_iso2_productive_VDJ_seqid_compare[
                stats_iso2_productive_VDJ_seqid_compare[f'Match_in_unmutated_Total_{iso2}_seqid'] == 'True']

            # Calculate Isotype 2 statistics
            unique_strings_mutated_Total_iso2_VDJ = stats_tab_mutated_Total_iso2_VDJ[f'{iso2}_VDJ'].nunique(dropna=True)
            print(f'Number of mutated Total {iso2} clones: {unique_strings_mutated_Total_iso2_VDJ}')

            unique_strings_unmutated_Total_iso2_VDJ = stats_tab_unmutated_Total_iso2_VDJ[f'{iso2}_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated Total {iso2} clones: {unique_strings_unmutated_Total_iso2_VDJ}')

            total_iso2 = unique_strings_mutated_Total_iso2_VDJ + unique_strings_unmutated_Total_iso2_VDJ
            if total_iso2 > 0:
                percentage_mutated_Total_iso2_clones = (unique_strings_mutated_Total_iso2_VDJ / total_iso2) * 100
                percentage_unmutated_Total_iso2_clones = (unique_strings_unmutated_Total_iso2_VDJ / total_iso2) * 100
            else:
                percentage_mutated_Total_iso2_clones = 0
                percentage_unmutated_Total_iso2_clones = 0
            print(f'Percentage of mutated Total {iso2} clones: {percentage_mutated_Total_iso2_clones:.2f}%')
            print(f'Percentage of unmutated Total {iso2} clones: {percentage_unmutated_Total_iso2_clones:.2f}%')

            mean_copynumber_mutated_Total_iso2_VDJ = stats_tab_mutated_Total_iso2_VDJ['onecopy'].mean()
            mean_copynumber_mutated_Total_iso2_VDJ = 0 if pd.isna(mean_copynumber_mutated_Total_iso2_VDJ) else mean_copynumber_mutated_Total_iso2_VDJ
            print(f'Mean copy number of mutated Total {iso2} clones: {mean_copynumber_mutated_Total_iso2_VDJ}')

            mean_copynumber_unmutated_Total_iso2_VDJ = stats_tab_unmutated_Total_iso2_VDJ['onecopy'].mean()
            mean_copynumber_unmutated_Total_iso2_VDJ = 0 if pd.isna(mean_copynumber_unmutated_Total_iso2_VDJ) else mean_copynumber_unmutated_Total_iso2_VDJ
            print(f'Mean copy number of unmutated Total {iso2} clones: {mean_copynumber_unmutated_Total_iso2_VDJ}')

            # Plot pie charts
            percentage_mutated_iso1 = percentage_mutated_Total_iso1_clones
            percentage_unmutated_iso1 = percentage_unmutated_Total_iso1_clones

            labels1 = [f'Mutated {iso1}', f'Unmutated {iso1}']
            sizes1 = [percentage_mutated_iso1, percentage_unmutated_iso1]
            colors1 = ['red', 'blue']
            explode1 = (0.1, 0)

            percentage_mutated_Total_iso2 = percentage_mutated_Total_iso2_clones
            percentage_unmutated_Total_iso2 = percentage_unmutated_Total_iso2_clones

            labels2 = [f'Mutated {iso2}', f'Unmutated {iso2}']
            sizes2 = [percentage_mutated_Total_iso2, percentage_unmutated_Total_iso2]
            colors2 = ['red', 'blue']
            explode2 = (0.1, 0)

            fig5, axis = plt.subplots(1, 2, figsize=(12, 6))

            axis[0].pie(sizes1, explode=explode1, labels=labels1, colors=colors1, autopct='%1.1f%%', startangle=90)
            axis[0].set_title(f'Percentage of Mutated and Unmutated {iso1}')

            axis[1].pie(sizes2, explode=explode2, labels=labels2, colors=colors2, autopct='%1.1f%%', startangle=90)
            axis[1].set_title(f'Percentage of Mutated and Unmutated {iso2}')

            plt.tight_layout()

            plot_path = os.path.join(plots_folder_base, f"Percentage_of_mutated_unmutated_{iso1}_and_{iso2}.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')

            # Handle column name variations for stats_iso2_mut_3
            if 'Sequence.ID' in stats_iso2_mut_3.columns:
                col_seqid_iso2 = 'Sequence.ID'
            elif 'Sequence ID' in stats_iso2_mut_3.columns:
                col_seqid_iso2 = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present.")
                
            stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid'] = stats_iso1_mut_3['Sequence ID'].isin(set(shared_iso1['seqid']))

            stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid'] = stats_iso2_mut_3[col_seqid_iso2].isin(set(shared_iso2['seqid']))

            stats_tab_Shared_iso1 = stats_iso1_mut_3[stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid']]
            stats_tab_Shared_iso2 = stats_iso2_mut_3[stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid']]

            stats_tab_mutated_Shared_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid']) & (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_unmutated_Shared_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid']) & (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] == 0)]

            stats_tab_mutated_Shared_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid']) & (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_unmutated_Shared_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid']) & (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] == 0)]

            mutated_seqids_Shared_iso1 = set(stats_tab_mutated_Shared_iso1['Sequence ID'])
            unmutated_seqids_Shared_iso1 = set(stats_tab_unmutated_Shared_iso1['Sequence ID'])

            mutated_seqids_Shared_iso2 = set(stats_tab_mutated_Shared_iso2[col_seqid_iso2])
            unmutated_seqids_Shared_iso2 = set(stats_tab_unmutated_Shared_iso2[col_seqid_iso2])

            stats_Shared_iso1_VDJ_seqid_compare = stats_iso1_productive_VDJ.copy()
            stats_Shared_iso2_VDJ_seqid_compare = stats_iso2_productive_VDJ.copy()

            stats_Shared_iso1_VDJ_seqid_compare[f'Match_in_mutated_shared_{iso1}_seqid'] = stats_Shared_iso1_VDJ_seqid_compare['seqid'].isin(mutated_seqids_Shared_iso1)
            stats_Shared_iso1_VDJ_seqid_compare[f'Match_in_unmutated_shared_{iso1}_seqid'] = stats_Shared_iso1_VDJ_seqid_compare['seqid'].isin(unmutated_seqids_Shared_iso1)

            stats_Shared_iso2_VDJ_seqid_compare[f'Match_in_mutated_shared_{iso2}_seqid'] = stats_Shared_iso2_VDJ_seqid_compare['seqid'].isin(mutated_seqids_Shared_iso2)
            stats_Shared_iso2_VDJ_seqid_compare[f'Match_in_unmutated_shared_{iso2}_seqid'] = stats_Shared_iso2_VDJ_seqid_compare['seqid'].isin(unmutated_seqids_Shared_iso2)

            for col in [f'Match_in_mutated_shared_{iso1}_seqid', f'Match_in_unmutated_shared_{iso1}_seqid']:
                stats_Shared_iso1_VDJ_seqid_compare[col] = stats_Shared_iso1_VDJ_seqid_compare[col].astype(str).str.strip()
            for col in [f'Match_in_mutated_shared_{iso2}_seqid', f'Match_in_unmutated_shared_{iso2}_seqid']:
                stats_Shared_iso2_VDJ_seqid_compare[col] = stats_Shared_iso2_VDJ_seqid_compare[col].astype(str).str.strip()

            stats_tab_mutated_shared_iso1_VDJ = stats_Shared_iso1_VDJ_seqid_compare[
                stats_Shared_iso1_VDJ_seqid_compare[f'Match_in_mutated_shared_{iso1}_seqid'] == 'True']
            stats_tab_unmutated_shared_iso1_VDJ = stats_Shared_iso1_VDJ_seqid_compare[
                stats_Shared_iso1_VDJ_seqid_compare[f'Match_in_unmutated_shared_{iso1}_seqid'] == 'True']

            stats_tab_mutated_shared_iso2_VDJ = stats_Shared_iso2_VDJ_seqid_compare[
                stats_Shared_iso2_VDJ_seqid_compare[f'Match_in_mutated_shared_{iso2}_seqid'] == 'True']
            stats_tab_unmutated_shared_iso2_VDJ = stats_Shared_iso2_VDJ_seqid_compare[
                stats_Shared_iso2_VDJ_seqid_compare[f'Match_in_unmutated_shared_{iso2}_seqid'] == 'True']

            # Plot copy number distribution of mutated and unmutated shared
            mutated_shared_iso1_copy_sort = stats_tab_mutated_shared_iso1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            mutated_shared_iso2_copy_sort = stats_tab_mutated_shared_iso2_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)

            mean1 = mutated_shared_iso1_copy_sort['onecopy'].mean() if not mutated_shared_iso1_copy_sort.empty else 0
            mean2 = mutated_shared_iso2_copy_sort['onecopy'].mean() if not mutated_shared_iso2_copy_sort.empty else 0

            sns.set(style="whitegrid")

            fig8, axes = plt.subplots(2, 2, figsize=(14, 12))

            sns.scatterplot(data=mutated_shared_iso1_copy_sort, x=mutated_shared_iso1_copy_sort.index, y='onecopy',
                            color='blue', s=100, ax=axes[0, 0], label=f'Mutated shared {iso1}', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[0, 0].set_title(f'Clonal Size Distribution - Mutated shared {iso1}', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()

            sns.scatterplot(data=mutated_shared_iso2_copy_sort, x=mutated_shared_iso2_copy_sort.index, y='onecopy',
                            color='red', s=100, ax=axes[0, 1], label=f'Mutated shared {iso2}', edgecolor='w')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[0, 1].set_title(f'Clonal Size Distribution - Mutated shared {iso2}', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()

            zoom_range_1 = (mean2 - 20, mean2 + 20)
            zoom_range_2 = (mean2 - 20, mean2 + 20)

            Mutated_iso1_zoom = mutated_shared_iso1_copy_sort[
                (mutated_shared_iso1_copy_sort['onecopy'] >= zoom_range_1[0]) &
                (mutated_shared_iso1_copy_sort['onecopy'] <= zoom_range_1[1])]
            Mutated_iso2_zoom = mutated_shared_iso2_copy_sort[
                (mutated_shared_iso2_copy_sort['onecopy'] >= zoom_range_2[0]) &
                (mutated_shared_iso2_copy_sort['onecopy'] <= zoom_range_2[1])]

            sns.scatterplot(data=Mutated_iso1_zoom, x=Mutated_iso1_zoom.index, y='onecopy', color='blue',
                            s=100, ax=axes[1, 0], label=f'Mutated shared {iso1}', edgecolor='w')
            axes[1, 0].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[1, 0].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[1, 0].set_title(f'Zoomed Clonal Size Distribution - Mutated shared {iso1}', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()

            sns.scatterplot(data=Mutated_iso2_zoom, x=Mutated_iso2_zoom.index, y='onecopy', color='red',
                            s=100, ax=axes[1, 1], label=f'Mutated shared {iso2}', edgecolor='w')
            axes[1, 1].axhline(mean1, color='orange', linestyle='--', label=f'{iso1} Mean: {mean1:.2f}')
            axes[1, 1].axhline(mean2, color='green', linestyle='--', label=f'{iso2} Mean: {mean2:.2f}')
            axes[1, 1].set_title(f'Zoomed Clonal Size Distribution - Mutated shared {iso2}', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()

            plt.tight_layout()

            plot_path = os.path.join(plots_folder_base, f"Clonal_size_distribution_of_mutated_shared_{iso1}_and_{iso2}.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')

            # Plot mutation distribution
            mutation_counts_iso1 = stats_tab_Shared_iso1['V-REGION Nb of mutations-ex'].tolist()
            mutation_counts_iso2 = stats_tab_Shared_iso2['V-REGION Nb of mutations-ex'].tolist()

            mutation_distribution_iso1 = Counter(mutation_counts_iso1)
            mutation_distribution_iso2 = Counter(mutation_counts_iso2)

            sorted_mutations_iso1 = sorted(mutation_distribution_iso1.items())
            sorted_mutations_iso2 = sorted(mutation_distribution_iso2.items())

            x_values_iso1 = [mutation for mutation, count in sorted_mutations_iso1]
            y_values_iso1 = [count for mutation, count in sorted_mutations_iso1]

            x_values_iso2 = [mutation for mutation, count in sorted_mutations_iso2]
            y_values_iso2 = [count for mutation, count in sorted_mutations_iso2]

            fig9 = plt.figure(figsize=(10, 6))

            plt.plot(x_values_iso1, y_values_iso1, 'bo-', markersize=8, label=f'Shared {iso1}', linewidth=1)
            plt.plot(x_values_iso2, y_values_iso2, 'rs-', markersize=8, label=f'Shared {iso2}', linewidth=1)

            plt.axvline(x=5, color='green', linestyle='--', linewidth=1)
            plt.axvline(x=11, color='orange', linestyle='--', linewidth=1)
            plt.axvline(x=12, color='red', linestyle='--', linewidth=1)

            max_y = max(y_values_iso1 + y_values_iso2) if (y_values_iso1 or y_values_iso2) else 1
            plt.text(2.5, max_y, 'Low Mutation (0-5)', color='green', fontsize=10.5, ha='right')
            plt.text(8.5, max_y, 'Moderate Mutation (6-11)', color='orange', fontsize=10.5, ha='center')
            plt.text(13, max_y, 'High Mutation (>=12)', color='red', fontsize=10.5, ha='left')

            plt.xlabel('VH region number of mutations', fontsize=12)
            plt.ylabel('Count', fontsize=12)
            plt.title('Mutation Distributions', fontsize=14)
            plt.legend()
            plt.grid(axis='both', linestyle='--', alpha=0.2)
            plt.tight_layout()

            plot_path = os.path.join(plots_folder_base, f"Mutation_distribution_of_shared_{iso1}_and_{iso2}.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')

            # Filter out low, moderate, and highly mutated shared clones
            stats_tab_mutated_low_mut_Shared_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid']) &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] >= 0) &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] <= 5)]

            stats_tab_mutated_low_mut_Shared_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid']) &
                (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] >= 0) &
                (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] <= 5)]

            stats_tab_mutated_mod_mut_Shared_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid']) &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] >= 6) &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] <= 11)]

            stats_tab_mutated_mod_mut_Shared_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid']) &
                (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] >= 6) &
                (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] <= 11)]

            stats_tab_mutated_high_mut_Shared_iso1 = stats_iso1_mut_3[
                (stats_iso1_mut_3[f'Match_in_shared_{iso1}_seqid']) &
                (stats_iso1_mut_3['V-REGION Nb of mutations-ex'] >= 12)]

            stats_tab_mutated_high_mut_Shared_iso2 = stats_iso2_mut_3[
                (stats_iso2_mut_3[f'Match_in_shared_{iso2}_seqid']) &
                (stats_iso2_mut_3['V-REGION Nb of mutations-ex'] >= 12)]

            # Create seqid sets for mutation levels
            low_mutated_seqids_Shared_iso1 = set(stats_tab_mutated_low_mut_Shared_iso1['Sequence ID'])
            moderate_mutated_seqids_Shared_iso1 = set(stats_tab_mutated_mod_mut_Shared_iso1['Sequence ID'])
            High_mutated_seqids_Shared_iso1 = set(stats_tab_mutated_high_mut_Shared_iso1['Sequence ID'])

            low_mutated_seqids_Shared_iso2 = set(stats_tab_mutated_low_mut_Shared_iso2[col_seqid_iso2])
            moderate_mutated_seqids_Shared_iso2 = set(stats_tab_mutated_mod_mut_Shared_iso2[col_seqid_iso2])
            High_mutated_seqids_Shared_iso2 = set(stats_tab_mutated_high_mut_Shared_iso2[col_seqid_iso2])

            stats_mut_Shared_iso1_VDJ_seqid_compare = stats_iso1_productive_VDJ.copy()
            stats_mut_Shared_iso2_VDJ_seqid_compare = stats_iso2_productive_VDJ.copy()

            stats_mut_Shared_iso1_VDJ_seqid_compare[f'Match_in_low_mutated_shared_{iso1}_seqid'] = stats_mut_Shared_iso1_VDJ_seqid_compare['seqid'].isin(low_mutated_seqids_Shared_iso1)
            stats_mut_Shared_iso1_VDJ_seqid_compare[f'Match_in_moderate_mutated_shared_{iso1}_seqid'] = stats_mut_Shared_iso1_VDJ_seqid_compare['seqid'].isin(moderate_mutated_seqids_Shared_iso1)
            stats_mut_Shared_iso1_VDJ_seqid_compare[f'Match_in_high_mutated_shared_{iso1}_seqid'] = stats_mut_Shared_iso1_VDJ_seqid_compare['seqid'].isin(High_mutated_seqids_Shared_iso1)

            stats_mut_Shared_iso2_VDJ_seqid_compare[f'Match_in_low_mutated_shared_{iso2}_seqid'] = stats_mut_Shared_iso2_VDJ_seqid_compare['seqid'].isin(low_mutated_seqids_Shared_iso2)
            stats_mut_Shared_iso2_VDJ_seqid_compare[f'Match_in_moderate_mutated_shared_{iso2}_seqid'] = stats_mut_Shared_iso2_VDJ_seqid_compare['seqid'].isin(moderate_mutated_seqids_Shared_iso2)
            stats_mut_Shared_iso2_VDJ_seqid_compare[f'Match_in_high_mutated_shared_{iso2}_seqid'] = stats_mut_Shared_iso2_VDJ_seqid_compare['seqid'].isin(High_mutated_seqids_Shared_iso2)

            for col in [f'Match_in_low_mutated_shared_{iso1}_seqid', f'Match_in_moderate_mutated_shared_{iso1}_seqid', f'Match_in_high_mutated_shared_{iso1}_seqid']:
                stats_mut_Shared_iso1_VDJ_seqid_compare[col] = stats_mut_Shared_iso1_VDJ_seqid_compare[col].astype(str).str.strip()
            for col in [f'Match_in_low_mutated_shared_{iso2}_seqid', f'Match_in_moderate_mutated_shared_{iso2}_seqid', f'Match_in_high_mutated_shared_{iso2}_seqid']:
                stats_mut_Shared_iso2_VDJ_seqid_compare[col] = stats_mut_Shared_iso2_VDJ_seqid_compare[col].astype(str).str.strip()

            stats_tab_low_mutated_shared_iso1_VDJ = stats_mut_Shared_iso1_VDJ_seqid_compare[
                stats_mut_Shared_iso1_VDJ_seqid_compare[f'Match_in_low_mutated_shared_{iso1}_seqid'] == 'True']
            stats_tab_moderate_mutated_shared_iso1_VDJ = stats_mut_Shared_iso1_VDJ_seqid_compare[
                stats_mut_Shared_iso1_VDJ_seqid_compare[f'Match_in_moderate_mutated_shared_{iso1}_seqid'] == 'True']
            stats_tab_high_mutated_shared_iso1_VDJ = stats_mut_Shared_iso1_VDJ_seqid_compare[
                stats_mut_Shared_iso1_VDJ_seqid_compare[f'Match_in_high_mutated_shared_{iso1}_seqid'] == 'True']

            stats_tab_low_mutated_shared_iso2_VDJ = stats_mut_Shared_iso2_VDJ_seqid_compare[
                stats_mut_Shared_iso2_VDJ_seqid_compare[f'Match_in_low_mutated_shared_{iso2}_seqid'] == 'True']
            stats_tab_moderate_mutated_shared_iso2_VDJ = stats_mut_Shared_iso2_VDJ_seqid_compare[
                stats_mut_Shared_iso2_VDJ_seqid_compare[f'Match_in_moderate_mutated_shared_{iso2}_seqid'] == 'True']
            stats_tab_high_mutated_shared_iso2_VDJ = stats_mut_Shared_iso2_VDJ_seqid_compare[
                stats_mut_Shared_iso2_VDJ_seqid_compare[f'Match_in_high_mutated_shared_{iso2}_seqid'] == 'True']

            # Plot low, moderate, and high mutated clonal size distribution
            stats_tab_low_mutated_shared_iso1_VDJ_sort = stats_tab_low_mutated_shared_iso1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            stats_tab_low_mutated_shared_iso2_VDJ_sort = stats_tab_low_mutated_shared_iso2_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)

            stats_tab_moderate_mutated_shared_iso1_VDJ_sort = stats_tab_moderate_mutated_shared_iso1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            stats_tab_moderate_mutated_shared_iso2_VDJ_sort = stats_tab_moderate_mutated_shared_iso2_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)

            stats_tab_high_mutated_shared_iso1_VDJ_sort = stats_tab_high_mutated_shared_iso1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            stats_tab_high_mutated_shared_iso2_VDJ_sort = stats_tab_high_mutated_shared_iso2_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)

            mean1 = stats_tab_low_mutated_shared_iso1_VDJ_sort['onecopy'].mean() if not stats_tab_low_mutated_shared_iso1_VDJ_sort.empty else 0
            mean2 = stats_tab_low_mutated_shared_iso2_VDJ_sort['onecopy'].mean() if not stats_tab_low_mutated_shared_iso2_VDJ_sort.empty else 0
            mean3 = stats_tab_moderate_mutated_shared_iso1_VDJ_sort['onecopy'].mean() if not stats_tab_moderate_mutated_shared_iso1_VDJ_sort.empty else 0
            mean4 = stats_tab_moderate_mutated_shared_iso2_VDJ_sort['onecopy'].mean() if not stats_tab_moderate_mutated_shared_iso2_VDJ_sort.empty else 0
            mean5 = stats_tab_high_mutated_shared_iso1_VDJ_sort['onecopy'].mean() if not stats_tab_high_mutated_shared_iso1_VDJ_sort.empty else 0
            mean6 = stats_tab_high_mutated_shared_iso2_VDJ_sort['onecopy'].mean() if not stats_tab_high_mutated_shared_iso2_VDJ_sort.empty else 0

            sns.set(style="whitegrid")

            fig10, axes = plt.subplots(3, 2, figsize=(14, 12))

            sns.scatterplot(data=stats_tab_low_mutated_shared_iso1_VDJ_sort, x=stats_tab_low_mutated_shared_iso1_VDJ_sort.index,
                            y='onecopy', color='blue', s=100, ax=axes[0, 0], label=f'low mutated shared {iso1}', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'mean low mutated {iso1}: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'mean low mutated {iso2}: {mean2:.2f}')
            axes[0, 0].set_title(f'low mutated shared {iso1}', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()

            sns.scatterplot(data=stats_tab_low_mutated_shared_iso2_VDJ_sort, x=stats_tab_low_mutated_shared_iso2_VDJ_sort.index,
                            y='onecopy', color='red', s=100, ax=axes[0, 1], label=f'low mutated shared {iso2}', edgecolor='w')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'mean low mutated {iso1}: {mean1:.2f}')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'mean low mutated {iso2}: {mean2:.2f}')
            axes[0, 1].set_title(f'low mutated shared {iso2}', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()

            sns.scatterplot(data=stats_tab_moderate_mutated_shared_iso1_VDJ_sort, x=stats_tab_moderate_mutated_shared_iso1_VDJ_sort.index,
                            y='onecopy', color='blue', s=100, ax=axes[1, 0], label=f'moderate mutated shared {iso1}', edgecolor='w')
            axes[1, 0].axhline(mean3, color='orange', linestyle='--', label=f'mean moderate mutated {iso1}: {mean3:.2f}')
            axes[1, 0].axhline(mean4, color='green', linestyle='--', label=f'mean moderate mutated {iso2}: {mean4:.2f}')
            axes[1, 0].set_title(f'moderate mutated shared {iso1}', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()

            sns.scatterplot(data=stats_tab_moderate_mutated_shared_iso2_VDJ_sort, x=stats_tab_moderate_mutated_shared_iso2_VDJ_sort.index,
                            y='onecopy', color='red', s=100, ax=axes[1, 1], label=f'moderate mutated shared {iso2}', edgecolor='w')
            axes[1, 1].axhline(mean3, color='orange', linestyle='--', label=f'mean moderate mutated {iso1}: {mean3:.2f}')
            axes[1, 1].axhline(mean4, color='green', linestyle='--', label=f'mean moderate mutated {iso2}: {mean4:.2f}')
            axes[1, 1].set_title(f'moderate mutated shared {iso2}', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()

            sns.scatterplot(data=stats_tab_high_mutated_shared_iso1_VDJ_sort, x=stats_tab_high_mutated_shared_iso1_VDJ_sort.index,
                            y='onecopy', color='blue', s=100, ax=axes[2, 0], label=f'high mutated shared {iso1}', edgecolor='w')
            axes[2, 0].axhline(mean5, color='orange', linestyle='--', label=f'mean high mutated {iso1}: {mean5:.2f}')
            axes[2, 0].axhline(mean6, color='green', linestyle='--', label=f'mean high mutated {iso2}: {mean6:.2f}')
            axes[2, 0].set_title(f'high mutated shared {iso1}', fontsize=14)
            axes[2, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[2, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[2, 0].legend()

            sns.scatterplot(data=stats_tab_high_mutated_shared_iso2_VDJ_sort, x=stats_tab_high_mutated_shared_iso2_VDJ_sort.index,
                            y='onecopy', color='red', s=100, ax=axes[2, 1], label=f'high mutated shared {iso2}', edgecolor='w')
            axes[2, 1].axhline(mean5, color='orange', linestyle='--', label=f'mean high mutated {iso1}: {mean5:.2f}')
            axes[2, 1].axhline(mean6, color='green', linestyle='--', label=f'mean high mutated {iso2}: {mean6:.2f}')
            axes[2, 1].set_title(f'high mutated shared {iso2}', fontsize=14)
            axes[2, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[2, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[2, 1].legend()

            plt.tight_layout()

            plot_path = os.path.join(plots_folder_base, f"Clonal_size_distribution_of_low_moderate_high_mutated_shared_{iso1}_and_{iso2}.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')

            print(f'mean copynumber low mutated {iso1}: {mean1}')
            print(f'mean copynumber low mutated {iso2}: {mean2}')
            print(f'mean copynumber moderate mutated {iso1}: {mean3}')
            print(f'mean copynumber moderate mutated {iso2}: {mean4}')
            print(f'mean copynumber high mutated {iso1}: {mean5}')
            print(f'mean copynumber high mutated {iso2}: {mean6}')

            # Filter top 10 most expanded clones
            def filter_top_ten(df, column):
                if df is not None and column in df.columns:
                    return df.sort_values(by=column, ascending=False).head(10)
                else:
                    return pd.DataFrame()

            stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten = filter_top_ten(stats_tab_low_mutated_shared_iso1_VDJ, 'onecopy')
            stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten = filter_top_ten(stats_tab_moderate_mutated_shared_iso1_VDJ, 'onecopy')
            stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten = filter_top_ten(stats_tab_high_mutated_shared_iso1_VDJ, 'onecopy')

            stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten = filter_top_ten(stats_tab_low_mutated_shared_iso2_VDJ, 'onecopy')
            stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten = filter_top_ten(stats_tab_moderate_mutated_shared_iso2_VDJ, 'onecopy')
            stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten = filter_top_ten(stats_tab_high_mutated_shared_iso2_VDJ, 'onecopy')

            # Calculate net charge
            def calculate_net_charge(df, sequence_column, pH=7.0, pKscale="Lehninger"):
                def net_charge(sequence):
                    try:
                        peptide = Peptide(sequence)
                        return peptide.charge(pH=pH, pKscale=pKscale)
                    except:
                        return 0
                return df[sequence_column].apply(net_charge)

            sequence_column = 'cdr3aa'

            if not stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten.empty:
                stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                    stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten, sequence_column)
            if not stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten.empty:
                stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                    stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten, sequence_column)
            if not stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten.empty:
                stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                    stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten, sequence_column)
            if not stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten.empty:
                stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                    stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten, sequence_column)
            if not stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten.empty:
                stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                    stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten, sequence_column)
            if not stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten.empty:
                stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                    stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten, sequence_column)

            # Calculate mean net charge
            mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten = 0
            mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten = 0
            mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten = 0
            mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten = 0
            mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten = 0
            mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten = 0

            if not stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten.empty and 'Net charge CDR3aa' in stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten.columns:
                mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten = stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                if pd.isna(mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten):
                    mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten = 0
            print(f'Mean CDR3 netcharge of low mutated {iso1}: {mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten}')

            if not stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten.empty and 'Net charge CDR3aa' in stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten.columns:
                mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten = stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                if pd.isna(mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten):
                    mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten = 0
            print(f'Mean CDR3 netcharge of moderate mutated {iso1}: {mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten}')

            if not stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten.empty and 'Net charge CDR3aa' in stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten.columns:
                mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten = stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                if pd.isna(mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten):
                    mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten = 0
            print(f'Mean CDR3 netcharge of high mutated {iso1}: {mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten}')

            if not stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten.empty and 'Net charge CDR3aa' in stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten.columns:
                mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten = stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                if pd.isna(mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten):
                    mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten = 0
            print(f'Mean CDR3 netcharge of low mutated {iso2}: {mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten}')

            if not stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten.empty and 'Net charge CDR3aa' in stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten.columns:
                mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten = stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                if pd.isna(mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten):
                    mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten = 0
            print(f'Mean CDR3 netcharge of moderate mutated {iso2}: {mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten}')

            if not stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten.empty and 'Net charge CDR3aa' in stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten.columns:
                mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten = stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                if pd.isna(mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten):
                    mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten = 0
            print(f'Mean CDR3 netcharge of high mutated {iso2}: {mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten}')

            # Calculate ratios with safe division
            if mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten != 0:
                ratio_stats_tab_low_mutated_shared_iso2_iso1_VDJ_sorted_topten = mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten / mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten
            else:
                ratio_stats_tab_low_mutated_shared_iso2_iso1_VDJ_sorted_topten = 0
            print(f'low mutated {iso2}/{iso1} netcharge ratio: {ratio_stats_tab_low_mutated_shared_iso2_iso1_VDJ_sorted_topten}')

            if mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten != 0:
                ratio_stats_tab_moderate_mutated_shared_iso2_iso1_VDJ_sorted_topten = mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten / mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten
            else:
                ratio_stats_tab_moderate_mutated_shared_iso2_iso1_VDJ_sorted_topten = 0
            print(f'moderate mutated {iso2}/{iso1} netcharge ratio: {ratio_stats_tab_moderate_mutated_shared_iso2_iso1_VDJ_sorted_topten}')

            if mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten != 0:
                ratio_stats_tab_high_mutated_shared_iso2_iso1_VDJ_sorted_topten = mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten / mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten
            else:
                ratio_stats_tab_high_mutated_shared_iso2_iso1_VDJ_sorted_topten = 0
            print(f'high mutated {iso2}/{iso1} netcharge ratio: {ratio_stats_tab_high_mutated_shared_iso2_iso1_VDJ_sorted_topten}')

            # Compare CDR3aa to Aminoacid tables for divergent analysis
            stats_iso1_Aminotab_productive = stats_iso1_Aminotab[stats_iso1_Aminotab['V-DOMAIN Functionality'] == "productive"]

            # Handle column variations for iso2 Aminotab
            if 'V.DOMAIN.Functionality' in stats_iso2_Aminotab.columns:
                col_functionality = 'V.DOMAIN.Functionality'
            elif 'V-DOMAIN Functionality' in stats_iso2_Aminotab.columns:
                col_functionality = 'V-DOMAIN Functionality'
            else:
                raise KeyError("No functionality column found in iso2 Aminotab")

            stats_iso2_Aminotab_productive = stats_iso2_Aminotab[stats_iso2_Aminotab[col_functionality] == "productive"]

            # Get CDR3 column names
            if 'CDR3.IMGT' in stats_iso2_Aminotab_productive.columns:
                col_cdr3_iso2 = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in stats_iso2_Aminotab_productive.columns:
                col_cdr3_iso2 = 'CDR3-IMGT'
            else:
                raise KeyError("No CDR3 column found")

            stats_iso1_Aminotab_productive_CDR3IMGT = set(stats_iso1_Aminotab_productive['CDR3-IMGT'])
            stats_iso2_Aminotab_productive_CDR3IMGT = set(stats_iso2_Aminotab_productive[col_cdr3_iso2])

            # Compare for divergent clones
            stats_tab_low_mutated_shared_iso1_VDJ_compare = stats_tab_low_mutated_shared_iso1_VDJ.copy()
            stats_tab_moderate_mutated_shared_iso1_VDJ_compare = stats_tab_moderate_mutated_shared_iso1_VDJ.copy()
            stats_tab_high_mutated_shared_iso1_VDJ_compare = stats_tab_high_mutated_shared_iso1_VDJ.copy()

            stats_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_{iso2}_CDR3IMGT'] = stats_tab_low_mutated_shared_iso1_VDJ_compare['cdr3aa'].isin(stats_iso2_Aminotab_productive_CDR3IMGT)
            stats_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_{iso2}_CDR3IMGT'] = stats_tab_moderate_mutated_shared_iso1_VDJ_compare['cdr3aa'].isin(stats_iso2_Aminotab_productive_CDR3IMGT)
            stats_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_{iso2}_CDR3IMGT'] = stats_tab_high_mutated_shared_iso1_VDJ_compare['cdr3aa'].isin(stats_iso2_Aminotab_productive_CDR3IMGT)

            for col in [f'Match_in_low_mutated_{iso2}_CDR3IMGT']:
                stats_tab_low_mutated_shared_iso1_VDJ_compare[col] = stats_tab_low_mutated_shared_iso1_VDJ_compare[col].astype(str).str.strip()
            for col in [f'Match_in_moderate_mutated_{iso2}_CDR3IMGT']:
                stats_tab_moderate_mutated_shared_iso1_VDJ_compare[col] = stats_tab_moderate_mutated_shared_iso1_VDJ_compare[col].astype(str).str.strip()
            for col in [f'Match_in_high_mutated_{iso2}_CDR3IMGT']:
                stats_tab_high_mutated_shared_iso1_VDJ_compare[col] = stats_tab_high_mutated_shared_iso1_VDJ_compare[col].astype(str).str.strip()

            # Filter divergent clones
            divergent_tab_low_mutated_shared_iso1_VDJ = stats_tab_low_mutated_shared_iso1_VDJ_compare[
                stats_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_{iso2}_CDR3IMGT'] == 'False']
            divergent_tab_moderate_mutated_shared_iso1_VDJ = stats_tab_moderate_mutated_shared_iso1_VDJ_compare[
                stats_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_{iso2}_CDR3IMGT'] == 'False']
            divergent_tab_high_mutated_shared_iso1_VDJ = stats_tab_high_mutated_shared_iso1_VDJ_compare[
                stats_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_{iso2}_CDR3IMGT'] == 'False']

            # Filter non-divergent clones
            non_divergent_tab_low_mutated_shared_iso1_VDJ = stats_tab_low_mutated_shared_iso1_VDJ_compare[
                stats_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_{iso2}_CDR3IMGT'] == 'True']
            non_divergent_tab_moderate_mutated_shared_iso1_VDJ = stats_tab_moderate_mutated_shared_iso1_VDJ_compare[
                stats_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_{iso2}_CDR3IMGT'] == 'True']
            non_divergent_tab_high_mutated_shared_iso1_VDJ = stats_tab_high_mutated_shared_iso1_VDJ_compare[
                stats_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_{iso2}_CDR3IMGT'] == 'True']

            # Create seqid sets for divergent analysis
            stats_iso1_low_mutated_seqid_CDR3_2 = set(divergent_tab_low_mutated_shared_iso1_VDJ['seqid'])
            stats_iso1_moderate_mutated_seqid_CDR3_2 = set(divergent_tab_moderate_mutated_shared_iso1_VDJ['seqid'])
            stats_iso1_high_mutated_seqid_CDR3_2 = set(divergent_tab_high_mutated_shared_iso1_VDJ['seqid'])

            non_divergent_stats_iso1_low_mutated_seqid_CDR3_2 = set(non_divergent_tab_low_mutated_shared_iso1_VDJ['seqid'])
            non_divergent_stats_iso1_moderate_mutated_seqid_CDR3_2 = set(non_divergent_tab_moderate_mutated_shared_iso1_VDJ['seqid'])
            non_divergent_stats_iso1_high_mutated_seqid_CDR3_2 = set(non_divergent_tab_high_mutated_shared_iso1_VDJ['seqid'])

            # Filter divergent with CDR3 mutations >= 2
            stats_iso1_mut_3_mut_level_compare = stats_iso1_mut_3.copy()

            stats_iso1_mut_3_mut_level_compare[f'divergent_low_mutated_shared_{iso1}_seqid'] = stats_iso1_mut_3_mut_level_compare['Sequence ID'].isin(stats_iso1_low_mutated_seqid_CDR3_2)
            stats_iso1_mut_3_mut_level_compare[f'divergent_moderate_mutated_shared_{iso1}_seqid'] = stats_iso1_mut_3_mut_level_compare['Sequence ID'].isin(stats_iso1_moderate_mutated_seqid_CDR3_2)
            stats_iso1_mut_3_mut_level_compare[f'divergent_high_mutated_shared_{iso1}_seqid'] = stats_iso1_mut_3_mut_level_compare['Sequence ID'].isin(stats_iso1_high_mutated_seqid_CDR3_2)

            stats_iso1_mut_3_mut_level_compare[f'non_divergent_low_mutated_shared_{iso1}_seqid'] = stats_iso1_mut_3_mut_level_compare['Sequence ID'].isin(non_divergent_stats_iso1_low_mutated_seqid_CDR3_2)
            stats_iso1_mut_3_mut_level_compare[f'non_divergent_moderate_mutated_shared_{iso1}_seqid'] = stats_iso1_mut_3_mut_level_compare['Sequence ID'].isin(non_divergent_stats_iso1_moderate_mutated_seqid_CDR3_2)
            stats_iso1_mut_3_mut_level_compare[f'non_divergent_high_mutated_shared_{iso1}_seqid'] = stats_iso1_mut_3_mut_level_compare['Sequence ID'].isin(non_divergent_stats_iso1_high_mutated_seqid_CDR3_2)

            for col in [f'divergent_low_mutated_shared_{iso1}_seqid', f'divergent_moderate_mutated_shared_{iso1}_seqid',
                        f'divergent_high_mutated_shared_{iso1}_seqid', f'non_divergent_low_mutated_shared_{iso1}_seqid',
                        f'non_divergent_moderate_mutated_shared_{iso1}_seqid', f'non_divergent_high_mutated_shared_{iso1}_seqid']:
                stats_iso1_mut_3_mut_level_compare[col] = stats_iso1_mut_3_mut_level_compare[col].astype(str).str.strip()

            # Filter divergent and non-divergent mutation tables with >= 2 CDR3 mutations
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2 = stats_iso1_mut_3_mut_level_compare[
                (stats_iso1_mut_3_mut_level_compare[f'divergent_low_mutated_shared_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3_thresh)]
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2 = stats_iso1_mut_3_mut_level_compare[
                (stats_iso1_mut_3_mut_level_compare[f'divergent_moderate_mutated_shared_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3_thresh)]
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2 = stats_iso1_mut_3_mut_level_compare[
                (stats_iso1_mut_3_mut_level_compare[f'divergent_high_mutated_shared_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3_thresh)]

            non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2 = stats_iso1_mut_3_mut_level_compare[
                (stats_iso1_mut_3_mut_level_compare[f'non_divergent_low_mutated_shared_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3_thresh)]
            non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2 = stats_iso1_mut_3_mut_level_compare[
                (stats_iso1_mut_3_mut_level_compare[f'non_divergent_moderate_mutated_shared_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3_thresh)]
            non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2 = stats_iso1_mut_3_mut_level_compare[
                (stats_iso1_mut_3_mut_level_compare[f'non_divergent_high_mutated_shared_{iso1}_seqid'] == 'True') &
                (stats_iso1_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3_thresh)]

            # Create seqid sets
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2_seqid = set(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2['Sequence ID'])
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2_seqid = set(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2['Sequence ID'])
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2_seqid = set(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2['Sequence ID'])

            non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2_seqid = set(non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2['Sequence ID'])
            non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2_seqid = set(non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2['Sequence ID'])
            non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2_seqid = set(non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2['Sequence ID'])

            # Update comparison DataFrames
            divergent_tab_low_mutated_shared_iso1_VDJ_compare = divergent_tab_low_mutated_shared_iso1_VDJ.copy()
            divergent_tab_moderate_mutated_shared_iso1_VDJ_compare = divergent_tab_moderate_mutated_shared_iso1_VDJ.copy()
            divergent_tab_high_mutated_shared_iso1_VDJ_compare = divergent_tab_high_mutated_shared_iso1_VDJ.copy()

            non_divergent_tab_low_mutated_shared_iso1_VDJ_compare = non_divergent_tab_low_mutated_shared_iso1_VDJ.copy()
            non_divergent_tab_moderate_mutated_shared_iso1_VDJ_compare = non_divergent_tab_moderate_mutated_shared_iso1_VDJ.copy()
            non_divergent_tab_high_mutated_shared_iso1_VDJ_compare = non_divergent_tab_high_mutated_shared_iso1_VDJ.copy()

            divergent_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] = divergent_tab_low_mutated_shared_iso1_VDJ_compare['seqid'].isin(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            divergent_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] = divergent_tab_moderate_mutated_shared_iso1_VDJ_compare['seqid'].isin(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            divergent_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] = divergent_tab_high_mutated_shared_iso1_VDJ_compare['seqid'].isin(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2_seqid)

            non_divergent_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] = non_divergent_tab_low_mutated_shared_iso1_VDJ_compare['seqid'].isin(non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            non_divergent_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] = non_divergent_tab_moderate_mutated_shared_iso1_VDJ_compare['seqid'].isin(non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            non_divergent_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] = non_divergent_tab_high_mutated_shared_iso1_VDJ_compare['seqid'].isin(non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2_seqid)

            for df in [divergent_tab_low_mutated_shared_iso1_VDJ_compare, divergent_tab_moderate_mutated_shared_iso1_VDJ_compare,
                       divergent_tab_high_mutated_shared_iso1_VDJ_compare, non_divergent_tab_low_mutated_shared_iso1_VDJ_compare,
                       non_divergent_tab_moderate_mutated_shared_iso1_VDJ_compare, non_divergent_tab_high_mutated_shared_iso1_VDJ_compare]:
                for col in df.columns:
                    if 'Match_in' in col and 'CDR3_2_seqid' in col:
                        df[col] = df[col].astype(str).str.strip()

            # Filter divergent clones
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_compare[
                divergent_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] == 'True']
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_compare[
                divergent_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] == 'True']
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_compare[
                divergent_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] == 'True']

            non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_iso1_VDJ_compare[
                non_divergent_tab_low_mutated_shared_iso1_VDJ_compare[f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] == 'True']
            non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_iso1_VDJ_compare[
                non_divergent_tab_moderate_mutated_shared_iso1_VDJ_compare[f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] == 'True']
            non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_iso1_VDJ_compare[
                non_divergent_tab_high_mutated_shared_iso1_VDJ_compare[f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3_2_seqid'] == 'True']

            # Calculate unique clone counts
            unique_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3 low mutated {iso1} clones: {unique_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            unique_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3 moderate mutated {iso1} clones: {unique_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            unique_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3 high mutated {iso1} clones: {unique_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            unique_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of non_divergent CDR3 low mutated {iso1} clones: {unique_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            unique_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of non_divergent CDR3 moderate mutated {iso1} clones: {unique_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            unique_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of non_divergent CDR3 high mutated {iso1} clones: {unique_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            # Calculate mean copy numbers for divergent and non-divergent clones
            mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone) else mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'mean copynumber of divergent CDR3 low mutated {iso1}: {mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone) else mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'mean copynumber of divergent CDR3 moderate mutated {iso1}: {mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone) else mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'mean copynumber of divergent CDR3 high mutated {iso1}: {mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            mean_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(mean_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone) else mean_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'mean copynumber of non_divergent CDR3 low mutated {iso1}: {mean_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            mean_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(mean_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone) else mean_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'mean copynumber of non_divergent CDR3 moderate mutated {iso1}: {mean_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            mean_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(mean_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone) else mean_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'mean copynumber of non_divergent CDR3 high mutated {iso1}: {mean_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            # Calculate sum copy numbers
            sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone) else sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'sum copynumber of divergent CDR3 low mutated {iso1}: {sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone) else sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'sum copynumber of divergent CDR3 moderate mutated {iso1}: {sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone) else sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'sum copynumber of divergent CDR3 high mutated {iso1}: {sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            sum_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(sum_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone) else sum_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'sum copynumber of non_divergent CDR3 low mutated {iso1}: {sum_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            sum_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(sum_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone) else sum_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'sum copynumber of non_divergent CDR3 moderate mutated {iso1}: {sum_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            sum_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone = 0 if pd.isna(sum_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone) else sum_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone
            print(f'sum copynumber of non_divergent CDR3 high mutated {iso1}: {sum_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')

            # Find clonally related iso2 counterparts
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone_seqid = set(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'])
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone_seqid = set(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'])
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone_seqid = set(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone[f'{iso1}_VDJ'])

            clonally_related_low_mod_high_shared_iso2 = shared_iso2.copy()

            clonally_related_low_mod_high_shared_iso2[f'divergent_low_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid'] = clonally_related_low_mod_high_shared_iso2[f'{iso2}_VDJ'].isin(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone_seqid)
            clonally_related_low_mod_high_shared_iso2[f'divergent_moderate_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid'] = clonally_related_low_mod_high_shared_iso2[f'{iso2}_VDJ'].isin(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone_seqid)
            clonally_related_low_mod_high_shared_iso2[f'divergent_high_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid'] = clonally_related_low_mod_high_shared_iso2[f'{iso2}_VDJ'].isin(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone_seqid)

            for col in [f'divergent_low_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid', 
                        f'divergent_moderate_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid',
                        f'divergent_high_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid']:
                clonally_related_low_mod_high_shared_iso2[col] = clonally_related_low_mod_high_shared_iso2[col].astype(str).str.strip()

            clonally_related_low_mutated_shared_iso2 = clonally_related_low_mod_high_shared_iso2[
                clonally_related_low_mod_high_shared_iso2[f'divergent_low_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid'] == 'True']
            clonally_related_moderate_mutated_shared_iso2 = clonally_related_low_mod_high_shared_iso2[
                clonally_related_low_mod_high_shared_iso2[f'divergent_moderate_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid'] == 'True']
            clonally_related_high_mutated_shared_iso2 = clonally_related_low_mod_high_shared_iso2[
                clonally_related_low_mod_high_shared_iso2[f'divergent_high_mutated_shared_{iso1}_VDJ_in_{iso2}_seqid'] == 'True']

            # Calculate mean copy numbers for clonally related
            mean_clonally_related_low_mutated_shared_iso2 = clonally_related_low_mutated_shared_iso2['onecopy'].mean()
            mean_clonally_related_low_mutated_shared_iso2 = 0 if pd.isna(mean_clonally_related_low_mutated_shared_iso2) else mean_clonally_related_low_mutated_shared_iso2
            print(f'mean copynumber of {iso2} clonally related to low mutated {iso1}: {mean_clonally_related_low_mutated_shared_iso2}')

            mean_clonally_related_moderate_mutated_shared_iso2 = clonally_related_moderate_mutated_shared_iso2['onecopy'].mean()
            mean_clonally_related_moderate_mutated_shared_iso2 = 0 if pd.isna(mean_clonally_related_moderate_mutated_shared_iso2) else mean_clonally_related_moderate_mutated_shared_iso2
            print(f'mean copynumber of {iso2} clonally related to moderate mutated {iso1}: {mean_clonally_related_moderate_mutated_shared_iso2}')

            mean_clonally_related_high_mutated_shared_iso2 = clonally_related_high_mutated_shared_iso2['onecopy'].mean()
            mean_clonally_related_high_mutated_shared_iso2 = 0 if pd.isna(mean_clonally_related_high_mutated_shared_iso2) else mean_clonally_related_high_mutated_shared_iso2
            print(f'mean copynumber of {iso2} clonally related to high mutated {iso1}: {mean_clonally_related_high_mutated_shared_iso2}')

            # Calculate sum copy numbers for clonally related
            sum_clonally_related_low_mutated_shared_iso2 = clonally_related_low_mutated_shared_iso2['onecopy'].sum()
            sum_clonally_related_low_mutated_shared_iso2 = 0 if pd.isna(sum_clonally_related_low_mutated_shared_iso2) else sum_clonally_related_low_mutated_shared_iso2
            print(f'sum copynumber of {iso2} clonally related to low mutated {iso1}: {sum_clonally_related_low_mutated_shared_iso2}')

            sum_clonally_related_moderate_mutated_shared_iso2 = clonally_related_moderate_mutated_shared_iso2['onecopy'].sum()
            sum_clonally_related_moderate_mutated_shared_iso2 = 0 if pd.isna(sum_clonally_related_moderate_mutated_shared_iso2) else sum_clonally_related_moderate_mutated_shared_iso2
            print(f'sum copynumber of {iso2} clonally related to moderate mutated {iso1}: {sum_clonally_related_moderate_mutated_shared_iso2}')

            sum_clonally_related_high_mutated_shared_iso2 = clonally_related_high_mutated_shared_iso2['onecopy'].sum()
            sum_clonally_related_high_mutated_shared_iso2 = 0 if pd.isna(sum_clonally_related_high_mutated_shared_iso2) else sum_clonally_related_high_mutated_shared_iso2
            print(f'sum copynumber of {iso2} clonally related to high mutated {iso1}: {sum_clonally_related_high_mutated_shared_iso2}')
            self.logger.log_memory("After divergent clone analysis")



            # ============================================
            # CDR3CDR2 DIVERGENT ANALYSIS SECTION
            # ============================================
            
            # Analyse the low, moderate and highly divergent based on the CDR3-CDR2 aminoacid sequence
            # Merge the CDR3 and CDR2 aminoacid sequence of Total Isotype 1
            stats_iso1_Aminotab_productive_CDR3CDR2 = stats_iso1_Aminotab_productive.copy()
            stats_iso1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['CDR2-IMGT'].astype(str) + stats_iso1_Aminotab_productive_CDR3CDR2['CDR3-IMGT'].astype(str)
            stats_iso1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].str.replace(' ', '', regex=False)
            
            # Merge the CDR3 and CDR2 aminoacid sequence of Total Isotype 2
            stats_iso2_Aminotab_productive_CDR3CDR2 = stats_iso2_Aminotab_productive.copy()
            
            if 'CDR2.IMGT' in stats_iso2_Aminotab_productive_CDR3CDR2.columns:
                column_cdr2_iso2_aa = 'CDR2.IMGT'
            elif 'CDR2-IMGT' in stats_iso2_Aminotab_productive_CDR3CDR2.columns:
                column_cdr2_iso2_aa = 'CDR2-IMGT'
            else:
                raise KeyError("Neither 'CDR2.IMGT' nor 'CDR2-IMGT' column is present in the DataFrame.")
            
            if 'CDR3.IMGT' in stats_iso2_Aminotab_productive_CDR3CDR2.columns:
                column_cdr3_iso2_aa = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in stats_iso2_Aminotab_productive_CDR3CDR2.columns:
                column_cdr3_iso2_aa = 'CDR3-IMGT'
            else:
                raise KeyError("Neither 'CDR3.IMGT' nor 'CDR3-IMGT' column is present in the DataFrame.")
            
            stats_iso2_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'] = \
                stats_iso2_Aminotab_productive_CDR3CDR2[column_cdr2_iso2_aa].astype(str) + stats_iso2_Aminotab_productive_CDR3CDR2[column_cdr3_iso2_aa].astype(str)
            stats_iso2_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'] = stats_iso2_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].str.replace(' ', '', regex=False)
            
            # Extract the high, moderate and low mutated CDR3CDR2
            stats_iso1_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['Sequence ID'].isin(low_mutated_seqids_Shared_iso1)
            stats_iso1_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['Sequence ID'].isin(moderate_mutated_seqids_Shared_iso1)
            stats_iso1_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['Sequence ID'].isin(High_mutated_seqids_Shared_iso1)
            
            stats_iso1_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-IMGT'].astype(str).str.strip()
            stats_iso1_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'].astype(str).str.strip()
            stats_iso1_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-IMGT'] = stats_iso1_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-IMGT'].astype(str).str.strip()
            
            stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_iso1_Aminotab_productive_CDR3CDR2[(stats_iso1_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-IMGT'] == 'True')]
            stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_iso1_Aminotab_productive_CDR3CDR2[(stats_iso1_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'] == 'True')]
            stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_iso1_Aminotab_productive_CDR3CDR2[(stats_iso1_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-IMGT'] == 'True')]
            
            # Compare the iso2 CDR3CDR2 aminoacid table to the mutated CDR3CDR2
            stats_iso2_Aminotab_productive_CDR3CDR2_IMGTAA = set(stats_iso2_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'])
            
            stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2.copy()
            stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2.copy()
            stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2.copy()
            
            stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2[f'low_mutated_CDR3CDR2-in_{iso2}_CDR3CDR2'] = stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].isin(stats_iso2_Aminotab_productive_CDR3CDR2_IMGTAA)
            stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2[f'moderate_mutated_CDR3CDR2-in_{iso2}_CDR3CDR2'] = stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].isin(stats_iso2_Aminotab_productive_CDR3CDR2_IMGTAA)
            stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2[f'high_mutated_CDR3CDR2-in_{iso2}_CDR3CDR2'] = stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].isin(stats_iso2_Aminotab_productive_CDR3CDR2_IMGTAA)
            
            divergent_stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2[~stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2[f'low_mutated_CDR3CDR2-in_{iso2}_CDR3CDR2']]
            divergent_stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2[~stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2[f'moderate_mutated_CDR3CDR2-in_{iso2}_CDR3CDR2']]
            divergent_stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2 = stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2[~stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2[f'high_mutated_CDR3CDR2-in_{iso2}_CDR3CDR2']]
            
            divergent_stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2_seqid = set(divergent_stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2['Sequence ID'])
            divergent_stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2_seqid = set(divergent_stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2['Sequence ID'])
            divergent_stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2_seqid = set(divergent_stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2['Sequence ID'])
            
            stats_iso1_mut_3_CDR3CDR2_CDR3_2 = stats_iso1_mut_3.copy()
            
            stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_low_mutated_{iso1}_seqid'] = stats_iso1_mut_3_CDR3CDR2_CDR3_2['Sequence ID'].isin(divergent_stats_low_mutated_iso1_Aminotab_productive_CDR3CDR2_seqid)
            stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_moderate_mutated_{iso1}_seqid'] = stats_iso1_mut_3_CDR3CDR2_CDR3_2['Sequence ID'].isin(divergent_stats_moderate_mutated_iso1_Aminotab_productive_CDR3CDR2_seqid)
            stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_high_mutated_{iso1}_seqid'] = stats_iso1_mut_3_CDR3CDR2_CDR3_2['Sequence ID'].isin(divergent_stats_high_mutated_iso1_Aminotab_productive_CDR3CDR2_seqid)
            
            stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_low_mutated_{iso1}_seqid'] = stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_low_mutated_{iso1}_seqid'].astype(str).str.strip()
            stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_moderate_mutated_{iso1}_seqid'] = stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_moderate_mutated_{iso1}_seqid'].astype(str).str.strip()
            stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_high_mutated_{iso1}_seqid'] = stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_high_mutated_{iso1}_seqid'].astype(str).str.strip()
            
            divergent_tab_low_mutated_shared_iso1_CDR3CDR2_mutab = stats_iso1_mut_3_CDR3CDR2_CDR3_2[(stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_low_mutated_{iso1}_seqid'] == 'True') & (stats_iso1_mut_3_CDR3CDR2_CDR3_2['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3cdr2_thresh)]
            divergent_tab_moderate_mutated_shared_iso1_CDR3CDR2_mutab = stats_iso1_mut_3_CDR3CDR2_CDR3_2[(stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_moderate_mutated_{iso1}_seqid'] == 'True') & (stats_iso1_mut_3_CDR3CDR2_CDR3_2['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3cdr2_thresh)]
            divergent_tab_high_mutated_shared_iso1_CDR3CDR2_mutab = stats_iso1_mut_3_CDR3CDR2_CDR3_2[(stats_iso1_mut_3_CDR3CDR2_CDR3_2[f'Match_in_divergent_high_mutated_{iso1}_seqid'] == 'True') & (stats_iso1_mut_3_CDR3CDR2_CDR3_2['CDR3-IMGT Nb of nonsilent mutations-ex'] >= cdr3cdr2_thresh)]
            
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2_2_seqid = set(divergent_tab_low_mutated_shared_iso1_CDR3CDR2_mutab['Sequence ID'])
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2_2_seqid = set(divergent_tab_moderate_mutated_shared_iso1_CDR3CDR2_mutab['Sequence ID'])
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2_2_seqid = set(divergent_tab_high_mutated_shared_iso1_CDR3CDR2_mutab['Sequence ID'])
            
            divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare = shared_iso1_main.copy()
            
            divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare['seqid'].isin(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2_2_seqid)
            divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare['seqid'].isin(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2_2_seqid)
            divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare['seqid'].isin(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2_2_seqid)
            
            for col in [f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid', f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid', f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid']:
                divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[col] = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[col].astype(str).str.strip()
            
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[(divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[f'Match_in_low_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid'] == 'True')]
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[(divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[f'Match_in_moderate_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid'] == 'True')]
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[(divergent_tab_low_moderate_high_mutated_shared_iso1_VDJ_CDR3CDR2_compare[f'Match_in_high_mutated_shared_{iso1}_VDJ_CDR3CDR2_2_seqid'] == 'True')]
            
            unique_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3CDR2 low mutated {iso1} clones: {unique_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            unique_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3CDR2 moderate mutated {iso1} clones: {unique_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            unique_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone[f'{iso1}_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3CDR2 high mutated {iso1} clones: {unique_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone['onecopy'].mean()
            mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone) else mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone
            print(f'mean copynumber of divergent CDR3CDR2 low mutated {iso1}: {mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone['onecopy'].mean()
            mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone) else mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone
            print(f'mean copynumber of divergent CDR3CDR2 moderate mutated {iso1}: {mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone['onecopy'].mean()
            mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone) else mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone
            print(f'mean copynumber of divergent CDR3CDR2 high mutated {iso1}: {mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone['onecopy'].sum()
            sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone) else sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone
            print(f'sum copynumber of divergent CDR3CDR2 low mutated {iso1}: {sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone['onecopy'].sum()
            sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone) else sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone
            print(f'sum copynumber of divergent CDR3CDR2 moderate mutated {iso1}: {sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone['onecopy'].sum()
            sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone) else sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone
            print(f'sum copynumber of divergent CDR3CDR2 high mutated {iso1}: {sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            # Filter clonally related CDR3CDR2 iso2 counterparts
            divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone_seqid = set(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone[f'{iso1}_VDJ'])
            divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone_seqid = set(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone[f'{iso1}_VDJ'])
            divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone_seqid = set(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone[f'{iso1}_VDJ'])
            
            clonally_related_low_moderate_high_shared_CDR3CDR2_iso2 = shared_iso2_main.copy()
            
            clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'divergent_low_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'{iso2}_VDJ'].isin(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone_seqid)
            clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'divergent_moderate_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'{iso2}_VDJ'].isin(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone_seqid)
            clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'divergent_high_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'{iso2}_VDJ'].isin(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone_seqid)
            
            for col in [f'divergent_low_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid', f'divergent_moderate_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid', f'divergent_high_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid']:
                clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[col] = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[col].astype(str).str.strip()
            
            clonally_related_low_mutated_shared_CDR3CDR2_iso2 = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[(clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'divergent_low_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid'] == 'True')]
            clonally_related_moderate_mutated_shared_CDR3CDR2_iso2 = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[(clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'divergent_moderate_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid'] == 'True')]
            clonally_related_high_mutated_shared_CDR3CDR2_iso2 = clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[(clonally_related_low_moderate_high_shared_CDR3CDR2_iso2[f'divergent_high_mutated_shared_{iso1}_VDJ_CDR3CDR2_in_{iso2}_seqid'] == 'True')]
            
            mean_clonally_related_low_mutated_shared_CDR3CDR2_iso2 = clonally_related_low_mutated_shared_CDR3CDR2_iso2['onecopy'].mean()
            mean_clonally_related_low_mutated_shared_CDR3CDR2_iso2 = 0 if pd.isna(mean_clonally_related_low_mutated_shared_CDR3CDR2_iso2) else mean_clonally_related_low_mutated_shared_CDR3CDR2_iso2
            print(f'mean copynumber of {iso2} clonally related to low mutated {iso1} (CDR3CDR2): {mean_clonally_related_low_mutated_shared_CDR3CDR2_iso2}')
            
            mean_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2 = clonally_related_moderate_mutated_shared_CDR3CDR2_iso2['onecopy'].mean()
            mean_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2 = 0 if pd.isna(mean_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2) else mean_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2
            print(f'mean copynumber of {iso2} clonally related to moderate mutated {iso1} (CDR3CDR2): {mean_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2}')
            
            mean_clonally_related_high_mutated_shared_CDR3CDR2_iso2 = clonally_related_high_mutated_shared_CDR3CDR2_iso2['onecopy'].mean()
            mean_clonally_related_high_mutated_shared_CDR3CDR2_iso2 = 0 if pd.isna(mean_clonally_related_high_mutated_shared_CDR3CDR2_iso2) else mean_clonally_related_high_mutated_shared_CDR3CDR2_iso2
            print(f'mean copynumber of {iso2} clonally related to high mutated {iso1} (CDR3CDR2): {mean_clonally_related_high_mutated_shared_CDR3CDR2_iso2}')


            sum_clonally_related_low_mutated_shared_CDR3CDR2_iso2 = clonally_related_low_mutated_shared_CDR3CDR2_iso2['onecopy'].sum()
            sum_clonally_related_low_mutated_shared_CDR3CDR2_iso2 = 0 if pd.isna(sum_clonally_related_low_mutated_shared_CDR3CDR2_iso2) else sum_clonally_related_low_mutated_shared_CDR3CDR2_iso2
            print(f'sum copynumber of {iso2} clonally related to low mutated {iso1} (CDR3CDR2): {sum_clonally_related_low_mutated_shared_CDR3CDR2_iso2}')
            
            sum_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2 = clonally_related_moderate_mutated_shared_CDR3CDR2_iso2['onecopy'].sum()
            sum_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2 = 0 if pd.isna(sum_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2) else sum_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2
            print(f'sum copynumber of {iso2} clonally related to moderate mutated {iso1} (CDR3CDR2): {sum_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2}')
            
            sum_clonally_related_high_mutated_shared_CDR3CDR2_iso2 = clonally_related_high_mutated_shared_CDR3CDR2_iso2['onecopy'].sum()
            sum_clonally_related_high_mutated_shared_CDR3CDR2_iso2 = 0 if pd.isna(sum_clonally_related_high_mutated_shared_CDR3CDR2_iso2) else sum_clonally_related_high_mutated_shared_CDR3CDR2_iso2
            print(f'sum copynumber of {iso2} clonally related to high mutated {iso1} (CDR3CDR2): {sum_clonally_related_high_mutated_shared_CDR3CDR2_iso2}')


            
            # ============================================
            # VH SEQUENCE EXTRACTION SECTION
            # ============================================
            
            # Filter out the entire VH gene of the divergent low, moderate and highly mutated iso1 using CDR3 AA change
            stats_iso1_Aminotab_productive_CDR3mut2 = stats_iso1_Aminotab_productive.copy()
            
            stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_low_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2['Sequence ID'].isin(divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_moderate_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2['Sequence ID'].isin(divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_high_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2['Sequence ID'].isin(divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            
            stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_low_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_low_mutated_shared_{iso1}_VH_AA'].astype(str).str.strip()
            stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_moderate_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_moderate_mutated_shared_{iso1}_VH_AA'].astype(str).str.strip()
            stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_high_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_high_mutated_shared_{iso1}_VH_AA'].astype(str).str.strip()
            
            divergent_low_mutated_shared_iso1_VH_AA_tab = stats_iso1_Aminotab_productive_CDR3mut2[(stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_low_mutated_shared_{iso1}_VH_AA'] == 'True')]
            divergent_moderate_mutated_shared_iso1_VH_AA_tab = stats_iso1_Aminotab_productive_CDR3mut2[(stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_moderate_mutated_shared_{iso1}_VH_AA'] == 'True')]
            divergent_high_mutated_shared_iso1_VH_AA_tab = stats_iso1_Aminotab_productive_CDR3mut2[(stats_iso1_Aminotab_productive_CDR3mut2[f'divergent_high_mutated_shared_{iso1}_VH_AA'] == 'True')]
            
            # Extract VH sequences
            divergent_low_mutated_VH_AAsequences = divergent_low_mutated_shared_iso1_VH_AA_tab["V-D-J-REGION"].tolist()
            for i, sequence in enumerate(divergent_low_mutated_VH_AAsequences, 1):
                print(f"Divergent low mutated {iso1} Row {i}: {sequence}")
            
            divergent_moderate_mutated_VH_AAsequences = divergent_moderate_mutated_shared_iso1_VH_AA_tab["V-D-J-REGION"].tolist()
            for i, sequence in enumerate(divergent_moderate_mutated_VH_AAsequences, 1):
                print(f"Divergent moderate mutated {iso1} Row {i}: {sequence}")
            
            divergent_high_mutated_VH_AAsequences = divergent_high_mutated_shared_iso1_VH_AA_tab["V-D-J-REGION"].tolist()
            for i, sequence in enumerate(divergent_high_mutated_VH_AAsequences, 1):
                print(f"Divergent high mutated {iso1} Row {i}: {sequence}")
            
            # Filter out the entire VH gene of the clonally related iso2
            clonally_related_low_mutated_shared_iso2_seqid = set(clonally_related_low_mutated_shared_iso2['seqid'])
            clonally_related_moderate_mutated_shared_iso2_seqid = set(clonally_related_moderate_mutated_shared_iso2['seqid'])
            clonally_related_high_mutated_shared_iso2_seqid = set(clonally_related_high_mutated_shared_iso2['seqid'])
            
            stats_iso2_Aminotab_productive_CDR3mut2 = stats_iso2_Aminotab_productive.copy()
            
            if 'Sequence.ID' in stats_iso2_Aminotab_productive_CDR3mut2.columns:
                column_seqid_iso2_aa = 'Sequence.ID'
            elif 'Sequence ID' in stats_iso2_Aminotab_productive_CDR3mut2.columns:
                column_seqid_iso2_aa = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_low_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] = stats_iso2_Aminotab_productive_CDR3mut2[column_seqid_iso2_aa].isin(clonally_related_low_mutated_shared_iso2_seqid)
            stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_moderate_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] = stats_iso2_Aminotab_productive_CDR3mut2[column_seqid_iso2_aa].isin(clonally_related_moderate_mutated_shared_iso2_seqid)
            stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_high_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] = stats_iso2_Aminotab_productive_CDR3mut2[column_seqid_iso2_aa].isin(clonally_related_high_mutated_shared_iso2_seqid)
            
            stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_low_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] = stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_low_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'].astype(str).str.strip()
            stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_moderate_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] = stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_moderate_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'].astype(str).str.strip()
            stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_high_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] = stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_high_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'].astype(str).str.strip()
            
            clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AA = stats_iso2_Aminotab_productive_CDR3mut2[(stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_low_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] == 'True')]
            clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AA = stats_iso2_Aminotab_productive_CDR3mut2[(stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_moderate_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] == 'True')]
            clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AA = stats_iso2_Aminotab_productive_CDR3mut2[(stats_iso2_Aminotab_productive_CDR3mut2[f'clonally_related_high_mutated_shared_{iso2}_seqid_in_{iso2}_seqid'] == 'True')]
            
            # Determine VDJ region column for iso2
            if 'V.D.J.REGION' in clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AA.columns:
                column_vdj_region_iso2 = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AA.columns:
                column_vdj_region_iso2 = 'V-D-J-REGION'
            else:
                column_vdj_region_iso2 = 'V-D-J-REGION'  # Default
            
            clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AAsequences = clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AA[column_vdj_region_iso2].tolist() if column_vdj_region_iso2 in clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AA.columns else []
            for i, sequence in enumerate(clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AAsequences, 1):
                print(f"Clonally related low mutated {iso2} Row {i}: {sequence}")
            
            clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AAsequences = clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AA[column_vdj_region_iso2].tolist() if column_vdj_region_iso2 in clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AA.columns else []
            for i, sequence in enumerate(clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AAsequences, 1):
                print(f"Clonally related moderate mutated {iso2} Row {i}: {sequence}")
            
            clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AAsequences = clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AA[column_vdj_region_iso2].tolist() if column_vdj_region_iso2 in clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AA.columns else []
            for i, sequence in enumerate(clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AAsequences, 1):
                print(f"Clonally related high mutated {iso2} Row {i}: {sequence}")

            # Non-divergent VH extraction
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv = stats_iso1_Aminotab_productive.copy()
            
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_low_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2_nondiv['Sequence ID'].isin(non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_moderate_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2_nondiv['Sequence ID'].isin(non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_high_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2_nondiv['Sequence ID'].isin(non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3_2_seqid)
            
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_low_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_low_mutated_shared_{iso1}_VH_AA'].astype(str).str.strip()
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_moderate_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_moderate_mutated_shared_{iso1}_VH_AA'].astype(str).str.strip()
            stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_high_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_high_mutated_shared_{iso1}_VH_AA'].astype(str).str.strip()
            
            non_divergent_low_mutated_shared_iso1_VH_AA_tab = stats_iso1_Aminotab_productive_CDR3mut2_nondiv[(stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_low_mutated_shared_{iso1}_VH_AA'] == 'True')]
            non_divergent_moderate_mutated_shared_iso1_VH_AA_tab = stats_iso1_Aminotab_productive_CDR3mut2_nondiv[(stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_moderate_mutated_shared_{iso1}_VH_AA'] == 'True')]
            non_divergent_high_mutated_shared_iso1_VH_AA_tab = stats_iso1_Aminotab_productive_CDR3mut2_nondiv[(stats_iso1_Aminotab_productive_CDR3mut2_nondiv[f'non_divergent_high_mutated_shared_{iso1}_VH_AA'] == 'True')]
            
            non_divergent_low_mutated_VH_AAsequences = non_divergent_low_mutated_shared_iso1_VH_AA_tab["V-D-J-REGION"].tolist()
            for i, sequence in enumerate(non_divergent_low_mutated_VH_AAsequences, 1):
                print(f"Non-divergent low mutated {iso1} Row {i}: {sequence}")
            
            non_divergent_moderate_mutated_VH_AAsequences = non_divergent_moderate_mutated_shared_iso1_VH_AA_tab["V-D-J-REGION"].tolist()
            for i, sequence in enumerate(non_divergent_moderate_mutated_VH_AAsequences, 1):
                print(f"Non-divergent moderate mutated {iso1} Row {i}: {sequence}")
            
            non_divergent_high_mutated_VH_AAsequences = non_divergent_high_mutated_shared_iso1_VH_AA_tab["V-D-J-REGION"].tolist()
            for i, sequence in enumerate(non_divergent_high_mutated_VH_AAsequences, 1):
                print(f"Non-divergent high mutated {iso1} Row {i}: {sequence}")

            # Top ten VH extraction for iso1
            stats_iso1_Aminotab_productive_topten = stats_iso1_Aminotab_productive.copy()
            
            stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_seqid = set(stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten['seqid'])
            stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_seqid = set(stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten['seqid'])
            stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_seqid = set(stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten['seqid'])
            
            stats_iso1_Aminotab_productive_topten[f'topten_low_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_topten['Sequence ID'].isin(stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_seqid)
            stats_iso1_Aminotab_productive_topten[f'topten_moderate_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_topten['Sequence ID'].isin(stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_seqid)
            stats_iso1_Aminotab_productive_topten[f'topten_high_mutated_shared_{iso1}_VH_AA'] = stats_iso1_Aminotab_productive_topten['Sequence ID'].isin(stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_seqid)
            
            stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_tab = stats_iso1_Aminotab_productive_topten[stats_iso1_Aminotab_productive_topten[f'topten_low_mutated_shared_{iso1}_VH_AA']]
            stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_tab = stats_iso1_Aminotab_productive_topten[stats_iso1_Aminotab_productive_topten[f'topten_moderate_mutated_shared_{iso1}_VH_AA']]
            stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_tab = stats_iso1_Aminotab_productive_topten[stats_iso1_Aminotab_productive_topten[f'topten_high_mutated_shared_{iso1}_VH_AA']]
            
            stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_AA = stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_tab["V-D-J-REGION"].tolist()
            stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_AA = stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_tab["V-D-J-REGION"].tolist()
            stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_AA = stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_tab["V-D-J-REGION"].tolist()
            
            print(f"\n top ten Low Mutated Shared {iso1} Sequences:")
            for i, sequence in enumerate(stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print(f"\n top ten Moderate Mutated Shared {iso1} Sequences:")
            for i, sequence in enumerate(stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print(f"\n top ten High Mutated Shared {iso1} Sequences:")
            for i, sequence in enumerate(stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")

            # Top ten VH extraction for iso2
            stats_iso2_Aminotab_productive_topten = stats_iso2_Aminotab_productive.copy()
            
            stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_seqid = set(stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten['seqid'])
            stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_seqid = set(stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten['seqid'])
            stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_seqid = set(stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten['seqid'])
            
            stats_iso2_Aminotab_productive_topten[f'topten_low_mutated_shared_{iso2}_VH_AA'] = stats_iso2_Aminotab_productive_topten[column_seqid_iso2_aa].isin(stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_seqid)
            stats_iso2_Aminotab_productive_topten[f'topten_moderate_mutated_shared_{iso2}_VH_AA'] = stats_iso2_Aminotab_productive_topten[column_seqid_iso2_aa].isin(stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_seqid)
            stats_iso2_Aminotab_productive_topten[f'topten_high_mutated_shared_{iso2}_VH_AA'] = stats_iso2_Aminotab_productive_topten[column_seqid_iso2_aa].isin(stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_seqid)
            
            stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_tab = stats_iso2_Aminotab_productive_topten[stats_iso2_Aminotab_productive_topten[f'topten_low_mutated_shared_{iso2}_VH_AA']]
            stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_tab = stats_iso2_Aminotab_productive_topten[stats_iso2_Aminotab_productive_topten[f'topten_moderate_mutated_shared_{iso2}_VH_AA']]
            stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_tab = stats_iso2_Aminotab_productive_topten[stats_iso2_Aminotab_productive_topten[f'topten_high_mutated_shared_{iso2}_VH_AA']]
            
            stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_AA = stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_tab[column_vdj_region_iso2].tolist() if column_vdj_region_iso2 in stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_tab.columns else []
            stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_AA = stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_tab[column_vdj_region_iso2].tolist() if column_vdj_region_iso2 in stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_tab.columns else []
            stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_AA = stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_tab[column_vdj_region_iso2].tolist() if column_vdj_region_iso2 in stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_tab.columns else []
            
            print(f"\n top ten Low Mutated Shared {iso2} Sequences:")
            for i, sequence in enumerate(stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print(f"\n top ten Moderate Mutated Shared {iso2} Sequences:")
            for i, sequence in enumerate(stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print(f"\n top ten High Mutated Shared {iso2} Sequences:")
            for i, sequence in enumerate(stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")

            # ============================================
            # FASTA ALIGNMENT AND VH SEQUENCE COMPLETION
            # ============================================
            
            from pathlib import Path
            from Bio.Align import PairwiseAligner
            
            def load_imgt_fasta(fasta_file):
                sequences = {}
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequences[record.id] = str(record.seq)
                return sequences
            
            def find_best_alignment(fasta_sequences, query_sequence):
                aligner = PairwiseAligner()
                aligner.mode = "global"
                best_score = -float('inf')
                best_reference = None
                best_seq_id = None
                for seq_id, reference_sequence in fasta_sequences.items():
                    score = aligner.score(reference_sequence, query_sequence)
                    if score > best_score:
                        best_score = score
                        best_reference = reference_sequence
                        best_seq_id = seq_id
                return best_seq_id, best_reference
            
            def calculate_similarity(seq1, seq2):
                matches = sum(a == b for a, b in zip(seq1, seq2))
                return matches / len(seq1) if len(seq1) > 0 else 0
            
            def extract_and_concatenate(query_sequence, reference_sequence):
                aligner = PairwiseAligner()
                aligner.mode = "global"
                alignment = aligner.align(reference_sequence, query_sequence)[0]
                aligned_reference = alignment.aligned[1]
                if len(aligned_reference) == 0:
                    return query_sequence
                start_index_reference = aligned_reference[-1][1]
                upstream_region = reference_sequence[:start_index_reference]
                overlap_length = 0
                for i in range(1, len(upstream_region) + 1):
                    overlap_candidate = upstream_region[-i:]
                    query_start_candidate = query_sequence[:i]
                    if len(overlap_candidate) == len(query_start_candidate):
                        similarity = calculate_similarity(overlap_candidate, query_start_candidate)
                        if similarity > 0.6:
                            overlap_length = i
                            break
                trimmed_upstream = upstream_region[:-overlap_length] if overlap_length > 0 else upstream_region
                return trimmed_upstream + query_sequence
            
            def keep_one_repeat(sequence, min_length=5):
                substrings = {}
                result = list(sequence)
                for i in range(len(sequence)):
                    for j in range(i + min_length, len(sequence) + 1):
                        substring = sequence[i:j]
                        if substring in substrings:
                            start_index = sequence.find(substring, substrings[substring] + 1)
                            if start_index != -1:
                                for k in range(start_index, start_index + len(substring)):
                                    result[k] = ""
                        else:
                            substrings[substring] = i
                return "".join(result)
            
            def get_vdj_region_column_name(df):
                if 'V.D.J.REGION' in df.columns:
                    return 'V.D.J.REGION'
                elif 'V-D-J-REGION' in df.columns:
                    return 'V-D-J-REGION'
                else:
                    raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' is present in the DataFrame.")
            
            def process_and_save_sequences(input_df, fasta_sequences, output_filename, iso1, iso2):
                if input_df.empty:
                    print(f"Skipping {output_filename} - empty DataFrame")
                    return
                try:
                    vdj_region_col = get_vdj_region_column_name(input_df)
                except KeyError:
                    vdj_region_col = 'V-D-J-REGION'
                    if vdj_region_col not in input_df.columns:
                        print(f"Skipping {output_filename} - no VDJ region column found")
                        return
                completed_sequences = []
                for _, row in input_df.iterrows():
                    query_sequence = row[vdj_region_col]
                    _, best_reference = find_best_alignment(fasta_sequences, query_sequence)
                    completed_sequence = extract_and_concatenate(query_sequence, best_reference) if best_reference else query_sequence
                    completed_sequences.append(completed_sequence)
                input_df = input_df.copy()
                input_df["Completed_VH_Sequences"] = completed_sequences
                input_df["cleaned_sequence"] = input_df["Completed_VH_Sequences"].apply(keep_one_repeat)
                # Use the specified output directory (MODIFIED)
                vh_sequences_folder = os.path.join(self.output_dir, "VH_Sequences")
                os.makedirs(vh_sequences_folder, exist_ok=True)
                output_csv = os.path.join(vh_sequences_folder, output_filename)
                input_df.to_csv(output_csv, index=False)
                print(f"Output saved to {output_csv}")
            
            # Process VH sequences
            fasta_sequences = load_imgt_fasta(fasta_file_path)
            
            process_and_save_sequences(divergent_low_mutated_shared_iso1_VH_AA_tab.copy(), fasta_sequences, f"output_divergent_low_mutated_shared_{iso1}_VH_AA_tab.csv", iso1, iso2)
            process_and_save_sequences(divergent_moderate_mutated_shared_iso1_VH_AA_tab.copy(), fasta_sequences, f"output_divergent_moderate_mutated_shared_{iso1}_VH_AA_tab.csv", iso1, iso2)
            process_and_save_sequences(divergent_high_mutated_shared_iso1_VH_AA_tab.copy(), fasta_sequences, f"output_divergent_high_mutated_shared_{iso1}_VH_AA_tab.csv", iso1, iso2)
            process_and_save_sequences(non_divergent_low_mutated_shared_iso1_VH_AA_tab.copy(), fasta_sequences, f"output_non_divergent_low_mutated_shared_{iso1}_VH_AA_tab.csv", iso1, iso2)
            process_and_save_sequences(non_divergent_moderate_mutated_shared_iso1_VH_AA_tab.copy(), fasta_sequences, f"output_non_divergent_moderate_mutated_shared_{iso1}_VH_AA_tab.csv", iso1, iso2)
            process_and_save_sequences(non_divergent_high_mutated_shared_iso1_VH_AA_tab.copy(), fasta_sequences, f"output_non_divergent_high_mutated_shared_{iso1}_VH_AA_tab.csv", iso1, iso2)
            process_and_save_sequences(stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten_tab.copy(), fasta_sequences, f"output_stats_tab_low_mutated_shared_{iso1}_VDJ_sorted_topten.csv", iso1, iso2)
            process_and_save_sequences(stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten_tab.copy(), fasta_sequences, f"output_stats_tab_moderate_mutated_shared_{iso1}_VDJ_sorted_topten.csv", iso1, iso2)
            process_and_save_sequences(stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten_tab.copy(), fasta_sequences, f"output_stats_tab_high_mutated_shared_{iso1}_VDJ_sorted_topten.csv", iso1, iso2)
            
            # Process iso2 VH sequences
            process_and_save_sequences(clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AA.copy(), fasta_sequences, f"output_{iso2}_clonally_related_to_low_mutated_divergent_shared_{iso1}_seqid_in_AA.csv", iso1, iso2)
            process_and_save_sequences(clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AA.copy(), fasta_sequences, f"output_{iso2}_clonally_related_to_moderate_mutated_divergent_shared_{iso1}_seqid_in_AA.csv", iso1, iso2)
            process_and_save_sequences(clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AA.copy(), fasta_sequences, f"output_{iso2}_clonally_related_to_high_mutated_divergent_shared_{iso1}_seqid_in_AA.csv", iso1, iso2)
            process_and_save_sequences(stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten_tab.copy(), fasta_sequences, f"output_stats_tab_low_mutated_shared_{iso2}_VDJ_sorted_topten.csv", iso1, iso2)
            process_and_save_sequences(stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten_tab.copy(), fasta_sequences, f"output_stats_tab_moderate_mutated_shared_{iso2}_VDJ_sorted_topten.csv", iso1, iso2)
            process_and_save_sequences(stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten_tab.copy(), fasta_sequences, f"output_stats_tab_high_mutated_shared_{iso2}_VDJ_sorted_topten.csv", iso1, iso2)
            self.logger.log_memory("After VH sequence extraction")

            # Amino acid composition analysis
            import re
            from collections import Counter

            def extract_and_process_aa(sequence):
                if not isinstance(sequence, str):
                    return '', ''
                segments = [seg.strip() for seg in sequence.split('|') if seg.strip()]

                def process_mutation(mut_str):
                    explicit_change = re.search(r'([A-Z])\d+>([A-Z])', mut_str)
                    if explicit_change:
                        return explicit_change.group(1), explicit_change.group(2)
                    repeated_aa = re.search(r'([A-Z])\d+;\s*([A-Z])\d+', mut_str)
                    if repeated_aa and repeated_aa.group(1) == repeated_aa.group(2):
                        return repeated_aa.group(1), repeated_aa.group(1)
                    return None, None

                germline = []
                mutated = []
                for segment in segments:
                    for mut_notation in [x.strip() for x in segment.split(',') if x.strip()]:
                        germ_aa, mut_aa = process_mutation(mut_notation)
                        if germ_aa is not None and mut_aa is not None:
                            germline.append(germ_aa)
                            mutated.append(mut_aa)
                return ''.join(germline), ''.join(mutated)

            # Define EMBOSS amino acid categories
            emboss_categories = {
                "Tiny": set("ACGST"),
                "Small": set("ABCDGNPSTV"),
                "Aliphatic": set("AILV"),
                "Aromatic": set("FHWY"),
                "Non-polar": set("ACFGILMPVWY"),
                "Polar": set("DEHKNQRSTZ"),
                "Charged": set("BDEHKRZ"),
                "Basic": set("KRH"),
                "Acidic": set("BDEZ")
            }

            def calculate_emboss_composition(sequence):
                if not sequence:
                    return {category: 0 for category in emboss_categories.keys()}
                total_aa = len(sequence)
                standard_aas = set("ACDEFGHIKLMNPQRSTVWYBZ")
                aa_count = Counter()
                category_counts = {category: 0 for category in emboss_categories.keys()}

                for aa in sequence:
                    aa_count[aa] += 1
                    if aa in standard_aas:
                        for category, aa_set in emboss_categories.items():
                            if aa in aa_set:
                                category_counts[category] += 1

                category_percentages = {
                    category: (count / total_aa) * 100 if total_aa > 0 else 0
                    for category, count in category_counts.items()
                }
                return category_percentages

            # Analyze amino acid composition for shared isotype 1
            mut_shared_iso1 = stats_tab_mutated_shared_iso1_VDJ.copy()
            shared_iso1_VDJ_extract = set(mut_shared_iso1['seqid'])

            iso1_AA_change_table1 = stats_iso1_Aminotab_change.copy()

            if 'Sequence.ID' in iso1_AA_change_table1.columns:
                col_seqid_aa = 'Sequence.ID'
            elif 'Sequence ID' in iso1_AA_change_table1.columns:
                col_seqid_aa = 'Sequence ID'
            else:
                col_seqid_aa = 'Sequence ID'

            iso1_AA_change_table1[f'Match_in_mut_shared_{iso1}_seqid'] = iso1_AA_change_table1[col_seqid_aa].isin(shared_iso1_VDJ_extract)
            iso1_AA_change_table1[f'Match_in_mut_shared_{iso1}_seqid'] = iso1_AA_change_table1[f'Match_in_mut_shared_{iso1}_seqid'].astype(str).str.strip()

            iso1_AA_change_table2 = iso1_AA_change_table1.loc[iso1_AA_change_table1[f'Match_in_mut_shared_{iso1}_seqid'] == 'True'].copy()

            # Get CDR column names for isotype 1
            if 'CDR3.IMGT' in iso1_AA_change_table2.columns:
                col_cdr3_aa = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in iso1_AA_change_table2.columns:
                col_cdr3_aa = 'CDR3-IMGT'
            else:
                col_cdr3_aa = 'CDR3-IMGT'

            if 'CDR2.IMGT' in iso1_AA_change_table2.columns:
                col_cdr2_aa = 'CDR2.IMGT'
            elif 'CDR2-IMGT' in iso1_AA_change_table2.columns:
                col_cdr2_aa = 'CDR2-IMGT'
            else:
                col_cdr2_aa = 'CDR2-IMGT'

            germline_column_cdr3 = []
            mutated_column_cdr3 = []
            germline_column_cdr2 = []
            mutated_column_cdr2 = []

            for cdr3, cdr2 in zip(iso1_AA_change_table2[col_cdr3_aa], iso1_AA_change_table2[col_cdr2_aa]):
                germline_cdr3, mutated_cdr3 = extract_and_process_aa(cdr3)
                germline_column_cdr3.append(germline_cdr3)
                mutated_column_cdr3.append(mutated_cdr3)
                germline_cdr2, mutated_cdr2 = extract_and_process_aa(cdr2)
                germline_column_cdr2.append(germline_cdr2)
                mutated_column_cdr2.append(mutated_cdr2)

            iso1_AA_change_table2.loc[:, 'Germline_CDR3'] = germline_column_cdr3
            iso1_AA_change_table2.loc[:, 'Mutated_CDR3'] = mutated_column_cdr3
            iso1_AA_change_table2.loc[:, 'Germline_CDR2'] = germline_column_cdr2
            iso1_AA_change_table2.loc[:, 'Mutated_CDR2'] = mutated_column_cdr2

            Germline_CDR3_concat = ''.join(iso1_AA_change_table2['Germline_CDR3'])
            Mutated_CDR3_concat = ''.join(iso1_AA_change_table2['Mutated_CDR3'])
            Germline_CDR2_concat = ''.join(iso1_AA_change_table2['Germline_CDR2'])
            Mutated_CDR2_concat = ''.join(iso1_AA_change_table2['Mutated_CDR2'])

            Germline_CDR3_composition = calculate_emboss_composition(Germline_CDR3_concat)
            Mutated_CDR3_composition = calculate_emboss_composition(Mutated_CDR3_concat)
            Germline_CDR2_composition = calculate_emboss_composition(Germline_CDR2_concat)
            Mutated_CDR2_composition = calculate_emboss_composition(Mutated_CDR2_concat)

            composition_results_iso1 = {
                'Germline_CDR3': Germline_CDR3_composition,
                'Mutated_CDR3': Mutated_CDR3_composition,
                'Germline_CDR2': Germline_CDR2_composition,
                'Mutated_CDR2': Mutated_CDR2_composition,
            }

            # Use the specified output directory
            data_folder = os.path.join(self.output_dir, "Data")
            os.makedirs(data_folder, exist_ok=True)
            output_file_path = os.path.join(data_folder, f"Shared_{iso1}_AminoAcid_Composition.xlsx")

            composition_df_iso1 = pd.DataFrame(composition_results_iso1)
            composition_df_iso1.to_excel(output_file_path, index=True)
            print(f"File saved successfully in {output_file_path}")

            # Analyze amino acid composition for shared isotype 2
            mut_shared_iso2 = stats_tab_mutated_shared_iso2_VDJ.copy()
            shared_iso2_VDJ_extract = set(mut_shared_iso2['seqid'])

            iso2_AA_change_table1 = stats_iso2_Aminotab_change.copy()

            if 'Sequence.ID' in iso2_AA_change_table1.columns:
                col_seqid_aa_iso2 = 'Sequence.ID'
            elif 'Sequence ID' in iso2_AA_change_table1.columns:
                col_seqid_aa_iso2 = 'Sequence ID'
            else:
                col_seqid_aa_iso2 = 'Sequence ID'

            iso2_AA_change_table1[f'Match_in_mut_shared_{iso2}_seqid'] = iso2_AA_change_table1[col_seqid_aa_iso2].isin(shared_iso2_VDJ_extract)
            iso2_AA_change_table1[f'Match_in_mut_shared_{iso2}_seqid'] = iso2_AA_change_table1[f'Match_in_mut_shared_{iso2}_seqid'].astype(str).str.strip()

            iso2_AA_change_table2 = iso2_AA_change_table1.loc[iso2_AA_change_table1[f'Match_in_mut_shared_{iso2}_seqid'] == 'True'].copy()

            if 'CDR3.IMGT' in iso2_AA_change_table2.columns:
                col_cdr3_aa_iso2 = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in iso2_AA_change_table2.columns:
                col_cdr3_aa_iso2 = 'CDR3-IMGT'
            else:
                col_cdr3_aa_iso2 = 'CDR3-IMGT'

            if 'CDR2.IMGT' in iso2_AA_change_table2.columns:
                col_cdr2_aa_iso2 = 'CDR2.IMGT'
            elif 'CDR2-IMGT' in iso2_AA_change_table2.columns:
                col_cdr2_aa_iso2 = 'CDR2-IMGT'
            else:
                col_cdr2_aa_iso2 = 'CDR2-IMGT'

            germline_column_cdr3_iso2 = []
            mutated_column_cdr3_iso2 = []
            germline_column_cdr2_iso2 = []
            mutated_column_cdr2_iso2 = []

            for cdr3, cdr2 in zip(iso2_AA_change_table2[col_cdr3_aa_iso2], iso2_AA_change_table2[col_cdr2_aa_iso2]):
                germline_cdr3, mutated_cdr3 = extract_and_process_aa(cdr3)
                germline_column_cdr3_iso2.append(germline_cdr3)
                mutated_column_cdr3_iso2.append(mutated_cdr3)
                germline_cdr2, mutated_cdr2 = extract_and_process_aa(cdr2)
                germline_column_cdr2_iso2.append(germline_cdr2)
                mutated_column_cdr2_iso2.append(mutated_cdr2)

            iso2_AA_change_table2.loc[:, 'Germline_CDR3'] = germline_column_cdr3_iso2
            iso2_AA_change_table2.loc[:, 'Mutated_CDR3'] = mutated_column_cdr3_iso2
            iso2_AA_change_table2.loc[:, 'Germline_CDR2'] = germline_column_cdr2_iso2
            iso2_AA_change_table2.loc[:, 'Mutated_CDR2'] = mutated_column_cdr2_iso2

            Germline_CDR3_concat_iso2 = ''.join(iso2_AA_change_table2['Germline_CDR3'])
            Mutated_CDR3_concat_iso2 = ''.join(iso2_AA_change_table2['Mutated_CDR3'])
            Germline_CDR2_concat_iso2 = ''.join(iso2_AA_change_table2['Germline_CDR2'])
            Mutated_CDR2_concat_iso2 = ''.join(iso2_AA_change_table2['Mutated_CDR2'])

            Germline_CDR3_composition_iso2 = calculate_emboss_composition(Germline_CDR3_concat_iso2)
            Mutated_CDR3_composition_iso2 = calculate_emboss_composition(Mutated_CDR3_concat_iso2)
            Germline_CDR2_composition_iso2 = calculate_emboss_composition(Germline_CDR2_concat_iso2)
            Mutated_CDR2_composition_iso2 = calculate_emboss_composition(Mutated_CDR2_concat_iso2)

            composition_results_iso2 = {
                'Germline_CDR3': Germline_CDR3_composition_iso2,
                'Mutated_CDR3': Mutated_CDR3_composition_iso2,
                'Germline_CDR2': Germline_CDR2_composition_iso2,
                'Mutated_CDR2': Mutated_CDR2_composition_iso2,
            }

            output_file_path_iso2 = os.path.join(data_folder, f"Shared_{iso2}_AminoAcid_Composition.xlsx")
            composition_df_iso2 = pd.DataFrame(composition_results_iso2)
            composition_df_iso2.to_excel(output_file_path_iso2, index=True)
            print(f"File saved successfully in {output_file_path_iso2}")
            self.logger.log_memory("After amino acid composition analysis")

           # Generate summary output
            class PrintCapture:
                def __init__(self):
                    self.output = StringIO()
                    self.saved_stdout = sys.stdout
                    sys.stdout = self.output

                def get_output(self):
                    return self.output.getvalue()

                def reset(self):
                    sys.stdout = self.saved_stdout

            capture = PrintCapture()

            # ============================================
            # COMPLETE SUMMARY STATISTICS - ALL SECTIONS
            # ============================================
            
            print("=" * 80)
            print(f"BCR ANALYSIS COMPLETE SUMMARY: {iso1} vs {iso2}")
            print("=" * 80)
            
            # --- SECTION 1: CLONE COUNTS ---
            print("\n" + "-" * 40)
            print("1. CLONE COUNTS")
            print("-" * 40)
            print(f'Number of {iso1} clones: {unique_strings_iso1_productive_VDJ}')
            print(f'Number of {iso2} clones: {unique_strings_iso2_productive_VDJ}')
            
            # --- SECTION 2: SHARED/UNIQUE CLONES ---
            print("\n" + "-" * 40)
            print("2. SHARED AND UNIQUE CLONES")
            print("-" * 40)
            print(f'Number of shared {iso1} clones: {unique_strings_shared_iso1}')
            print(f'Percentage of shared {iso1} clones: {percentage_shared_iso1_clones:.2f}%')
            print(f'Number of unique {iso1} clones: {unique_strings_unique_iso1}')
            print(f'Percentage of unique {iso1} clones: {percentage_unique_iso1_clones:.2f}%')
            print(f'Number of shared {iso2} clones: {unique_strings_shared_iso2}')
            print(f'Percentage of shared {iso2} clones: {percentage_shared_iso2_clones:.2f}%')
            print(f'Number of unique {iso2} clones: {unique_strings_unique_iso2}')
            print(f'Percentage of unique {iso2} clones: {percentage_unique_iso2_clones:.2f}%')
            
            # --- SECTION 3: CDR MUTATION PERCENTAGES ---
            print("\n" + "-" * 40)
            print("3. CDR MUTATION ANALYSIS")
            print("-" * 40)
            print(f'Mutated {iso1}-biased avg % CDR2 non-silent mutation: {mean_CDR2_stats_tab_mutated_iso1_bias:.4f}')
            print(f'Mutated {iso1}-biased avg % CDR3 non-silent mutation: {mean_CDR3_stats_tab_mutated_iso1_bias:.4f}')
            print(f'Mutated {iso2}-biased avg % CDR2 non-silent mutation: {mean_CDR2_stats_tab_mutated_iso2_bias:.4f}')
            print(f'Mutated {iso2}-biased avg % CDR3 non-silent mutation: {mean_CDR3_stats_tab_mutated_iso2_bias:.4f}')
            print(f'Mutated unique {iso1} avg % CDR2 non-silent mutation: {mean_CDR2_stats_tab_mutated_Unique_iso1:.4f}')
            print(f'Mutated unique {iso1} avg % CDR3 non-silent mutation: {mean_CDR3_stats_tab_mutated_Unique_iso1:.4f}')
            print(f'Unmutated {iso1}-biased avg % CDR2 non-silent mutation: {mean_CDR2_stats_tab_unmutated_iso1_bias:.4f}')
            print(f'Unmutated {iso1}-biased avg % CDR3 non-silent mutation: {mean_CDR3_stats_tab_unmutated_iso1_bias:.4f}')
            print(f'Unmutated {iso2}-biased avg % CDR2 non-silent mutation: {mean_CDR2_stats_tab_unmutated_iso2_bias:.4f}')
            print(f'Unmutated {iso2}-biased avg % CDR3 non-silent mutation: {mean_CDR3_stats_tab_unmutated_iso2_bias:.4f}')
            print(f'Unmutated unique {iso1} avg % CDR2 non-silent mutation: {mean_CDR2_stats_tab_unmutated_Unique_iso1:.4f}')
            print(f'Unmutated unique {iso1} avg % CDR3 non-silent mutation: {mean_CDR3_stats_tab_unmutated_Unique_iso1:.4f}')
            
            # --- SECTION 4: MUTATED CLONE COUNTS ---
            print("\n" + "-" * 40)
            print("4. MUTATED CLONE COUNTS")
            print("-" * 40)
            print(f'Number of mutated {iso1}-biased clones: {unique_strings_mutated_iso1_bias_VDJ}')
            print(f'Percentage of mutated {iso1}-biased clones: {percentage_mutated_iso1_bias_clones:.2f}%')
            print(f'Mean copy number mutated {iso1}-biased: {mean_copynumber_mutated_iso1_bias_VDJ:.2f}')
            print(f'Number of mutated {iso2}-biased clones: {unique_strings_mutated_iso2_bias_VDJ}')
            print(f'Percentage of mutated {iso2}-biased clones: {percentage_mutated_iso2_bias_clones:.2f}%')
            print(f'Mean copy number mutated {iso2}-biased: {mean_copynumber_mutated_iso2_bias_VDJ:.2f}')
            print(f'Number of mutated unique {iso1} clones: {unique_strings_mutated_unique_iso1_VDJ}')
            print(f'Percentage of mutated unique {iso1} clones: {percentage_mutated_unique_iso1_clones:.2f}%')
            print(f'Mean copy number mutated unique {iso1}: {mean_copynumber_mutated_unique_iso1_VDJ:.2f}')
            print(f'Number of mutated Total {iso1} clones: {unique_strings_mutated_Total_iso1_VDJ}')
            print(f'Percentage of mutated Total {iso1} clones: {percentage_mutated_Total_iso1_clones:.2f}%')
            print(f'Mean copy number mutated Total {iso1}: {mean_copynumber_mutated_Total_iso1_VDJ:.2f}')
            print(f'Number of mutated Total {iso2} clones: {unique_strings_mutated_Total_iso2_VDJ}')
            print(f'Percentage of mutated Total {iso2} clones: {percentage_mutated_Total_iso2_clones:.2f}%')
            print(f'Mean copy number mutated Total {iso2}: {mean_copynumber_mutated_Total_iso2_VDJ:.2f}')
            
            # --- SECTION 5: UNMUTATED CLONE COUNTS ---
            print("\n" + "-" * 40)
            print("5. UNMUTATED CLONE COUNTS")
            print("-" * 40)
            print(f'Number of unmutated {iso1}-biased clones: {unique_strings_unmutated_iso1_bias_VDJ}')
            print(f'Percentage of unmutated {iso1}-biased clones: {percentage_unmutated_iso1_bias_clones:.2f}%')
            print(f'Mean copy number unmutated {iso1}-biased: {mean_copynumber_unmutated_iso1_bias_VDJ:.2f}')
            print(f'Number of unmutated {iso2}-biased clones: {unique_strings_unmutated_iso2_bias_VDJ}')
            print(f'Percentage of unmutated {iso2}-biased clones: {percentage_unmutated_iso2_bias_clones:.2f}%')
            print(f'Mean copy number unmutated {iso2}-biased: {mean_copynumber_unmutated_iso2_bias_VDJ:.2f}')
            print(f'Number of unmutated unique {iso1} clones: {unique_strings_unmutated_unique_iso1_VDJ}')
            print(f'Percentage of unmutated unique {iso1} clones: {percentage_unmutated_unique_iso1_clones:.2f}%')
            print(f'Mean copy number unmutated unique {iso1}: {mean_copynumber_unmutated_unique_iso1_VDJ:.2f}')
            print(f'Number of unmutated Total {iso1} clones: {unique_strings_unmutated_Total_iso1_VDJ}')
            print(f'Percentage of unmutated Total {iso1} clones: {percentage_unmutated_Total_iso1_clones:.2f}%')
            print(f'Mean copy number unmutated Total {iso1}: {mean_copynumber_unmutated_Total_iso1_VDJ:.2f}')
            print(f'Number of unmutated Total {iso2} clones: {unique_strings_unmutated_Total_iso2_VDJ}')
            print(f'Percentage of unmutated Total {iso2} clones: {percentage_unmutated_Total_iso2_clones:.2f}%')
            print(f'Mean copy number unmutated Total {iso2}: {mean_copynumber_unmutated_Total_iso2_VDJ:.2f}')
            
            # --- SECTION 6: MUTATION LEVEL ANALYSIS ---
            print("\n" + "-" * 40)
            print("6. MUTATION LEVEL ANALYSIS (Low/Moderate/High)")
            print("-" * 40)
            print(f'Mean copy number low mutated {iso1}: {mean1:.2f}')
            print(f'Mean copy number low mutated {iso2}: {mean2:.2f}')
            print(f'Mean copy number moderate mutated {iso1}: {mean3:.2f}')
            print(f'Mean copy number moderate mutated {iso2}: {mean4:.2f}')
            print(f'Mean copy number high mutated {iso1}: {mean5:.2f}')
            print(f'Mean copy number high mutated {iso2}: {mean6:.2f}')
            
            # --- SECTION 7: CDR3 NET CHARGE ---
            print("\n" + "-" * 40)
            print("7. CDR3 NET CHARGE ANALYSIS (Top 10 Expanded)")
            print("-" * 40)
            print(f'Mean CDR3 net charge low mutated {iso1}: {mean_stats_tab_low_mutated_shared_iso1_VDJ_sorted_topten:.4f}')
            print(f'Mean CDR3 net charge moderate mutated {iso1}: {mean_stats_tab_moderate_mutated_shared_iso1_VDJ_sorted_topten:.4f}')
            print(f'Mean CDR3 net charge high mutated {iso1}: {mean_stats_tab_high_mutated_shared_iso1_VDJ_sorted_topten:.4f}')
            print(f'Mean CDR3 net charge low mutated {iso2}: {mean_stats_tab_low_mutated_shared_iso2_VDJ_sorted_topten:.4f}')
            print(f'Mean CDR3 net charge moderate mutated {iso2}: {mean_stats_tab_moderate_mutated_shared_iso2_VDJ_sorted_topten:.4f}')
            print(f'Mean CDR3 net charge high mutated {iso2}: {mean_stats_tab_high_mutated_shared_iso2_VDJ_sorted_topten:.4f}')
            print(f'Low mutated {iso2}/{iso1} net charge ratio: {ratio_stats_tab_low_mutated_shared_iso2_iso1_VDJ_sorted_topten:.4f}')
            print(f'Moderate mutated {iso2}/{iso1} net charge ratio: {ratio_stats_tab_moderate_mutated_shared_iso2_iso1_VDJ_sorted_topten:.4f}')
            print(f'High mutated {iso2}/{iso1} net charge ratio: {ratio_stats_tab_high_mutated_shared_iso2_iso1_VDJ_sorted_topten:.4f}')
            
            # --- SECTION 8: DIVERGENT CDR3 ANALYSIS ---
            print("\n" + "-" * 40)
            print(f"8. DIVERGENT CDR3 ANALYSIS (>= {cdr3_thresh} CDR3 mutations)")
            print("-" * 40)
            print(f'Number of divergent CDR3 low mutated {iso1} clones: {unique_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Number of divergent CDR3 moderate mutated {iso1} clones: {unique_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Number of divergent CDR3 high mutated {iso1} clones: {unique_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Number of non-divergent CDR3 low mutated {iso1} clones: {unique_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Number of non-divergent CDR3 moderate mutated {iso1} clones: {unique_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Number of non-divergent CDR3 high mutated {iso1} clones: {unique_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            
            # --- SECTION 9: DIVERGENT COPY NUMBERS ---
            print("\n" + "-" * 40)
            print("9. DIVERGENT CLONE COPY NUMBERS")
            print("-" * 40)
            print(f'Mean copy number divergent CDR3 low mutated {iso1}: {mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone:.2f}')
            print(f'Mean copy number divergent CDR3 moderate mutated {iso1}: {mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone:.2f}')
            print(f'Mean copy number divergent CDR3 high mutated {iso1}: {mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone:.2f}')
            print(f'Mean copy number non-divergent CDR3 low mutated {iso1}: {mean_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone:.2f}')
            print(f'Mean copy number non-divergent CDR3 moderate mutated {iso1}: {mean_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone:.2f}')
            print(f'Mean copy number non-divergent CDR3 high mutated {iso1}: {mean_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone:.2f}')
            print(f'Sum copy number divergent CDR3 low mutated {iso1}: {sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Sum copy number divergent CDR3 moderate mutated {iso1}: {sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Sum copy number divergent CDR3 high mutated {iso1}: {sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Sum copy number non-divergent CDR3 low mutated {iso1}: {sum_non_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Sum copy number non-divergent CDR3 moderate mutated {iso1}: {sum_non_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            print(f'Sum copy number non-divergent CDR3 high mutated {iso1}: {sum_non_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3mut_clone}')
            
            # --- SECTION 10: CLONALLY RELATED ISO2 ---
            print("\n" + "-" * 40)
            print(f"10. CLONALLY RELATED {iso2} TO DIVERGENT {iso1}")
            print("-" * 40)
            print(f'Mean copy number {iso2} clonally related to low mutated {iso1}: {mean_clonally_related_low_mutated_shared_iso2:.2f}')
            print(f'Mean copy number {iso2} clonally related to moderate mutated {iso1}: {mean_clonally_related_moderate_mutated_shared_iso2:.2f}')
            print(f'Mean copy number {iso2} clonally related to high mutated {iso1}: {mean_clonally_related_high_mutated_shared_iso2:.2f}')
            print(f'Sum copy number {iso2} clonally related to low mutated {iso1}: {sum_clonally_related_low_mutated_shared_iso2}')
            print(f'Sum copy number {iso2} clonally related to moderate mutated {iso1}: {sum_clonally_related_moderate_mutated_shared_iso2}')
            print(f'Sum copy number {iso2} clonally related to high mutated {iso1}: {sum_clonally_related_high_mutated_shared_iso2}')
            
            # --- SECTION 11: CDR3CDR2 DIVERGENT ANALYSIS ---
            print("\n" + "-" * 40)
            print(f"11. CDR3+CDR2 DIVERGENT ANALYSIS (>= {cdr3cdr2_thresh} mutations)")
            print("-" * 40)
            print(f'Number of divergent CDR3CDR2 low mutated {iso1} clones: {unique_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            print(f'Number of divergent CDR3CDR2 moderate mutated {iso1} clones: {unique_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            print(f'Number of divergent CDR3CDR2 high mutated {iso1} clones: {unique_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            print(f'Mean copy number divergent CDR3CDR2 low mutated {iso1}: {mean_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone:.2f}')
            print(f'Mean copy number divergent CDR3CDR2 moderate mutated {iso1}: {mean_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone:.2f}')
            print(f'Mean copy number divergent CDR3CDR2 high mutated {iso1}: {mean_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone:.2f}')
            print(f'Sum copy number divergent CDR3CDR2 low mutated {iso1}: {sum_divergent_tab_low_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            print(f'Sum copy number divergent CDR3CDR2 moderate mutated {iso1}: {sum_divergent_tab_moderate_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            print(f'Sum copy number divergent CDR3CDR2 high mutated {iso1}: {sum_divergent_tab_high_mutated_shared_iso1_VDJ_CDR3CDR2mut_clone}')
            
            # --- SECTION 12: CDR3CDR2 CLONALLY RELATED ---
            print("\n" + "-" * 40)
            print(f"12. CDR3+CDR2 CLONALLY RELATED {iso2}")
            print("-" * 40)
            print(f'Mean copy number {iso2} clonally related to low mutated {iso1} (CDR3CDR2): {mean_clonally_related_low_mutated_shared_CDR3CDR2_iso2:.2f}')
            print(f'Mean copy number {iso2} clonally related to moderate mutated {iso1} (CDR3CDR2): {mean_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2:.2f}')
            print(f'Mean copy number {iso2} clonally related to high mutated {iso1} (CDR3CDR2): {mean_clonally_related_high_mutated_shared_CDR3CDR2_iso2:.2f}')
            print(f'Sum copy number {iso2} clonally related to low mutated {iso1} (CDR3CDR2): {sum_clonally_related_low_mutated_shared_CDR3CDR2_iso2}')
            print(f'Sum copy number {iso2} clonally related to moderate mutated {iso1} (CDR3CDR2): {sum_clonally_related_moderate_mutated_shared_CDR3CDR2_iso2}')
            print(f'Sum copy number {iso2} clonally related to high mutated {iso1} (CDR3CDR2): {sum_clonally_related_high_mutated_shared_CDR3CDR2_iso2}')
            
            # --- SECTION 13: SEQUENCE COUNTS ---
            print("\n" + "-" * 40)
            print("13. VH SEQUENCE EXTRACTION COUNTS")
            print("-" * 40)
            print(f'Divergent low mutated {iso1} VH sequences: {len(divergent_low_mutated_VH_AAsequences)}')
            print(f'Divergent moderate mutated {iso1} VH sequences: {len(divergent_moderate_mutated_VH_AAsequences)}')
            print(f'Divergent high mutated {iso1} VH sequences: {len(divergent_high_mutated_VH_AAsequences)}')
            print(f'Non-divergent low mutated {iso1} VH sequences: {len(non_divergent_low_mutated_VH_AAsequences)}')
            print(f'Non-divergent moderate mutated {iso1} VH sequences: {len(non_divergent_moderate_mutated_VH_AAsequences)}')
            print(f'Non-divergent high mutated {iso1} VH sequences: {len(non_divergent_high_mutated_VH_AAsequences)}')
            print(f'Clonally related low mutated {iso2} VH sequences: {len(clonally_related_low_mutated_shared_iso2_seqid_in_iso2_AAsequences)}')
            print(f'Clonally related moderate mutated {iso2} VH sequences: {len(clonally_related_moderate_mutated_shared_iso2_seqid_in_iso2_AAsequences)}')
            print(f'Clonally related high mutated {iso2} VH sequences: {len(clonally_related_high_mutated_shared_iso2_seqid_in_iso2_AAsequences)}')
            
            print("\n" + "=" * 80)
            print("END OF SUMMARY")
            print("=" * 80)

            captured_output = capture.get_output()
            capture.reset()

            output_lines = captured_output.splitlines()
            df_output = pd.DataFrame(output_lines, columns=["Analysis Summary"])
            # Save summary to output directory
            summary_folder = os.path.join(self.output_dir, "Summary")
            os.makedirs(summary_folder, exist_ok=True)
            output_excel_file = os.path.join(summary_folder, f"Analysis_Summary_{iso1}_{iso2}.xlsx")
            df_output.to_excel(output_excel_file, index=False)
            print(f"Captured output has been saved to {output_excel_file}")

            # Final logging
            self.logger.log_memory("Final memory state")
            self.logger.log_info("Analysis completed successfully!")
            
            # Save log file to output directory
            logs_folder = os.path.join(self.output_dir, "Logs")
            os.makedirs(logs_folder, exist_ok=True)
            log_file_path = os.path.join(logs_folder, f"BCR_Analysis_Log_{iso1}_{iso2}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
            with open(log_file_path, 'w') as f:
                f.write(self.logger.get_full_log())
            self.logger.log_info(f"Log saved to: {log_file_path}")
            
            # Emit log content
            if self.enable_logging:
                self.log_updated.emit(self.logger.get_full_log())
            
            self.analysis_completed.emit("Analysis completed successfully!")

        except Exception as e:
            import traceback
            error_msg = f"{str(e)}\n{traceback.format_exc()}"
            self.logger.log_error(f"Analysis failed: {str(e)}")
            self.logger.log_error(traceback.format_exc())
            self.error_occurred.emit(error_msg)


class Screen2(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Plots Display")
        self.setGeometry(200, 200, 800, 600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.current_index = 0
        self.figures = []
        self.canvas = None

        self.load_all_plots()

        next_button = QPushButton("Next Plot")
        next_button.clicked.connect(self.show_next_plot)
        self.layout.addWidget(next_button)

    def load_all_plots(self):
        for fignum in plt.get_fignums():
            fig = plt.figure(fignum)
            self.figures.append(fig)
        if self.figures:
            self.display_plot(0)

    def display_plot(self, index):
        if self.canvas:
            self.layout.removeWidget(self.canvas)
            self.canvas.deleteLater()
        self.canvas = FigureCanvas(self.figures[index])
        self.layout.insertWidget(0, self.canvas)

    def show_next_plot(self):
        if not self.figures:
            return
        self.current_index = (self.current_index + 1) % len(self.figures)
        self.display_plot(self.current_index)


class LogViewerDialog(QDialog):
    """Dialog window to display analysis logs."""
    def __init__(self, log_content="", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Analysis Log Viewer")
        self.setGeometry(150, 150, 1000, 700)
        self.log_content = log_content
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Header
        header = QLabel("Analysis Log Viewer")
        header.setFont(QFont("Courier", 16, QFont.Bold))
        header.setStyleSheet("color: #cccccc;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # Filter buttons
        filter_layout = QHBoxLayout()
        
        self.show_all_btn = QPushButton("Show All")
        self.show_all_btn.clicked.connect(lambda: self.filter_logs("ALL"))
        
        self.show_errors_btn = QPushButton("Errors Only")
        self.show_errors_btn.clicked.connect(lambda: self.filter_logs("ERROR"))
        
        self.show_warnings_btn = QPushButton("Warnings")
        self.show_warnings_btn.clicked.connect(lambda: self.filter_logs("WARNING"))
        
        self.show_exclusions_btn = QPushButton("Exclusions")
        self.show_exclusions_btn.clicked.connect(lambda: self.filter_logs("EXCLUSION"))
        
        self.show_validation_btn = QPushButton("Validation")
        self.show_validation_btn.clicked.connect(lambda: self.filter_logs("VALIDATION"))
        
        self.show_memory_btn = QPushButton("Memory")
        self.show_memory_btn.clicked.connect(lambda: self.filter_logs("MEMORY"))
        
        for btn in [self.show_all_btn, self.show_errors_btn, self.show_warnings_btn, 
                    self.show_exclusions_btn, self.show_validation_btn, self.show_memory_btn]:
            btn.setStyleSheet("""
                QPushButton {
                    background-color: #444444;
                    color: white;
                    padding: 8px 16px;
                    border-radius: 5px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #666666;
                }
            """)
            filter_layout.addWidget(btn)
        
        layout.addLayout(filter_layout)
        
        # Log text area
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Courier", 10))
        self.log_text.setStyleSheet("""
            QTextEdit {
                background-color: #1a1a1a;
                color: #00ff00;
                border: 1px solid #444444;
                border-radius: 5px;
                padding: 10px;
            }
        """)
        self.log_text.setText(self.log_content)
        layout.addWidget(self.log_text)
        
        # Bottom buttons
        bottom_layout = QHBoxLayout()
        
        save_btn = QPushButton("Save Log")
        save_btn.clicked.connect(self.save_log)
        save_btn.setStyleSheet("""
            QPushButton {
                background-color: #555555;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #777777;
            }
        """)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #333333;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #555555;
            }
        """)
        
        bottom_layout.addWidget(save_btn)
        bottom_layout.addWidget(close_btn)
        layout.addLayout(bottom_layout)
        
        self.setLayout(layout)
        
        # Apply dark theme
        self.setStyleSheet("""
            QDialog {
                background-color: #2a2a2a;
            }
        """)
    
    def filter_logs(self, filter_type):
        if filter_type == "ALL":
            self.log_text.setText(self.log_content)
        else:
            filtered_lines = [line for line in self.log_content.split('\n') if filter_type in line]
            self.log_text.setText('\n'.join(filtered_lines))
    
    def save_log(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Log File", "", "Text Files (*.txt);;All Files (*)")
        if file_path:
            with open(file_path, 'w') as f:
                f.write(self.log_text.toPlainText())
            QMessageBox.information(self, "Saved", f"Log saved to {file_path}")
    
    def update_log(self, content):
        self.log_content = content
        self.log_text.setText(content)

# =============================================================================
# FASTBCR PLOT VIEWER
# =============================================================================

class FastBCRPlotViewer(QDialog):
    """Dialog to view generated plots with navigation."""
    
    def __init__(self, output_dir, parent=None):
        super().__init__(parent)
        self.output_dir = output_dir
        self.current_index = 0
        self.plot_files = []
        self.load_plot_files()
        self.init_ui()
        
    def load_plot_files(self):
        """Load available plot files, preferring PNG over PDF."""
        if not os.path.exists(self.output_dir):
            return
            
        all_files = os.listdir(self.output_dir)
        
        # Get PNG files (preferred for display)
        png_files = sorted([f for f in all_files if f.lower().endswith('.png')])
        
        # Get PDF files that don't have PNG versions
        # EXCLUDE Rplots.pdf and similar R default plot files
        pdf_files = sorted([f for f in all_files 
                           if f.lower().endswith('.pdf') 
                           and not f.lower().startswith('rplots')
                           and not f.lower().startswith('r.plots')])
        
        pdf_without_png = [f for f in pdf_files 
                          if f.replace('.pdf', '.png').replace('.PDF', '.png') not in png_files]
        
        # Combine: PNG files first, then PDFs without PNG versions
        self.plot_files = png_files + pdf_without_png
        
    def init_ui(self):
        self.setWindowTitle("FastBCR Plot Viewer")
        self.setMinimumSize(900, 700)
        
        layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("FastBCR Analysis Plots")
        title.setStyleSheet("font-size: 18px; font-weight: bold; color: #E74C3C;")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # Plot display area
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setStyleSheet("QScrollArea { background-color: #2b2b2b; border: none; }")
        
        self.plot_label = QLabel()
        self.plot_label.setAlignment(Qt.AlignCenter)
        self.plot_label.setStyleSheet("QLabel { background-color: #2b2b2b; color: white; padding: 20px; }")
        self.plot_label.setMinimumSize(800, 500)
        
        self.scroll_area.setWidget(self.plot_label)
        layout.addWidget(self.scroll_area)
        
        # File name label
        self.filename_label = QLabel()
        self.filename_label.setAlignment(Qt.AlignCenter)
        self.filename_label.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.filename_label)
        
        # Navigation
        nav_layout = QHBoxLayout()
        
        self.prev_button = QPushButton("◄ Previous")
        self.prev_button.clicked.connect(self.show_previous)
        self.prev_button.setStyleSheet("""
            QPushButton {
                background-color: #555;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #666; }
            QPushButton:disabled { background-color: #333; color: #666; }
        """)
        nav_layout.addWidget(self.prev_button)
        
        self.page_label = QLabel()
        self.page_label.setAlignment(Qt.AlignCenter)
        self.page_label.setStyleSheet("color: white; font-size: 14px;")
        nav_layout.addWidget(self.page_label)
        
        self.next_button = QPushButton("Next ►")
        self.next_button.clicked.connect(self.show_next)
        self.next_button.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #2980b9; }
            QPushButton:disabled { background-color: #333; color: #666; }
        """)
        nav_layout.addWidget(self.next_button)
        
        layout.addLayout(nav_layout)
        
        # Bottom buttons
        button_layout = QHBoxLayout()
        
        open_folder_btn = QPushButton("📁 Open Output Folder")
        open_folder_btn.clicked.connect(self.open_folder)
        open_folder_btn.setStyleSheet("""
            QPushButton {
                background-color: #f39c12;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #e67e22; }
        """)
        button_layout.addWidget(open_folder_btn)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #555;
                color: white;
                border: none;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #666; }
        """)
        button_layout.addWidget(close_btn)
        
        layout.addLayout(button_layout)
        
        # Show first plot
        self.show_current_plot()
        
    def show_current_plot(self):
        """Display the current plot."""
        if not self.plot_files:
            self.plot_label.setText("No plot files found in output directory.")
            self.page_label.setText("0 / 0")
            self.prev_button.setEnabled(False)
            self.next_button.setEnabled(False)
            return
            
        file_name = self.plot_files[self.current_index]
        file_path = os.path.join(self.output_dir, file_name)
        
        # Update navigation
        self.page_label.setText(f"{self.current_index + 1} / {len(self.plot_files)}")
        self.prev_button.setEnabled(self.current_index > 0)
        self.next_button.setEnabled(self.current_index < len(self.plot_files) - 1)
        self.filename_label.setText(file_name)
        
        if file_path.endswith('.png'):
            self.display_image(file_path)
        elif file_path.endswith('.pdf'):
            # Try to find PNG version
            png_path = file_path.replace('.pdf', '.png')
            if os.path.exists(png_path):
                self.display_image(png_path)
            else:
                self.plot_label.setText(
                    f"PDF File: {file_name}\n\n"
                    "PDF files cannot be displayed directly.\n"
                    "Click 'Open Output Folder' to view.\n\n"
                    "Tip: Re-run analysis to generate PNG versions."
                )
        else:
            self.plot_label.setText(f"Unsupported file type: {file_name}")
            
    def display_image(self, file_path):
        """Display an image file."""
        pixmap = QPixmap(file_path)
        if not pixmap.isNull():
            # Scale to fit the label while maintaining aspect ratio
            label_size = self.plot_label.size()
            scaled_pixmap = pixmap.scaled(
                label_size.width() - 40,
                label_size.height() - 40,
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )
            self.plot_label.setPixmap(scaled_pixmap)
        else:
            self.plot_label.setText(f"Error loading image: {os.path.basename(file_path)}")
            
    def show_previous(self):
        if self.current_index > 0:
            self.current_index -= 1
            self.show_current_plot()
            
    def show_next(self):
        if self.current_index < len(self.plot_files) - 1:
            self.current_index += 1
            self.show_current_plot()
            
    def open_folder(self):
        """Open the output folder in file explorer."""
        import subprocess
        import platform
        
        if platform.system() == "Windows":
            os.startfile(self.output_dir)
        elif platform.system() == "Darwin":  # macOS
            subprocess.run(["open", self.output_dir])
        else:  # Linux
            subprocess.run(["xdg-open", self.output_dir])


# =============================================================================
# FASTBCR LOG VIEWER
# =============================================================================

class FastBCRLogViewer(QDialog):
    """Dialog window to display FastBCR analysis logs."""
    
    def __init__(self, log_content="", parent=None):
        super().__init__(parent)
        self.setWindowTitle("FastBCR Analysis Log Viewer")
        self.setGeometry(150, 150, 1000, 700)
        self.log_content = log_content
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Header
        header = QLabel("FastBCR Analysis Log Viewer")
        header.setFont(QFont("Arial", 16, QFont.Bold))
        header.setStyleSheet("color: #4CAF50;")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # Filter buttons
        filter_layout = QHBoxLayout()
        
        filter_buttons = [
            ("Show All", "ALL"),
            ("Clusters", "cluster"),
            ("Exclusions", "exclud"),
            ("Errors", "Error"),
            ("Warnings", "Warning"),
            ("Memory", "MEMORY"),
            ("Summary", "SUMMARY"),
        ]
        
        for btn_text, filter_key in filter_buttons:
            btn = QPushButton(btn_text)
            btn.clicked.connect(lambda checked, f=filter_key: self.filter_logs(f))
            btn.setStyleSheet("""
                QPushButton {
                    background-color: #444444;
                    color: white;
                    padding: 6px 12px;
                    border-radius: 4px;
                    font-size: 11px;
                }
                QPushButton:hover {
                    background-color: #666666;
                }
            """)
            filter_layout.addWidget(btn)
        
        layout.addLayout(filter_layout)
        
        # Statistics summary panel
        self.stats_label = QLabel()
        self.stats_label.setStyleSheet("""
            QLabel {
                background-color: #1a3a1a;
                color: #00ff00;
                padding: 10px;
                border: 1px solid #004400;
                border-radius: 5px;
                font-family: Courier;
                font-size: 11px;
            }
        """)
        self.update_stats_summary()
        layout.addWidget(self.stats_label)
        
        # Log text area
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Courier", 10))
        self.log_text.setStyleSheet("""
            QTextEdit {
                background-color: #1a1a1a;
                color: #00ff00;
                border: 1px solid #444444;
                border-radius: 5px;
                padding: 10px;
            }
        """)
        self.log_text.setText(self.log_content)
        layout.addWidget(self.log_text)
        
        # Bottom buttons
        bottom_layout = QHBoxLayout()
        
        save_btn = QPushButton("💾 Save Log")
        save_btn.clicked.connect(self.save_log)
        save_btn.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #1976D2; }
        """)
        
        copy_btn = QPushButton("📋 Copy to Clipboard")
        copy_btn.clicked.connect(self.copy_to_clipboard)
        copy_btn.setStyleSheet("""
            QPushButton {
                background-color: #FF9800;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #F57C00; }
        """)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #666666;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover { background-color: #888888; }
        """)
        
        bottom_layout.addWidget(save_btn)
        bottom_layout.addWidget(copy_btn)
        bottom_layout.addStretch()
        bottom_layout.addWidget(close_btn)
        layout.addLayout(bottom_layout)
        
        self.setLayout(layout)
        
        # Apply dark theme
        self.setStyleSheet("QDialog { background-color: #2a2a2a; }")
    
    def update_stats_summary(self):
        """Parse log content and display key statistics."""
        stats = {
            'iso1_clusters': self.extract_value(r'(\d+)\s*clusters found', 0),
            'iso2_clusters': self.extract_value(r'(\d+)\s*clusters found', 1),
            'shared_pairs': self.extract_value(r'Shared cluster pairs found:\s*(\d+)'),
            'singletons': self.count_occurrences('singleton'),
            'non_productive': self.count_occurrences('productive'),
            'errors': self.count_occurrences('Error'),
            'warnings': self.count_occurrences('Warning'),
        }
        
        # Extract runtime if available
        import re
        runtime_match = re.search(r'End time:.*?(\d{2}:\d{2}:\d{2})', self.log_content)
        start_match = re.search(r'Start time:.*?(\d{2}:\d{2}:\d{2})', self.log_content)
        
        summary_text = f"""
┌─────────────────────────────────────────────────────────────────────────────┐
│  FASTBCR ANALYSIS STATISTICS SUMMARY                                        │
├─────────────────────────────────────────────────────────────────────────────┤
│  Isotype 1 Clusters: {str(stats['iso1_clusters']).ljust(10)}  │  Isotype 2 Clusters: {str(stats['iso2_clusters']).ljust(10)}  │
│  Shared Cluster Pairs: {str(stats['shared_pairs']).ljust(8)}  │  Singletons Mentioned: {str(stats['singletons']).ljust(8)}  │
│  Non-productive Refs: {str(stats['non_productive']).ljust(9)}  │  Errors: {str(stats['errors']).ljust(5)}  Warnings: {str(stats['warnings']).ljust(5)}  │
└─────────────────────────────────────────────────────────────────────────────┘
        """.strip()
        
        self.stats_label.setText(summary_text)
    
    def extract_value(self, pattern, occurrence=0):
        """Extract a numeric value from log using regex pattern."""
        import re
        matches = re.findall(pattern, self.log_content, re.IGNORECASE)
        if matches and len(matches) > occurrence:
            return matches[occurrence]
        return "N/A"
    
    def count_occurrences(self, keyword):
        """Count occurrences of a keyword in the log."""
        return self.log_content.lower().count(keyword.lower())
    
    def filter_logs(self, filter_type):
        """Filter log content by keyword."""
        if filter_type == "ALL":
            self.log_text.setText(self.log_content)
        else:
            filtered_lines = [
                line for line in self.log_content.split('\n') 
                if filter_type.lower() in line.lower()
            ]
            if filtered_lines:
                self.log_text.setText('\n'.join(filtered_lines))
            else:
                self.log_text.setText(f"No lines found containing '{filter_type}'")
    
    def save_log(self):
        """Save log to file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Log File", 
            f"FastBCR_Log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
            "Text Files (*.txt);;All Files (*)"
        )
        if file_path:
            with open(file_path, 'w') as f:
                f.write(self.log_text.toPlainText())
            QMessageBox.information(self, "Saved", f"Log saved to:\n{file_path}")
    
    def copy_to_clipboard(self):
        """Copy log content to clipboard."""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.log_text.toPlainText())
        QMessageBox.information(self, "Copied", "Log content copied to clipboard!")




# ===========================================================================
# COMPIGS TAB (Placeholder - Import from existing code)
# =============================================================================

class CompIgSTab(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("CompIgS - Comparative Immunoglobulin Sequence Analysis")
        self.setGeometry(100, 100, 1100, 700)

        main_layout = QVBoxLayout()
        form_layout = QVBoxLayout()

        # Apply palette
        palette = QPalette()
        gradient = QLinearGradient(0, 0, 1, 1)
        gradient.setColorAt(0.0, QColor(10, 10, 10))
        gradient.setColorAt(1.0, QColor(40, 40, 40))
        palette.setBrush(QPalette.Window, QBrush(gradient))
        palette.setColor(QPalette.WindowText, QColor(255, 255, 255))
        palette.setColor(QPalette.Base, QColor(20, 20, 20))
        palette.setColor(QPalette.Text, QColor(200, 200, 200))
        self.setPalette(palette)

        # Header
        header_label = QLabel("CompIgS")
        header_label.setAlignment(Qt.AlignCenter)
        header_label.setFont(QFont("Orbitron", 22, QFont.Bold))
        header_label.setStyleSheet("QLabel { color: #cccccc; }")
        main_layout.addWidget(header_label)

        # Isotype Name Input Section
        isotype_group = QGroupBox("Isotype Names")
        isotype_group.setStyleSheet("""
            QGroupBox {
                color: #cccccc;
                border: 1px solid #666666;
                border-radius: 5px;
                margin-top: 10px;
                padding: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        isotype_layout = QHBoxLayout()

        iso1_label = QLabel("Isotype 1 Name (e.g., IgE, IgA, IgM):")
        iso1_label.setFont(QFont("Courier", 12))
        iso1_label.setStyleSheet("color: #999999;")
        self.iso1_input = QLineEdit()
        self.iso1_input.setText("IgE")
        self.iso1_input.setStyleSheet("""
            QLineEdit {
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px;
                border-radius: 5px;
            }
        """)

        iso2_label = QLabel("Isotype 2 Name (e.g., IgG1, IgG2, IgD):")
        iso2_label.setFont(QFont("Courier", 12))
        iso2_label.setStyleSheet("color: #999999;")
        self.iso2_input = QLineEdit()
        self.iso2_input.setText("IgG1")
        self.iso2_input.setStyleSheet("""
            QLineEdit {
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px;
                border-radius: 5px;
            }
        """)

        isotype_layout.addWidget(iso1_label)
        isotype_layout.addWidget(self.iso1_input)
        isotype_layout.addWidget(iso2_label)
        isotype_layout.addWidget(self.iso2_input)
        isotype_group.setLayout(isotype_layout)
        main_layout.addWidget(isotype_group)

        # Output Directory Section (INSERT THIS ENTIRE BLOCK)
        output_group = QGroupBox("Output Directory")
        output_group.setStyleSheet("""
            QGroupBox {
                color: #cccccc;
                border: 1px solid #666666;
                border-radius: 5px;
                margin-top: 10px;
                padding: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        output_layout = QHBoxLayout()
        
        output_label = QLabel("Save Results To:")
        output_label.setFont(QFont("Courier", 12))
        output_label.setStyleSheet("color: #999999;")
        
        self.output_dir_input = QLineEdit()
        # Set default to Downloads/CompIgS_Analysis
        default_output = os.path.join(os.path.expanduser("~"), "Downloads", "CompIgS_Analysis")
        self.output_dir_input.setText(default_output)
        self.output_dir_input.setStyleSheet("""
            QLineEdit {
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px;
                border-radius: 5px;
            }
            QLineEdit:hover {
                border: 1px solid #999999;
            }
        """)
        
        output_browse_button = QPushButton("Browse")
        output_browse_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(100, 100, 100, 0.2);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px 12px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: rgba(150, 150, 150, 0.4);
            }
        """)
        output_browse_button.clicked.connect(self.browse_output_directory)
        
        open_folder_button = QPushButton("Open Folder")
        open_folder_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(100, 100, 100, 0.2);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px 12px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: rgba(150, 150, 150, 0.4);
            }
        """)
        open_folder_button.clicked.connect(self.open_output_directory)
        
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_dir_input)
        output_layout.addWidget(output_browse_button)
        output_layout.addWidget(open_folder_button)
        output_group.setLayout(output_layout)
        main_layout.addWidget(output_group)

        # Divergent Clone Threshold Selection Section
        threshold_group = QGroupBox("Divergent Clone Mutation Thresholds")
        threshold_group.setStyleSheet("""
            QGroupBox {
                color: #cccccc;
                border: 1px solid #666666;
                border-radius: 5px;
                margin-top: 10px;
                padding: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        threshold_layout = QHBoxLayout()

        # CDR3 threshold selector
        cdr3_label = QLabel("CDR3 Non-silent Mutations (≥):")
        cdr3_label.setFont(QFont("Courier", 11))
        cdr3_label.setStyleSheet("color: #999999;")
        
        from PyQt5.QtWidgets import QComboBox
        self.cdr3_threshold_combo = QComboBox()
        self.cdr3_threshold_combo.addItems(["1", "2", "3"])
        self.cdr3_threshold_combo.setCurrentIndex(1)  # Default to 2
        self.cdr3_threshold_combo.setStyleSheet("""
            QComboBox {
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px;
                border-radius: 5px;
                min-width: 60px;
            }
            QComboBox:hover {
                border: 1px solid #999999;
            }
            QComboBox::drop-down {
                border: none;
            }
            QComboBox::down-arrow {
                image: none;
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid #cccccc;
                margin-right: 5px;
            }
            QComboBox QAbstractItemView {
                background-color: rgba(30, 30, 30, 0.95);
                color: #cccccc;
                selection-background-color: #666666;
            }
        """)

        # CDR3CDR2 threshold selector
        cdr3cdr2_label = QLabel("CDR3+CDR2 Non-silent Mutations (≥):")
        cdr3cdr2_label.setFont(QFont("Courier", 11))
        cdr3cdr2_label.setStyleSheet("color: #999999;")
        
        self.cdr3cdr2_threshold_combo = QComboBox()
        self.cdr3cdr2_threshold_combo.addItems(["1", "2", "3"])
        self.cdr3cdr2_threshold_combo.setCurrentIndex(1)  # Default to 2
        self.cdr3cdr2_threshold_combo.setStyleSheet("""
            QComboBox {
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px;
                border-radius: 5px;
                min-width: 60px;
            }
            QComboBox:hover {
                border: 1px solid #999999;
            }
            QComboBox::drop-down {
                border: none;
            }
            QComboBox::down-arrow {
                image: none;
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid #cccccc;
                margin-right: 5px;
            }
            QComboBox QAbstractItemView {
                background-color: rgba(30, 30, 30, 0.95);
                color: #cccccc;
                selection-background-color: #666666;
            }
        """)

        threshold_layout.addWidget(cdr3_label)
        threshold_layout.addWidget(self.cdr3_threshold_combo)
        threshold_layout.addStretch()
        threshold_layout.addWidget(cdr3cdr2_label)
        threshold_layout.addWidget(self.cdr3cdr2_threshold_combo)
        threshold_layout.addStretch()
        
        threshold_group.setLayout(threshold_layout)
        main_layout.addWidget(threshold_group)

        # Divergent Clone Information Note (ADD THIS SECTION)
        divergent_info_frame = QFrame()
        divergent_info_frame.setStyleSheet("""
            QFrame {
                background-color: rgba(255, 152, 0, 0.15);
                border: 1px solid #FF9800;
                border-radius: 5px;
                padding: 5px;
            }
        """)
        divergent_info_layout = QHBoxLayout()
        divergent_info_layout.setContentsMargins(10, 8, 10, 8)
        
        warning_icon = QLabel("⚠️")
        warning_icon.setFont(QFont("Arial", 14))
        
        divergent_note = QLabel(
            "<b>Note:</b> Divergent clones are identified for <b style='color: #FF9800;'>Isotype 1</b> only. "
            "These are sequences from Isotype 1 that differ significantly from all Isotype 2 sequences."
        )
        divergent_note.setWordWrap(True)
        divergent_note.setStyleSheet("color: #cccccc; font-size: 11px;")
        
        divergent_info_layout.addWidget(warning_icon)
        divergent_info_layout.addWidget(divergent_note, 1)
        divergent_info_frame.setLayout(divergent_info_layout)
        main_layout.addWidget(divergent_info_frame)

        # Logging Options
        logging_group = QGroupBox("Logging Options")
        logging_group.setStyleSheet("""
            QGroupBox {
                color: #cccccc;
                border: 1px solid #666666;
                border-radius: 5px;
                margin-top: 10px;
                padding: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        logging_layout = QHBoxLayout()
        
        self.enable_logging_checkbox = QCheckBox("Enable Detailed Logging")
        self.enable_logging_checkbox.setChecked(True)
        self.enable_logging_checkbox.setFont(QFont("Courier", 11))
        self.enable_logging_checkbox.setEnabled(True)  # Ensure it's enabled
        self.enable_logging_checkbox.setStyleSheet("""
            QCheckBox {
                color: #999999;
                spacing: 8px;
            }
            QCheckBox::indicator {
                width: 18px;
                height: 18px;
            }
            QCheckBox::indicator:unchecked {
                background-color: #333333;
                border: 2px solid #666666;
                border-radius: 3px;
            }
            QCheckBox::indicator:checked {
                background-color: #00aa00;
                border: 2px solid #00cc00;
                border-radius: 3px;
            }
            QCheckBox::indicator:checked:hover {
                background-color: #00cc00;
            }
            QCheckBox::indicator:unchecked:hover {
                border: 2px solid #999999;
            }
        """)
        
        logging_info = QLabel("(Tracks: exclusions, errors, validation, runtime, memory usage)")
        logging_info.setFont(QFont("Courier", 9))
        logging_info.setStyleSheet("color: #777777;")
        
        self.view_log_button = QPushButton("View Log")
        self.view_log_button.setEnabled(False)
        self.view_log_button.setStyleSheet("""
            QPushButton {
                background-color: rgba(100, 100, 100, 0.2);
                color: #cccccc;
                border: 1px solid #666666;
                padding: 6px 12px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: rgba(150, 150, 150, 0.4);
            }
            QPushButton:disabled {
                color: #555555;
                border: 1px solid #444444;
            }
        """)
        self.view_log_button.clicked.connect(self.show_log_viewer)
        
        logging_layout.addWidget(self.enable_logging_checkbox)
        logging_layout.addWidget(logging_info)
        logging_layout.addStretch()
        logging_layout.addWidget(self.view_log_button)
        logging_group.setLayout(logging_layout)
        main_layout.addWidget(logging_group)
        
        # Store log content
        self.current_log_content = ""


        

        # File input fields - labels will update based on isotype names
        self.file_inputs = []
        self.file_labels = []
        
        file_label_templates = [
            "(Ig1){iso1} IMGT StatClonotype Output File:",
            "(Ig2){iso2} IMGT StatClonotype Output File:",
            "(Ig1){iso1} 8-V-region-nt-mutation-statistics table:",
            "(Ig2){iso2} 8-V-region-nt-mutation-statistics table:",
            "(Ig1){iso1} 5-AA-sequences table:",
            "(Ig2){iso2} 5-AA-sequences table:",
            "(Ig1){iso1} 7-V-REGION-mutation-and-AA-change-table:",
            "(Ig2){iso2} 7-V-REGION-mutation-and-AA-change-table:",
            "IMGT mouse VH reference:"
        ]

        for i, label_template in enumerate(file_label_templates):
            row_layout = QHBoxLayout()
            label_text = label_template.format(iso1=self.iso1_input.text(), iso2=self.iso2_input.text())
            label = QLabel(label_text)
            label.setFont(QFont("Courier", 11))
            label.setStyleSheet("color: #999999;")
            label.setMinimumWidth(400)
            self.file_labels.append(label)
            
            line_edit = QLineEdit()
            line_edit.setStyleSheet("""
                QLineEdit {
                    background-color: rgba(20, 20, 20, 0.9);
                    color: #cccccc;
                    border: 1px solid #666666;
                    padding: 6px;
                    border-radius: 5px;
                }
                QLineEdit:hover {
                    border: 1px solid #999999;
                }
            """)
            browse_button = QPushButton("Browse")
            browse_button.setStyleSheet("""
                QPushButton {
                    background-color: rgba(100, 100, 100, 0.2);
                    color: #cccccc;
                    border: 1px solid #666666;
                    padding: 6px 12px;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: rgba(150, 150, 150, 0.4);
                }
            """)
            browse_button.clicked.connect(lambda _, le=line_edit: self.browse_file(le))

            row_layout.addWidget(label)
            row_layout.addWidget(line_edit)
            row_layout.addWidget(browse_button)

            form_layout.addLayout(row_layout)
            self.file_inputs.append(line_edit)

        # Connect isotype name changes to update labels
        self.iso1_input.textChanged.connect(self.update_labels)
        self.iso2_input.textChanged.connect(self.update_labels)

        main_layout.addLayout(form_layout)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("background-color: #666666; height: 2px;")
        main_layout.addWidget(separator)

        # Progress bar
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #666666;
                border-radius: 6px;
                text-align: center;
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
            }
            QProgressBar::chunk {
                background-color: #666666;
                width: 20px;
            }
        """)

        self.status_label = QLabel("Status: Waiting for input")
        self.status_label.setAlignment(Qt.AlignCenter)
        self.status_label.setFont(QFont("Courier", 12))
        self.status_label.setStyleSheet("color: #999999;")

        main_layout.addWidget(self.progress_bar)
        main_layout.addWidget(self.status_label)

        # Buttons
        button_layout = QHBoxLayout()
        
        start_button = QPushButton("Start Analysis")
        start_button.setStyleSheet("""
            QPushButton {
                background-color: #666666;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #999999;
            }
        """)
        start_button.clicked.connect(self.start_analysis)

        self.goto_screen2_button = QPushButton("View Plots", self)
        self.goto_screen2_button.setStyleSheet("""
            QPushButton {
                background-color: #444444;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #777777;
            }
        """)
        self.goto_screen2_button.clicked.connect(self.gotoScreen2)

        clear_all_button = QPushButton("Clear All")
        clear_all_button.setStyleSheet("""
            QPushButton {
                background-color: #555555;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #888888;
            }
        """)
        clear_all_button.clicked.connect(self.clear_all_inputs)

        exit_button = QPushButton("Exit")
        exit_button.setStyleSheet("""
            QPushButton {
                background-color: #333333;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #666666;
            }
        """)
        exit_button.clicked.connect(self.close)

        button_layout.addWidget(start_button)
        button_layout.addWidget(self.goto_screen2_button)
        button_layout.addWidget(clear_all_button)
        button_layout.addWidget(exit_button)
        main_layout.addLayout(button_layout)

        self.setLayout(main_layout)

    def update_labels(self):
        """Update file input labels when isotype names change."""
        iso1 = self.iso1_input.text() or "Iso1"
        iso2 = self.iso2_input.text() or "Iso2"
        
        label_templates = [
            f"(Ig1){iso1} IMGT StatClonotype Output File:",
            f"(Ig2){iso2} IMGT StatClonotype Output File:",
            f"(Ig1){iso1} 8-V-region-nt-mutation-statistics table:",
            f"(Ig2){iso2} 8-V-region-nt-mutation-statistics table:",
            f"(Ig1){iso1} 5-AA-sequences table:",
            f"(Ig2){iso2} 5-AA-sequences table:",
            f"(Ig1){iso1} 7-V-REGION-mutation-and-AA-change-table:",
            f"(Ig2){iso2} 7-V-REGION-mutation-and-AA-change-table:",
            "IMGT mouse VH reference:"
        ]
        
        for i, label in enumerate(self.file_labels):
            label.setText(label_templates[i])

    def browse_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if file_path:
            line_edit.setText(file_path)

    def start_analysis(self):
        file_paths = [le.text() for le in self.file_inputs]
        iso1_name = self.iso1_input.text().strip()
        iso2_name = self.iso2_input.text().strip()
        
        # Get threshold values from combo boxes
        cdr3_threshold = int(self.cdr3_threshold_combo.currentText())
        cdr3cdr2_threshold = int(self.cdr3cdr2_threshold_combo.currentText())
        
        # Get logging preference
        enable_logging = self.enable_logging_checkbox.isChecked()
        
        # Get output directory (ADD THIS SECTION)
        output_dir = self.output_dir_input.text().strip()
        if not output_dir:
            output_dir = os.path.join(os.path.expanduser("~"), "Downloads", "CompIgS_Analysis")
        
        # Create output directory if it doesn't exist
        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            QMessageBox.warning(self, "Directory Error", f"Could not create output directory: {e}")
            return

        if not iso1_name or not iso2_name:
            QMessageBox.warning(self, "Input Error", "Please provide names for both isotypes.")
            return

        if all(file_paths):
            self.thread = AnalysisThread(file_paths, iso1_name, iso2_name, cdr3_threshold, cdr3cdr2_threshold, enable_logging, output_dir)
            self.thread.progress_updated.connect(self.update_progress)
            self.thread.analysis_completed.connect(self.on_analysis_complete)
            self.thread.error_occurred.connect(self.on_error)
            self.thread.log_updated.connect(self.update_log_content)

            self.thread.start()
            log_status = " [Logging enabled]" if enable_logging else ""
            self.status_label.setText(f"Status: Analyzing {iso1_name} vs {iso2_name} (CDR3≥{cdr3_threshold}, CDR3CDR2≥{cdr3cdr2_threshold}){log_status}...")
        else:
            QMessageBox.warning(self, "Input Error", "Please provide all required files.")
            self.status_label.setText("Error: Missing required files.")

    def update_progress(self, value):
        self.progress_bar.setValue(value)

    def on_analysis_complete(self, message):
        output_dir = self.output_dir_input.text()
        self.status_label.setText(f"{message} - Results saved to: {output_dir}")
        QMessageBox.information(self, "Success", f"{message}\n\nResults saved to:\n{output_dir}")


## Summary of Output Directory Structure

    def on_error(self, error_message):
        self.status_label.setText(f"Error: {error_message[:50]}...")
        QMessageBox.critical(self, "Error", error_message)

    def show_log_viewer(self):
        """Open the log viewer dialog."""
        if self.current_log_content:
            dialog = LogViewerDialog(self.current_log_content, self)
            dialog.exec_()
        else:
            QMessageBox.information(self, "No Log", "No log data available. Run an analysis first.")
    
    def update_log_content(self, log_content):
        """Update the stored log content and enable view button."""
        self.current_log_content = log_content
        self.view_log_button.setEnabled(True)

    def gotoScreen2(self):
        self.screen2 = Screen2()
        self.screen2.exec_()

    def clear_all_inputs(self):
        for line_edit in self.file_inputs:
            line_edit.clear()
        plt.close('all')
        self.status_label.setText("Status: All input fields cleared.")
    
    def browse_output_directory(self):
        """Browse for output directory."""
        dir_path = QFileDialog.getExistingDirectory(
            self, 
            "Select Output Directory",
            self.output_dir_input.text()
        )
        if dir_path:
            self.output_dir_input.setText(dir_path)
    
    def open_output_directory(self):
        """Open the output directory in file explorer."""
        output_dir = self.output_dir_input.text()
        
        if not output_dir:
            QMessageBox.warning(self, "No Directory", "Please specify an output directory first.")
            return
        
        # Create directory if it doesn't exist
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Could not create directory: {e}")
                return
        
        # Open in file explorer
        import subprocess
        import platform
        
        try:
            if platform.system() == "Windows":
                os.startfile(output_dir)
            elif platform.system() == "Darwin":  # macOS
                subprocess.run(["open", output_dir])
            else:  # Linux
                subprocess.run(["xdg-open", output_dir])
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not open directory: {e}")

# =============================================================================
# METHOD COMPARISON TAB
# =============================================================================

class ComparisonTab(QWidget):
    """Tab showing comparison between the two analysis methods"""
    
    def __init__(self):
        super().__init__()
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Title
        title = QLabel("Method Comparison")
        title.setFont(QFont("Arial", 18, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("color: #9C27B0;")
        layout.addWidget(title)
        
        # Comparison content
        content = QTextEdit()
        content.setReadOnly(True)
        content.setStyleSheet("""
            QTextEdit {
                background-color: #1a1a1a;
                color: #cccccc;
                border: 1px solid #444444;
                border-radius: 5px;
                padding: 15px;
            }
        """)
        content.setHtml("""
<style>
table { border-collapse: collapse; width: 100%; margin: 10px 0; }
th, td { border: 1px solid #555; padding: 10px; text-align: left; }
th { background-color: #333; color: #4CAF50; }
tr:nth-child(even) { background-color: #252525; }
h2 { color: #4CAF50; border-bottom: 2px solid #4CAF50; padding-bottom: 5px; }
h3 { color: #2196F3; }
.highlight { background-color: #3a3a00; padding: 2px 5px; }
.vdj { color: #2196F3; font-weight: bold; }
.vj { color: #4CAF50; font-weight: bold; }
</style>

<h3>Use <span class="vdj">CompIgS (VDJ Matching)</span> when:</h3>
<ul>
<li>You have IMGT HighV-QUEST output files</li>
<li>You need <b>clonotype/cluster definitions</b> based on exact VDJ gene usage</li>
</ul>

<h3>Use <span class="vj">fastBCR (VJ + k-mer)</span> when:</h3>
<ul>
<li>You have AIRR-formatted data </li>
<li>You want <b>within isotype clustering</b> using VJ genes and k-mer seeding</li>
<li>You need to identify <b>shared clusters between isotypes</b> using junction similarity</li>
</ul>

<h2>📊 Feature Comparison</h2>

<table>
<tr>
<th>Feature</th>
<th class="vdj">CompIgS (VDJ)</th>
<th class="vj">fastBCR-based (VJ + k-mer)</th>
</tr>
<tr>
<td><b>Within-Isotype Clustering</b></td>
<td>Exact V+D+J gene match</td>
<td>VJ genes + k-mer seeding + consensus</td>
</tr>
<tr>
<td><b>Between-Isotype Sharing</b></td>
<td>Exact V+D+J match</td>
<td>Shared VJ + % junction similarity</td>
</tr>
<tr>
<td><b>Input Format</b></td>
<td>IMGT StatClonotype + mutation tables</td>
<td>AIRR TSV (db-pass.tsv)</td>
</tr>
<tr>
<td><b>D Gene Requirement</b></td>
<td>Required</td>
<td>Not required</td>
</tr>
<tr>
<td><b>Indel Handling</b></td>
<td>Low (strict matching)</td>
<td>Higher (k-mer and similarity-based)</td>
</tr>
<tr>
<td><b>Mutation Analysis</b></td>
<td>Detailed CDR2/CDR3 statistics</td>
<td>SHM % calculation</td>
</tr>
<tr>
<td><b>Divergent Clone ID</b></td>
<td>CDR3 AA differences</td>
<td>Junction AA differences between isotypes</td>
</tr>
<tr>
<td><b>Visualizations</b></td>
<td>Bar plots, pie charts, distributions</td>
<td>+ Lineage trees, MSA, heatmaps</td>
</tr>
<tr>
<td><b>Language</b></td>
<td>Python</td>
<td>R (fastBCR package)</td>
</tr>
</table>

<h2>🧬Method</h2>

<h3><span class="vdj">CompIgS</span></h3>
<p><b>Within-Isotype Clustering:</b> Same V+D+J genes</p>
<p><b>Between-Isotype Sharing:</b> Same V+D+J genes (exact match)</p>

<h3><span class="vj">fastBCR-based</span></h3>
<p><b>Within-Isotype Clustering:</b> VJ genes + k-mer seeding (six 5-mers) + consensus ≥0.8</p>
<p><b>Between-Isotype Sharing:</b> Shared VJ + junction similarity % (configurable)</p>

<h2>🎯 Divergent Clone Detection</h2>

<p>Both methods identify <b>divergent clones</b> - sequences from an isotype that have diverged significantly 
from another isotype (e.g., IgE clones divergent from IgG1).</p>

<table>
<tr>
<th>Parameter</th>
<th>CompIgS</th>
<th>fastBCR-based</th>
</tr>
<tr>
<td>Shared cluster definition</td>
<td>Exact V+D+J match</td>
<td>Shared V+J + junction similarity %</td>
</tr>
<tr>
<td>Divergence threshold</td>
<td>CDR3 non-silent mutations ≥ N</td>
<td>Junction AA differences ≥ N</td>
</tr>
<tr>
<td>Default value</td>
<td>≥ 2 mutations</td>
<td>≥ 2 AA differences</td>
</tr>
<tr>
<td>Configurable range</td>
<td>1-3</td>
<td>1-10</td>
</tr>
</table>
""")
        layout.addWidget(content)
        
        self.setLayout(layout)


# =============================================================================
# MAIN INTEGRATED APPLICATION
# =============================================================================

class BCRAnalysisSuite(QMainWindow):
    """Main application window combining both analysis methods"""
    
    def __init__(self):
        super().__init__()
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle(f"{Config.APP_NAME} v{Config.VERSION}")
        self.setGeometry(100, 100, 1200, 900)
        
        # Apply dark theme
        self.apply_dark_theme()
        
        # Central widget
        central_widget = QWidget()
        main_layout = QVBoxLayout()
        
        
        # Tab widget
        self.tabs = QTabWidget()
        self.tabs.setStyleSheet("""
            QTabWidget::pane {
                border: 1px solid #444444;
                background-color: #2a2a2a;
                border-radius: 5px;
            }
            QTabBar::tab {
                background-color: #3a3a3a;
                color: #cccccc;
                padding: 12px 25px;
                margin-right: 3px;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                font-weight: bold;
            }
            QTabBar::tab:selected {
                background-color: #4CAF50;
                color: white;
            }
            QTabBar::tab:hover:!selected {
                background-color: #4a4a4a;
            }
        """)
        
        # Add tabs
        self.fastbcr_tab = FastBCRTab()
        self.tabs.addTab(self.fastbcr_tab, "📊 FastBCR (VJ + Junction)")
        
        self.compigs_tab = CompIgSTab()
        self.tabs.addTab(self.compigs_tab, "🧬 CompIgS (VDJ)")
        
        self.comparison_tab = ComparisonTab()
        self.tabs.addTab(self.comparison_tab, "📖 Method Comparison")
        
        main_layout.addWidget(self.tabs)
        
        # Bottom Button Layout (ADD THIS SECTION)
        bottom_button_layout = QHBoxLayout()
        bottom_button_layout.addStretch()  # Push button to the right
        
        # Open Output Folder Button
        open_folder_button = QPushButton("📁 Open Output Folder")
        open_folder_button.clicked.connect(self.open_default_output)
        open_folder_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                padding: 10px 20px;
                font-size: 12px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
        """)
        bottom_button_layout.addWidget(open_folder_button)
        
        # About Button
        about_button = QPushButton("ℹ️ About")
        about_button.clicked.connect(self.show_about)
        about_button.setStyleSheet("""
            QPushButton {
                background-color: #9C27B0;
                color: white;
                padding: 10px 20px;
                font-size: 12px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #7B1FA2;
            }
        """)
        bottom_button_layout.addWidget(about_button)
        
        # Exit Button
        exit_button = QPushButton("❌ Exit")
        exit_button.clicked.connect(self.close_application)
        exit_button.setStyleSheet("""
            QPushButton {
                background-color: #f44336;
                color: white;
                padding: 10px 20px;
                font-size: 12px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #d32f2f;
            }
        """)
        bottom_button_layout.addWidget(exit_button)
        
        main_layout.addLayout(bottom_button_layout)
        
        # Status bar
        self.statusBar().showMessage("Ready - Select an analysis method to begin")
        self.statusBar().setStyleSheet("background-color: #1a1a1a; color: #888888; padding: 5px;")
        
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)
    
    def apply_dark_theme(self):
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(42, 42, 42))
        palette.setColor(QPalette.WindowText, QColor(200, 200, 200))
        palette.setColor(QPalette.Base, QColor(30, 30, 30))
        palette.setColor(QPalette.AlternateBase, QColor(42, 42, 42))
        palette.setColor(QPalette.Text, QColor(200, 200, 200))
        palette.setColor(QPalette.Button, QColor(60, 60, 60))
        palette.setColor(QPalette.ButtonText, QColor(200, 200, 200))
        palette.setColor(QPalette.Highlight, QColor(76, 175, 80))
        palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
        self.setPalette(palette)
    
    def close_application(self):
        """Close the application with confirmation."""
        reply = QMessageBox.question(
            self,
            "Exit Application",
            "Are you sure you want to exit BCR Analysis Suite?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            # Stop any running threads first
            if hasattr(self, 'fastbcr_tab') and self.fastbcr_tab.r_thread:
                self.fastbcr_tab.r_thread.stop()
            
            # Force close
            QApplication.instance().quit()
            sys.exit(0)
    
    def open_default_output(self):
        """Open the default output folder."""
        import subprocess
        import platform
        
        # Default output location
        output_dir = os.path.join(os.path.expanduser("~"), "Downloads")
        
        # Check if CompIgS or FastBCR specific folders exist
        compigs_dir = os.path.join(output_dir, "CompIgS_Analysis")
        fastbcr_dir = os.path.join(output_dir, "FastBCR_Analysis")
        
        if os.path.exists(compigs_dir):
            output_dir = compigs_dir
        elif os.path.exists(fastbcr_dir):
            output_dir = fastbcr_dir
        
        try:
            if platform.system() == "Windows":
                os.startfile(output_dir)
            elif platform.system() == "Darwin":
                subprocess.run(["open", output_dir])
            else:
                subprocess.run(["xdg-open", output_dir])
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not open folder: {e}")
    
    def show_about(self):
        """Show about dialog."""
        about_text = f"""
        <h2 style="color: #4CAF50;">BCR Analysis Suite</h2>
        <p><b>Version:</b> {Config.VERSION}</p>
        <p><b>Description:</b><br>
        Integrated BCR Repertoire Analysis Platform combining two independent analysis methods:</p>
        <ul>
            <li><b>CompIgS (VDJ Matching):</b> Python-based analysis using exact V+D+J gene matching</li>
            <li><b>FastBCR (VJ + Junction):</b> R-based analysis using V+J matching with junction similarity</li>
        </ul>
        <p><b>Features:</b></p>
        <ul>
            <li>Shared/Unique clone identification</li>
            <li>Divergent clone detection</li>
            <li>Mutation analysis</li>
            <li>Comprehensive visualizations</li>
        </ul>
        <p style="color: #888888; font-size: 10px;">© 2024 BCR Analysis Suite</p>
        """
        
        msg = QMessageBox(self)
        msg.setWindowTitle("About BCR Analysis Suite")
        msg.setTextFormat(Qt.RichText)
        msg.setText(about_text)
        msg.setStyleSheet("""
            QMessageBox {
                background-color: #2a2a2a;
            }
            QMessageBox QLabel {
                color: #cccccc;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                padding: 8px 20px;
                border-radius: 4px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        msg.exec_()


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    
    # Set application-wide font
    font = QFont("Segoe UI", 10)
    app.setFont(font)
    
    window = BCRAnalysisSuite()
    window.show()
    
    sys.exit(app.exec_())


# In[ ]:





# In[ ]:




