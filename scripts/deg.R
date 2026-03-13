########################################################
# differential gene expression analysis using pseudobulk
# requires running annotate.R first
########################################################

# libraries  ----
library(Seurat)
library(tidyverse)
library(qs)
library(limma)
library(DESeq2)
library(viridis)
library(Libra)
library(EnhancedVolcano)
library(readxl)
library(purrr)
library(writexl)

# functions ----
source(file.path("scripts", "dge_helper.R"))

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

lookup <- qs::qread(file.path("objects", "lookup.qs"), nthread = 4)

# Run DE analysis for CIDP vs CTRL in CSF
de_cidp_ctrl_csf <- performDEAnalysis(sc_merge, "CIDP", "CTRL", "CSF")

# Run DE analysis for GBS vs CTRL in CSF
de_gbs_ctrl_csf <- performDEAnalysis(sc_merge, "GBS", "CTRL", "CSF")

# Run DE analysis for CIDP vs CTRL in PBMC
de_cidp_ctrl_pbmc <- performDEAnalysis(sc_merge, "CIDP", "CTRL", "PBMC")

# Run DE analysis for GBS vs CTRL in PBMC
de_gbs_ctrl_pbmc <- performDEAnalysis(sc_merge, "GBS", "CTRL", "PBMC")
