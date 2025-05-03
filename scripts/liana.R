# Load required libraries ----
library(liana)
library(tidyverse)
library(Seurat)
library(qs)

# Interaction analysis setup ----
# Get consensus interactions and filter for BNB-related pairs
consensus_interactions <- select_resource("Consensus")$Consensus

markers <- read_csv(file.path("lookup", "markers.csv"))
bnb_genes <- markers$BNB[!is.na(markers$BNB)]

bnb_pairs <- consensus_interactions |>
  filter(source_genesymbol %in% bnb_genes | target_genesymbol %in% bnb_genes)

write_csv(
  bnb_pairs,
  file.path("results", "liana", "bnb_pairs.csv")
)

