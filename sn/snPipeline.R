library(magrittr)
library(stringr)
library(dplyr)
library(Seurat)
library(Signac)

set.seed(1234)

source("sc_data/functionsSeuratSignac.R")
outDir <- "scratch/ulb/BCTL/apapagia/sc_RNA_ATAC/Robj"
prefix <- "neoS"

# create seurat RNA and ATAC objects from cellranger output
neoS <- list(
    RNA <- mergeRNA(cellrangerDirectory = "/scratch/ulb/BCTL/apapagia/sc_RNA_ATAC/cellrangerOut/", outDir, prefix),
    ATAC <- mergeATAC(cellrangerDirectory = "/scratch/ulb/BCTL/apapagia/sc_RNA_ATAC/cellrangerOut/", outDir, prefix)
)

# set metadata and rename cells for both objects
neoS %<>% lapply(setMetadata)
neoS %<>% lapply(renameCells)

# process both assays
neoS$RNA %<>% processRNA(outDir, prefix)
neoS$ATAC %<>% processATAC(outDir, prefix)

# combine in one object
neoS$Joint <- combineRNAandATAC(neoS$RNA, neoS$ATAC, outDir, prefix)

# add chromatin motif activity
neoS$ATAC %<>% chromVar(outDir, prefix)