library(tidyr)
library(stringr)
library(magrittr)
library(tximport)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

set.seed(1234)

# read dfIds and initialize a files list
dfIds <- readRDS("data/dfIds.RDS")

# get salmon output files for oct samples
files <- list(
    oct = list(
        normal = dir("data/salmon/salmon_oct" , pattern = "quant.sf", recursive = TRUE, full.names = TRUE),
        lnc = dir("data/salmon/salmon_lncRNA_neoRhea" , pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
    ),
    ffpe = list(
        normal = dir("data/salmon/salmon_ffpe" , pattern = "quant.sf", recursive = TRUE, full.names = TRUE),
        lnc = dir("data/salmon/salmon_lncRNA_neoForma" , pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
    )
)

# add names
files$oct %<>% lapply(function(set) {
    names(set) <- dfIds$oct$sample[match(set %>% word(2, sep = "RHEA|Rhea-") %>% word(1, sep = "_") %>% paste0("NeoRhea_", .), dfIds$oct$id)]
    return(set)
})
files$ffpe %<>% lapply(function(set) {
    names(set) <- dfIds$ffpe$sample[match(set %>% word(4, sep = "/") %>% word(7, sep = "-") %>% word(1, sep = "_") %>% paste("NEOFORMA", .), dfIds$ffpe$id)]
    return(set)
})

# make an ensembl transcripts to entrez ids to gene symbols dataframe
k <- keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "TXNAME")
tx2entrez <- select(TxDb.Hsapiens.UCSC.hg38.knownGene, k, "GENEID", "TXNAME")
entrez2symbol <- select(org.Hs.eg.db, tx2entrez$GENEID, "SYMBOL", "ENTREZID")
tx2symbol <- data.frame(tx = tx2entrez$TXNAME, symbol = entrez2symbol$SYMBOL) %>% drop_na()

# read all salmon files
tx <- lapply(files, function(set) lapply(set, tximport, type = "salmon", tx2gene = tx2symbol))
saveRDS(tx, "data/Robj/salmon.RDS")
