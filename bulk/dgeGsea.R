library(magrittr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

rm(list = ls())

calcBulkDGE <- function(
    cm,
    df,
    varId,
    varDiff,
    varConf = NULL,
    preFilter = TRUE,
    convertCmToInt = TRUE
) {
    require(DESeq2)
    # gather samples to include
    samplesToIncl <- df[[varId]][!is.na(df[[varDiff]])]
    # filter cm and df
    df <- df[match(samplesToIncl, df[[varId]]),]
    cm <- cm[, match(samplesToIncl, colnames(cm))]
    # convert cm to integer
    if (convertCmToInt) cm %<>% apply(2, round)
    # calc DDS
    ddsDesign <- ifelse(is.null(varConf), paste0("~ ", varDiff), paste0("~ ", varConf, "+", varDiff))
    dds <- DESeqDataSetFromMatrix(countData = cm, colData = df, design = formula(ddsDesign))
    # low count filtering: keep genes with a count of 3 or higher in at least 1/10 samples 
    if (preFilter) dds <- dds[rowSums(counts(dds) >= 3) >= ncol(dds)/10,]
    # find differentially expressed genes 
    dds <- DESeq(dds)
    # get DGE table
    resLFC <- lfcShrink(dds, coef = resultsNames(dds)[resultsNames(dds) %>% length], type = "apeglm")
    return(resLFC)
}

calcGSEA <- function(
    dge,
    filterGenelist = F
) {
    require(clusterProfiler)
    require(msigdbr)
    # create gene list
    resSig <- dge %>% data.frame %>% tidyr::drop_na()
    if (filterGenelist) resSig %<>% filter(padj < 0.05)
    geneList <- setNames(resSig$log2FoldChange, rownames(resSig)) %>% sort(decreasing = TRUE)
    # get hallmark genesets
    geneSets <- msigdbr(species = "human", category = "H") %>% select(gs_name, gene_symbol, gs_description)
    hallGsea <- GSEA(gene = geneList, TERM2GENE = geneSets, pvalueCutoff = 1, minGSSize = 0, eps = 1e-50,)
    hallGsea@result$desc <- lapply(hallGsea@result$ID, function(gsID) geneSets$gs_description[geneSets$gs_name == gsID][1])
    return(hallGsea)
}

clinD <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/clinical/data/clinDf.RDS")
# gather oct vs ffpe info
salmon <- lapply(readRDS("data/Robj/salmon.RDS"), function(tissue) {
    rbind(tissue$normal$counts, tissue$lnc$counts[!rownames(tissue$lnc$counts) %in% rownames(tissue$normal$counts),])
})
# get oct samples and ffpe samples
samples <- sapply(salmon, colnames)
# laod cms
cms <- readRDS("data/Robj/cms.RDS")
# make a dataframe annotating tissue of origin
samplesDf <- data.frame(id = cms$nc$raw %>% colnames) %>%
    mutate(
        time = id %>% str_sub(1, 1) %>% recode("B" = "Baseline", "S" = "Surgery") %>% factor,
        patient = id %>% str_sub(3, 8) %>% factor,
        tissue = if_else(id %in% samples$ffpe, "FFPE", "OCT") %>% factor
    )
# add CCCA and UR
samplesDf$CCCA <- sapply(samplesDf$patient, function(pat) clinD$`Ki67 Response`[clinD$`Patient ID` == pat])
samplesDf$US <- sapply(samplesDf$patient, function(pat) clinD$`Ultrasound Response`[clinD$`Patient ID` == pat])
samplesDf$CCCA %<>% recode("Responder" = "Yes", "Non responder" = "No")
samplesDf$US %<>% recode("Responder" = "R", "Non responder" = "NR")
samplesDf %>% summary

dges <- list(
    time = list(
        All = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf,
            varId = "id",
            varDiff = "time",
            varConf = "tissue"
        ),
        CCCA_No = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf %>% filter(CCCA == "No"),
            varId = "id",
            varDiff = "time",
            varConf = "tissue"
        ),
        CCCA_Yes = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf %>% filter(CCCA == "Yes"),
            varId = "id",
            varDiff = "time",
            varConf = "tissue"
        ),
        US_NR = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf %>% filter(US == "NR"),
            varId = "id",
            varDiff = "time",
            varConf = "tissue"
        ),
        US_R = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf %>% filter(US == "R"),
            varId = "id",
            varDiff = "time",
            varConf = "tissue"
        )
    ),
    response = list(
        CCCA = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf,
            varId = "id",
            varDiff = "CCCA",
            varConf = "tissue"
        ),
        US = calcBulkDGE(
            cm = cms$nc$raw,
            df = samplesDf,
            varId = "id",
            varDiff = "US",
            varConf = "tissue"
        )
    )
)

saveRDS(dges, "data/Robj/dges.RDS")

gseas <- lapply(dges, function(dfSet) lapply(dfSet, calcGSEA))
gseas$time$All@result %>% names

lapply(gseas$response %>% names, function(set) {
    df <- gseas$response[[set]]@result %>% select(ID, NES, pvalue, p.adjust, qvalue)
    xlsx::write.xlsx(df, "results/gsea_resp.xlsx", sheetName = set, row.names = FALSE, append = TRUE)
})
