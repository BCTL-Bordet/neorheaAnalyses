library(magrittr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(edgeR)
library(minfi)
library(ELMER)

set.seed(1234)
rm(list = ls())

# read expression data
exp <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/bulk_RNA/data/cm_combo.RDS")
genes <- exp$raw %>% rownames
exp <- exp$`vst, patient correction`
rownames(exp) <- genes; rm(genes)
# read methylation data
mSet <- readRDS("Robj/mSet.RDS")
# make variable same as exp colnames
mSet$sample <- paste0(mSet$Patient %>% str_sub(-5, -1), "_", mSet$Timepoint %>% str_sub(1, 1))
# patient NRHBE0025007 has both oct and ffpe samples, keep the ffpe one (oct is not timepoint matched)
mSet <- mSet[,colnames(mSet) != "NeoDia_63"]

# subset and sort using matched patients
commonSamples <- intersect(mSet$sample, exp %>% colnames) %>% sort
mPats <- commonSamples[commonSamples %>% str_sub(1,5) %>% duplicated] %>% str_sub(1,5)
commonSamples <- commonSamples[commonSamples %>% str_sub(1,5) %in% mPats]
exp <- exp[, match(commonSamples, colnames(exp))]
mSet <- mSet[, match(commonSamples, mSet$sample)]

dim(exp)
dim(mSet)

# get distal probes
distal.probes <- get.feature.probe(genome = "hg19", met.platform = "EPIC", promoter = F)
# get beta methylation values
met <- getBeta(mSet)
# set colnames same as exp
colnames(met) <- mSet$sample

all.equal(colnames(exp), colnames(met))

# transform gene symbols to ENSEMBL IDs
geneSymbToEnsembl <- function(cm){
  require(biomaRt)
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  genedesc <- getBM(attributes = c('external_gene_name','ensembl_gene_id'),
                  filters = 'external_gene_name', values = rownames(cm), mart = ensembl)
  unique(genedesc$external_gene_name) %>% length
  # filter out duplicated genes
  cm <- cm[rownames(cm) %in% unique(genedesc$external_gene_name),]
  # match order
  cm <- cm[match(genedesc$external_gene_name, rownames(cm)),]
  rownames(cm) <- genedesc$ensembl_gene_id
  return(cm)
}
exp[1:3,1:3]
exp %<>% geneSymbToEnsembl
exp[1:3,1:3]

annot <- data.frame(
    primary = mSet$sample,
    patient = mSet$sample %>% str_sub(1, 5),
    timepoint = mSet$Timepoint,
    tissue = mSet$Tissue
)

annot %>% head

# create and save elmer mae object
elmer.mae <- createMAE(
    exp = exp, met =  met,
    met.platform = "EPIC", genome = "hg19",
    filter.probes = distal.probes,
    colData = annot,
    linearize.exp = F, TCGA = F, save = F
)
# add responder resistor by US in colData
clinD <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/clinicoPathological/data/clinDf.RDS")
elmer.mae$response <- sapply(elmer.mae$patient, function(pat) clinD$`DR Result G`[clinD$`Patient ID` == pat]) %>% unname
elmer.mae$response <- ifelse(elmer.mae$response == "Response", T, F)
# save
saveRDS(elmer.mae, "Robj/elmer.mae.RDS")

elmer.mae_allprobes <- createMAE(
    exp = exp, met =  met,
    met.platform = "EPIC", genome = "hg19",
    colData = annot,
    linearize.exp = F, TCGA = F, save = F
)
# add responder resistor by US in colData
elmer.mae_allprobes$response <- sapply(elmer.mae_allprobes$patient, function(pat) clinD$`DR Result G`[clinD$`Patient ID` == pat]) %>% unname
elmer.mae_allprobes$response <- ifelse(elmer.mae_allprobes$response == "Response", T, F)
# save
saveRDS(elmer.mae_allprobes, "Robj/elmer.mae_allprobes.RDS")

rm(annot, clinD, commonSamples, distal.probes, exp, met, geneSymbToEnsembl, mSet)

# rebegin now with precreated mae object
elmer.mae_2 <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/bulk_RNA/ELMER/elmer.mae.RDS")
cm_combo <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/bulk_RNA/ELMER/cm_combo.RDS")
cm_combo %>% names

elmer.mae_2@colData
# change response from ki to ur
elmer.mae_2$response <- sapply(elmer.mae_2$Patient_ID %>% str_sub(-5, -1), function(pat) clinD$`DR Result G`[clinD$`Patient ID` == pat]) %>% unname

# ELMER mega-function ----
doELMER <- function(elmer.mae, variable, group1, group2, patientsVar,
                    probePval = 0.05, probeDeltaB = 0.1, methDiffMethod = "edgeR", 
                    pairPval = 0.05, enrichPval = 0.05, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.05,
                    resultsDir, cores = 4, comment = NULL, onlyOCT = F,
                    probeAnnot){
  # housekeeping ----
  plotsDir <- paste0(resultsDir, "/plots")
  dir.create(resultsDir, recursive = T)
  dir.create(plotsDir)
  # keep track of parameters
  log_df <- data.frame(parameter = c("variable", "group1", "group2", "patientsVar",
                                     "probePval", "probeDeltaB", "methDiffMethod",
                                     "pairPval", "enrichPval", "enrichLOR", "enrichMinInd",
                                     "topTFper", "resultDir", "cores", "comment"),
                       value = c(variable, group1, group2, patientsVar,
                                 probePval, probeDeltaB, methDiffMethod,
                                 pairPval, enrichPval, enrichLOR, enrichMinInd,
                                 topTFper, resultsDir, cores, comment))
  write.csv(log_df, file = paste0(resultsDir, "/call.csv"))
  # clinical data inspection ----
  df <- elmer.mae@colData %>% data.frame
  df$Sample_type <- df$primary <- NULL
  df <- df[!(duplicated(df)),]
  names(df) %<>% str_replace_all("\\.", " ")
  # identify differentially methylated probes ----
  if (methDiffMethod == "ELMER"){
    diff_hype <- get.diff.meth(data = elmer.mae, group.col = variable, group1 =  group1, group2 = group2,
                               minSubgroupFrac = 1, sig.dif = probeDeltaB, pvalue = probePval,
                               diff.dir = "hyper", # Search for hypermethylated probes in group 1
                               cores = cores, dir.out = paste0(resultsDir, "/diffMeth_SgHyper"))
    diff_hypo <- get.diff.meth(data = elmer.mae, group.col = variable, group1 =  group1, group2 = group2,
                               minSubgroupFrac = 1, sig.dif = probeDeltaB, pvalue = probePval,
                               diff.dir = "hypo", # Search for hypomethylated probes in group 1
                               cores = cores, dir.out = paste0(resultsDir, "/diffMeth_SgHypo"))
    # get near genes for each differentially methylated probe (10 upstream and 10 dowstream genes)
    nearGenes_hype <- GetNearGenes(data = elmer.mae, probes = diff_hype$probe, numFlankingGenes = 20)
    nearGenes_hypo <- GetNearGenes(data = elmer.mae, probes = diff_hypo$probe, numFlankingGenes = 20)
  } else {
    # verify using edgeR
    message("edgeR - differential probe analysis")
    getResTable <- function(betaMeth, annot, individual, factorOfInterest, contrast){
      # this is the factor of interest
      contrastVar <- annot[[factorOfInterest]] %>% factor
      # this is the individual effect that we need to account for
      individual <- annot[[individual]] %>% factor
      # use the above to create a design matrix
      design <- model.matrix(~ 0 + contrastVar + individual, data = annot)
      colnames(design) <- c(levels(contrastVar), levels(individual)[-1])
      # fit the linear model
      fit <- lmFit(betaMeth, design)
      # create a contrast matrix for specific comparisons
      cmd <- paste("contMatrix <- makeContrasts(", contrast, ", levels = design)", sep = '"')
      eval(parse(text = cmd)); rm(cmd)
      # fit the contrasts
      fit2 <- contrasts.fit(fit, contMatrix)
      fit2 <- eBayes(fit2)
      # get the table of results for the first contrast
      DMPs <- topTable(fit2, num = Inf, coef = 1) %>% S4Vectors::DataFrame()
      DMPs$FDR <- p.adjust(DMPs$P.Value, method = "fdr")
      # add delta Beta values for each probe
      f1 <- contrast %>% str_split("-") %>% unlist %>% magrittr::extract(1)
      f2 <- contrast %>% str_split("-") %>% unlist %>% magrittr::extract(2)
      met1 <- betaMeth[, annot[[factorOfInterest]] == f1] %>% rowMeans
      met2 <- betaMeth[, annot[[factorOfInterest]] == f2] %>% rowMeans
      metDif <- met1 - met2
      metDif <- metDif[match(rownames(DMPs), names(metDif))]
      DMPs$betaD <- metDif
      return(DMPs)
    }
    
    # get differential probe results correcting for individual patient effects
    met <- assays(getMet(elmer.mae))[[1]]
    if (!onlyOCT){
      met_OCT <- met[,colnames(met) %in% elmer.mae$primary[elmer.mae$tissue == "OCT"]]
      met_FFPE <- met[,colnames(met) %in% elmer.mae$primary[elmer.mae$tissue == "FFPE"]]
    }
    contExpr <- paste0(group1, "-", group2)
    DMPs <- list()
    DMPs$All <- getResTable(met, elmer.mae@colData, patientsVar, variable, contExpr)
    if (!onlyOCT){
      DMPs$OCT <- getResTable(met_OCT, elmer.mae@colData[elmer.mae$tissue == "OCT",], patientsVar, variable, contExpr)
      DMPs$FFPE <- getResTable(met_FFPE,  elmer.mae@colData[elmer.mae$tissue == "FFPE",],  patientsVar, variable, contExpr)
    }
    # save
    saveRDS(DMPs, paste0(resultsDir, "/DMPs.RDS"))
    message("edgeR done - get near genes & plot results")
    
    # plot differentially methylated probe volcanos ----
    plotVolc <- function(resDF, xThres = 0.3, yThres = probePval, yAxis = "FDR", yMin = 0,
                         xMin = -1, xMax = 1){
      resDF %<>% data.frame
      filteredGenes <- resDF %>% complete.cases %>% table %>% magrittr::extract(1) %>% unname
      resDF <- resDF[complete.cases(resDF),]
      upreg <- sum(resDF[[yAxis]] < yThres & resDF$betaD > xThres)
      upregP <- round(upreg/nrow(resDF)*100, 2)
      downreg <- sum(resDF[[yAxis]] < yThres & resDF$betaD < - xThres)
      downregP <- round(downreg/nrow(resDF)*100, 2)
      # add colors
      resDF$color <- ifelse(abs(resDF$betaD) > xThres & resDF[[yAxis]] < yThres, "red", "black") %>% unname
      # add y axis limits
      yMin <- ifelse(yMin == 0, min(resDF[[yAxis]]), yMin)
      # plot fig
      fig <- ggplot(resDF, aes(x = betaD, y = .data[[yAxis]])) +
        geom_point(color = resDF$color, alpha = 0.2, size = 0.9) +
        scale_x_continuous(limits = c(xMin, xMax)) +
        scale_y_continuous(limits = c(1, yMin), trans = ggforce::trans_reverser('log10')) +
        labs(x = expression(Delta ~ "Beta"), y = expression(-log[10]~FDR),
             caption = paste0(nrow(resDF), " distal probes\n",
                              "Hyper-methylated: ", upreg, " (", upregP, "%)",
                              "\nHypo-methylated: ", downreg, " (", downregP, "%)")) +
        theme_classic() + theme(plot.title = element_text(hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5),
                                legend.position = "none")
      return(fig)
    }

    yMin <- min(min(DMPs$All$FDR), min(DMPs$OCT$FDR), min(DMPs$FFPE$FDR))

    A <- plotVolc(DMPs$All, xThres = 0.1, yMin = yMin, xMin = -0.25, xMax = 0.25) +
      labs(x = "", title = "All samples", subtitle = paste("n =", ncol(met)))
    if (!onlyOCT){
      B <- plotVolc(DMPs$OCT, xThres = 0.1, yMin = yMin, xMin = -0.25, xMax = 0.25) +
        labs(y = "", x = expression(Delta ~ "Beta"), title = "OCT samples",
             subtitle = paste("n =", ncol(met_OCT)))
      C <- plotVolc(DMPs$FFPE, xThres = 0.1, yMin = yMin, xMin = -0.25, xMax = 0.25) +
        labs(y = "", x = "", title = "FFPE samples", subtitle = paste("n =", ncol(met_FFPE)))
      pVolc <- annotate_figure(ggarrange(A, B, C, ncol = 3), bottom = comment)
    } else pVolc <- A
    ggsave(paste0(plotsDir, "/volc_ALL.png"), pVolc, width = 15, height = 6.5)
    # get near genes for each differentially methylated probe (10 upstream and 10 dowstream genes)
    probes_hype <- DMPs$All %>% data.frame %>% filter(FDR < probePval & betaD > probeDeltaB) %>% rownames
    probes_hypo <- DMPs$All %>% data.frame %>% filter(FDR < probePval & betaD < -probeDeltaB) %>% rownames
    # get near genes for each differentially methylated probe (10 upstream and 10 dowstream genes)
    nearGenes_hype <- GetNearGenes(data = elmer.mae, probes = probes_hype, numFlankingGenes = 20)
    nearGenes_hypo <- GetNearGenes(data = elmer.mae, probes = probes_hypo, numFlankingGenes = 20)
  }
  
  # plot probe annotation pie ----
  probeAnnot %>% data.frame %>% dplyr::filter(Name %in% probes_hype) %>%
    dplyr::count(Relation_to_Island) %>%
    mutate(text = paste0(Relation_to_Island, " (n = ", n, ", ", round(n/length(probes_hype), 4)*100, "%)")) %>%
    mutate(text = factor(text)) %>%
    mutate(text = forcats::lvls_reorder(text, c(1, 3, 2, 4, 5, 6))) %>%
    ggplot(aes(x = "", y = n, fill = text)) +
    geom_bar(stat = "identity", position = position_fill()) +
    scale_fill_manual(values = wesanderson::wes_palette("Zissou1", 6, type = "continuous") %>% rev) +
    coord_polar(theta = "y", start = 0) +
    labs(fill = "", caption = "Shore = 0-2 kb from island, Shelf = 2-4 kb from island,\n
       N = upstream (5') of CpG island, S = downstream (3') of CpG island\n",
         title = "Hypermethylated distal probes relation to CpG islands") +
    theme_void() + theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotsDir, "/probesHype_relation_to_island.pdf"), width = 9, height = 6)
  
  probeAnnot %>% data.frame %>% dplyr::filter(Name %in% probes_hypo) %>%
    dplyr::count(Relation_to_Island) %>%
    mutate(., text = paste0(Relation_to_Island, " (n = ", n, ", ", round(n/length(probes_hypo), 4)*100, "%)")) %>%
    mutate(text = factor(text)) %>%
    mutate(text = forcats::lvls_reorder(text, c(1, 3, 2, 4, 5, 6))) %>%
    ggplot(aes(x = "", y = n, fill = text)) +
    geom_bar(stat = "identity", position = position_fill()) +
    scale_fill_manual(values = wesanderson::wes_palette("Zissou1", 6, type = "continuous") %>% rev) +
    coord_polar(theta = "y", start = 0) +
    labs(fill = "", caption = "Shore = 0-2 kb from island, Shelf = 2-4 kb from island,\n
       N = upstream (5') of CpG island, S = downstream (3') of CpG island\n",
         title = "Hypomethylated distal probes relation to CpG islands") +
    theme_void() + theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotsDir, "/probesHypo_relation_to_island.pdf"), width = 9, height = 6)
  
  # identify putative probe-gene pairs, enriched motifs and regulatory TFs ----
  pairs_hype <- get.pair(data = elmer.mae, group.col = variable, group1 =  group1, group2 = group2,
                         nearGenes = nearGenes_hype, mode = "supervised",
                         correlation = "negative", diff.dir = "hyper", raw.pvalue = pairPval,
                         filter.probes = T, filter.percentage = 0.05, filter.portion = 0.3,
                         dir.out = paste0(resultsDir, "/pairs_hype"), label = "hyper", cores = cores)
  
  pairs_hypo <- get.pair(data = elmer.mae, group.col = variable, group1 =  group1, group2 = group2,
                         nearGenes = nearGenes_hypo, mode = "supervised",
                         correlation = "negative", diff.dir = "hypo", raw.pvalue = pairPval,
                         filter.probes = T, filter.percentage = 0.05, filter.portion = 0.3,
                         dir.out = paste0(resultsDir, "/pairs_hypo"), label = "hypo", cores = cores)
  
  # plot distance from TSS for each case and filter out probe-genes pairs that are proximal
  rm_hype_pairs <- which(abs(pairs_hype$Distance) < 2000)
  rm_hypo_pairs <- which(abs(pairs_hypo$Distance) < 2000)
  pairs_hypo <- pairs_hypo[-rm_hypo_pairs,]
  pairs_hype <- pairs_hype[-rm_hype_pairs,]
  
  dfPair <- data.frame(label = c(rep("Hypomethylated", nrow(pairs_hypo)),
                                 rep("Hypermethylated", nrow(pairs_hype))),
                       dist = c(pairs_hypo$Distance, pairs_hype$Distance))
  A <- ggplot(dfPair, aes(y = dist, x = label)) + geom_violin() + geom_jitter() +
    labs(x = "Type of probe", y = "Distance in bases from paired gene") + theme_classic()
  B <- ggplot(dfPair, aes(y = dist, x = label)) + geom_violin() + geom_jitter() +
    labs(x = "Type of probe", y = "Distance in bases from paired gene (clipped)") +
    theme_classic() + scale_y_continuous(limits = c(-1e+04, 1e+04))
  ggarrange(A, B)
  ggsave(paste0(plotsDir, "/pairsDistance.pdf"), width = 15, height = 7)
  # identify enriched motif for significantly hypomethylated probes which have putative target genes
  enrich_hype <- get.enriched.motif(data = elmer.mae, probes = pairs_hype$Probe, 
                                    min.incidence = enrichMinInd, lower.OR = enrichLOR, pvalue = enrichPval,
                                    dir.out = paste0(resultsDir, "/enrich_hype"), label = "hyper")
  enrich_hypo <- get.enriched.motif(data = elmer.mae, probes = pairs_hypo$Probe, 
                                    min.incidence = enrichMinInd, lower.OR = enrichLOR, pvalue = enrichPval,
                                    dir.out = paste0(resultsDir, "/enrich_hypo"), label = "hypo")
  
  # identify regulatory TF for the enriched motifs
  TF_hype <- get.TFs(data = elmer.mae, group.col = variable, group1 =  group1, group2 = group2,
                     mode = "supervised", diff.dir = "hyper", correlation = "negative",
                     enriched.motif = enrich_hype, minSubgroupFrac = 1, topTFper = topTFper,
                     save.plots = T, dir.out = paste0(resultsDir, "/TF_hype"), label = "hyper", cores = cores)
  TF_hypo <- get.TFs(data = elmer.mae, group.col = variable, group1 =  group1, group2 = group2,
                     mode = "supervised", diff.dir = "hypo", correlation = "negative",
                     enriched.motif = enrich_hypo, minSubgroupFrac = 1, topTFper = topTFper,
                     save.plots = T, dir.out = paste0(resultsDir, "/TF_hypo"), label = "hypo", cores = cores)
  
  # additional plots ----
  # pairs heatmap plots
  message("Plot - pairs heatmap")
  heatmapPairs(data = elmer.mae, pairs = pairs_hype,
               group.col = variable, group1 =  group1, group2 = group2,
               filename = paste0(plotsDir, "/heat_hyperMeth_pairs.pdf"))
  heatmapPairs(data = elmer.mae, pairs = pairs_hypo,
               group.col = variable, group1 =  group1, group2 = group2,
               filename = paste0(plotsDir, "/heat_hypoMeth_pairs.pdf"))
  
  # motif enrichment plots
  message("Plot - enriched motifs OR")
  lol <- read.csv(paste0(resultsDir, "/enrich_hypo/getMotif.hypo.motif.enrichment.csv"))
  lol$motif %<>% str_split("_") %>% unlist %>% magrittr::extract(seq(1, nrow(lol)*2, 2))
  lol <- lol[,1:8] %>% arrange(desc(OR)) %>% filter(!duplicated(motif)) %>% head(10)
  lol$motif %<>% factor(levels = rev(lol$motif))
  ggplot(lol, aes(x = OR, y = motif)) +
    geom_errorbar(aes(xmin = lowerOR, xmax = upperOR), width = 0.2, color = "grey") +
    geom_point() + labs(y = "", x = "Odds ratio", title = "Top 10 motifs enriched\nin hypo-methylated probes",
                        caption = comment) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotsDir, "/motifs_hypoMeth.pdf"), width = 4, height = 6)
  
  lol <- read.csv(paste0(resultsDir, "/enrich_hype/getMotif.hyper.motif.enrichment.csv"))
  lol$motif %<>% str_split("_") %>% unlist %>% magrittr::extract(seq(1, nrow(lol)*2, 2))
  lol <- lol[,1:8] %>% arrange(desc(OR)) %>% filter(!duplicated(motif)) %>% head(10)
  lol$motif %<>% factor(levels = rev(lol$motif))
  ggplot(lol, aes(x = OR, y = motif)) +
    geom_errorbar(aes(xmin = lowerOR, xmax = upperOR), width = 0.2, color = "grey") +
    geom_point() + labs(y = "", x = "Odds ratio", title = "Top 10 motifs enriched\nin hyper-methylated probes",
                        caption = comment) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotsDir, "/motifs_hypeMeth.pdf"), width = 4, height = 6)
  
  # plot most sig transcription factors
  message("Plot - TF frequencies")
  # get all potential families
  families <- TF_hype$potential.TF.family
  families <- families[!is.na(families)]
  families %<>% str_split(";") %>% unlist %>% table %>% data.frame %>% arrange(Freq)
  colnames(families) = c("TF", "freq")
  families$freq %<>% factor
  
  ggplot(families, aes(x = TF, y = freq)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = families$TF) +
    coord_flip() +
    labs(x = "Transcription factor family", y = "Frequency of being potential master TF family for motif",
         title = "Master regulator TFs of motifs\nenriched in hyper-methylated probes",
         caption = comment) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotsDir, "/masterTF_hyper.pdf"), width = 5, height = 6)
  topFamsHype <- families$TF %>% as.character
  
  # plot most sig transcription factors
  # get all potential families
  families <- TF_hypo$potential.TF.family
  families <- families[!is.na(families)]
  families %<>% str_split(";") %>% unlist %>% table %>% data.frame %>% arrange(Freq)
  colnames(families) = c("TF", "freq")
  families$freq %<>% factor
  
  ggplot(families, aes(x = TF, y = freq)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = families$TF) +
    coord_flip() +
    labs(x = "Transcription factor family", y = "Frequency of being potential master TF family for motif",
         title = "Master regulator TFs of motifs\nenriched in hypo-methylated probes",
         caption = comment) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plotsDir, "/masterTF_hypo.pdf"), width = 5, height = 6)
  topFamsHypo <- families$TF %>% as.character
  
  # plot TF by expression for top families for hypermethylated motifs
  message("Plot - TF exp/meth scatter plots")
  dir.create(paste0(plotsDir, "/TF_hype"))
  dir.create(paste0(plotsDir, "/TF_hypo"))
  # plot TF by expression by motif that it regulates (hypo)
  hypo_plots <- lapply(names(enrich_hypo), function(motif){
    motifName <- motif %>% str_split("_") %>% sapply(., `[`, 1)
    # get probes of that motif
    probes <- enrich_hypo[[motif]]
    # find which TF is regulating that motif
    TF <- TF_hypo$top.potential.TF.family[TF_hypo$motif == motif]
    if (!is.na(TF) & !is.null(probes)){
      scatter.plot(data = elmer.mae, byTF = list(TF = TF, probe = probes), 
                   category = variable, save = F, lm_line = T) +
        labs(caption = comment, x = paste0("Avg DNA meththylation of probes in motif: ", motif))
      ggsave(paste0(plotsDir, "/TF_hypo/TF_", TF, "_", motif, ".pdf"), width = 6, height = 6)
    }
  })
  # plot TF by expression by motif that it regulates (hyper)
  hype_plots <- lapply(names(enrich_hype), function(motif){
    motifName <- motif %>% str_split("_") %>% sapply(., `[`, 1)
    # get probes of that motif
    probes <- enrich_hype[[motif]]
    # find which TF is regulating that motif
    TF <- TF_hype$top.potential.TF.family[TF_hype$motif == motif]
    if (!is.na(TF) & !is.null(probes)){
      scatter.plot(data = elmer.mae, byTF = list(TF = TF, probe = probes), 
                   category = variable, save = F, lm_line = T) +
        labs(caption = comment, x = paste0("Avg DNA meththylation of probes in motif: ", motif))
      ggsave(paste0(plotsDir, "/TF_hype/TF_", TF, "_", motif, ".pdf"), width = 6, height = 6)
    }
  })
  # heatmap of TFs and motifs ----
  message("Plot - TF/motif heatmaps")
  # get correlation between median enriched probe and exp of TF
  met <- getMet(elmer.mae) %>% assay
  exp <- getExp(elmer.mae) %>% assay
  
  motifProbesBetaHypo <- sapply(enrich_hypo, function(motifProbes){
    met[met %>% rownames %in% motifProbes,] %>% colMedians
  })
  motifProbesBetaHype <- sapply(enrich_hype, function(motifProbes){
    met[met %>% rownames %in% motifProbes,] %>% colMedians
  })
  
  # get expression of TFs
  ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  genedescHypo = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "external_gene_name", values = topFamsHypo, mart = ensembl)
  genedescHype = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "external_gene_name", values = topFamsHype, mart = ensembl)
  topFamsHypo <- genedescHypo$ensembl_gene_id[genedescHypo$ensembl_gene_id %in% rownames(exp)]
  topFamsHype <- genedescHype$ensembl_gene_id[genedescHype$ensembl_gene_id %in% rownames(exp)]
  
  expOfTFsHypo <- sapply(topFamsHypo, function(gene) exp[rownames(exp) == gene,])
  expOfTFsHype <- sapply(topFamsHype, function(gene) exp[rownames(exp) == gene,])
  
  # calcultate correlations between TF exp and motif met
  temp <- cbind(expOfTFsHypo, motifProbesBetaHypo)
  corMat <- cor(temp)
  corMat <- corMat[1:ncol(expOfTFsHypo), -c(1:ncol(expOfTFsHypo))] %>% t
  colnames(corMat) <- genedescHypo$external_gene_name[genedescHypo$ensembl_gene_id %in% rownames(exp)]
  pheatmap::pheatmap(corMat, main = "Hypomethylated motifs TFs heatmap",
           color = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(100),
           cluster_rows = F, cluster_cols = F,
           filename = paste0(plotsDir, "/heat_hypo_TFs.pdf"))
  
  temp <- cbind(expOfTFsHype, motifProbesBetaHype)
  corMat <- cor(temp)
  corMat <- corMat[1:ncol(expOfTFsHype), -c(1:ncol(expOfTFsHype))] %>% t
  colnames(corMat) <- genedescHype$external_gene_name[genedescHype$ensembl_gene_id %in% rownames(exp)]
  pheatmap::pheatmap(corMat, main = "Hypermethylated motifs TFs heatmap",
           color = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(100),
           cluster_rows = F, cluster_cols = F,
           filename = paste0(plotsDir, "/heat_hype_TFs.pdf"))
  
  # summary df ----
  # keep track of statistics
  sum_df <- data.frame(variable = c("Number of dysregulated probes",
                                    "Number of gene-probe pairs",
                                    "Number of enriched motifs",
                                    "Number of potential TF families"),
                       Hypomethylated = c(probes_hypo %>% length,
                                          pairs_hypo %>% nrow,
                                          enrich_hypo %>% length,
                                          topFamsHypo %>% length),
                       Hypermethylated = c(probes_hype %>% length,
                                           pairs_hype %>% nrow,
                                           enrich_hype %>% length,
                                           topFamsHype %>% length))
  write.csv(sum_df, file = paste0(resultsDir, "/summary.csv"))
}

ann <- getAnnotation(mSet)

# stringent criteria analysis
elmer.mae$patient %<>% paste0("p_", .) # hack for edgeR
doELMER(elmer.mae, variable = "timepoint", group1 = "Surgery", "Baseline", patientsVar = "patient",
        probePval = 0.01, probeDeltaB = 0.1, methDiffMethod = "edgeR",
        pairPval = 0.01, enrichPval = 0.01, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.01,
        resultsDir = "ELMER_results/strict", cores = 4, probeAnnot = ann,
        comment = "All 45 matched patients, strict p-values")

# responders analysis
elmer.mae_2_f <- elmer.mae[,elmer.mae$response]
doELMER(elmer.mae_2_f, variable = "timepoint", group1 = "Surgery", "Baseline", patientsVar = "patient",
        probePval = 0.01, probeDeltaB = 0.1, methDiffMethod = "edgeR",
        pairPval = 0.01, enrichPval = 0.01, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.01,
        resultsDir = "ELMER_results/strict_responders", cores = 4, probeAnnot = ann,
        comment = "Responder patients, strict p-values")

# non responders analysis
elmer.mae_2_f <- elmer.mae[,!elmer.mae$response]
doELMER(elmer.mae_2_f, variable = "timepoint", group1 = "Surgery", "Baseline", patientsVar = "patient",
        probePval = 0.01, probeDeltaB = 0.1, methDiffMethod = "edgeR",
        pairPval = 0.01, enrichPval = 0.01, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.01,
        resultsDir = "ELMER_results/strict_non_responders", cores = 4, probeAnnot = ann,
        comment = "Non responder patients, strict p-values")

# same but for all probes
elmer.mae_allprobes$patient %<>% paste0("p_", .) # hack for edgeR
doELMER(elmer.mae_allprobes, variable = "timepoint", group1 = "Surgery", "Baseline", patientsVar = "patient",
        probePval = 0.01, probeDeltaB = 0.1, methDiffMethod = "edgeR",
        pairPval = 0.01, enrichPval = 0.01, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.01,
        resultsDir = "ELMER_results/strict_allprobes", cores = 4, probeAnnot = ann,
        comment = "All 45 matched patients, strict p-values, all probes")

elmer.mae_2_f <- elmer.mae_allprobes[,elmer.mae_allprobes$response]
doELMER(elmer.mae_2_f, variable = "timepoint", group1 = "Surgery", "Baseline", patientsVar = "patient",
        probePval = 0.01, probeDeltaB = 0.1, methDiffMethod = "edgeR",
        pairPval = 0.01, enrichPval = 0.01, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.01,
        resultsDir = "ELMER_results/strict_responders_allprobes", cores = 4, probeAnnot = ann,
        comment = "Responder patients, strict p-values, all probes")

elmer.mae_2_f <- elmer.mae_allprobes[,!elmer.mae_allprobes$response]
doELMER(elmer.mae_2_f, variable = "timepoint", group1 = "Surgery", "Baseline", patientsVar = "patient",
        probePval = 0.01, probeDeltaB = 0.1, methDiffMethod = "edgeR",
        pairPval = 0.01, enrichPval = 0.01, enrichLOR = 1.1, enrichMinInd = 10, topTFper = 0.01,
        resultsDir = "ELMER_results/strict_non_responders_allprobes", cores = 4, probeAnnot = ann,
        comment = "Non responder patients, strict p-values, all probes")
