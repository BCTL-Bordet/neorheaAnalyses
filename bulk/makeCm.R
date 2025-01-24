library(dplyr)
library(stringr)
library(magrittr)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(patchwork)

set.seed(1234)

# read salmon output and make a combined count matrix
cm <- lapply(readRDS("data/Robj/salmon.RDS"), function(tissue) {
    rbind(tissue$normal$counts, tissue$lnc$counts[!rownames(tissue$lnc$counts) %in% rownames(tissue$normal$counts),])
})

# find intersecting samples (all of them are baseline)
inter_samples <- intersect(colnames(cm$ffpe), colnames(cm$oct))
# for each of them select the tissue with both timepoints, or 
for (sample in inter_samples) {
    message(sample, "-------------")
    pat <- sample %>% str_sub(3, 8)
    in_cm <- lapply(cm, function(c) colnames(c)[colnames(c) %>% grepl(pat, .)])
    # if it has both timepoints in ffpe cm, delete the oct entry else delete the ffpe one
    if (length(in_cm$ffpe) > length(in_cm$oct)) {
        cm$oct <- cm$oct[, !colnames(cm$oct) %in% in_cm$oct]
        message("Deleting from OCT cm")
    } else {
        cm$ffpe <- cm$ffpe[, !colnames(cm$ffpe) %in% in_cm$ffpe]
        message("Deleting from FFPE cm")
    }
}
# combine both oct and ffpe in one cm
cms <- list(nc = list(raw = cbind(cm$oct, cm$ffpe)))

# cpm normalisation (edgeR)
cms$nc$cpm <- cpm(cms$nc$raw)
# vst normalisation (DESeq)
annData <- data.frame(sample = colnames(cms$nc$raw)) %>% 
    mutate(
        patient = sample %>% str_sub(3, 8),
        time = sample %>% str_sub(1, 1) %>% recode("B" = "Baseline", "S" = "Surgery") %>% factor,
        tissue = if_else(sample %in% colnames(cm$oct), "OCT", "FFPE") 
    )
cms$nc$vst <- cms$nc$raw %>% apply(2, as.integer) %>% DESeqDataSetFromMatrix(annData, design = ~ time) %>% vst %>% assay
rownames(cms$nc$vst) <- rownames(cms$nc$raw)

# remove batch effects for tissue and patient
removeBatch <- function(cm, annot, batchEffectVar, contrastVar){
    require(limma)
    mm <- model.matrix(~ annot[[contrastVar]], annot)
    mat <- removeBatchEffect(cm, batch = annot[[batchEffectVar]], design = mm)
    return(mat)  
}
cms$tc <- lapply(cms$nc, removeBatch, annData, "tissue", "time")
cms$pc <- lapply(cms$nc, removeBatch, annData, "patient", "time")
saveRDS(cms, "data/Robj/cms.RDS")

# plot some stuff to get an idea
figUmap <- lapply(cms %>% names, function(cor){
    lapply(cms[[cor]] %>% names, function(nor) {
        require(umap)
        cm <- cms[[cor]][[nor]]
        custom.config <- umap.defaults
        custom.config$random_state <- 0
        umap_res <- umap(data.frame(t(cm)), config = custom.config)
        df = data.frame(sample = annData$sample,
            umapX = umap_res$layout[,1], umapY = umap_res$layout[,2],
            sample_type = annData$time, wet = annData$tissue)
        ggplot(df, aes(x = umapX, y = umapY, color = wet, shape = sample_type)) + 
            geom_point() +
            labs(x = "", y = "", color = "", shape = "", title = paste(cor, nor)) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
                axis.text = element_blank(),
                axis.ticks = element_blank())
    }) 
})
figUmap <- wrap_plots(lapply(pUmap, wrap_plots), nrow = 3) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")

pMeanSD <- lapply(cms %>% names, function(cor){
    lapply(cms[[cor]] %>% names, function(nor) {
        cm <- cms[[cor]][[nor]]
        fig <- vsn::meanSdPlot(cm, plot = F)$gg +
            theme_classic() +
            labs(title =  paste(cor, nor), x = "", y = "") +
            theme(plot.title = element_text(hjust = 0.5))
    })
})
figMeanSD <- wrap_plots(lapply(pMeanSD, wrap_plots), nrow = 3) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")

# vst normalization looks always better
# patient and correction look very similar, maybe meadSd in bettern in tissue
