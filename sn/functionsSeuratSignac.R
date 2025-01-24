mergeRNA <- function(
    cellrangerDirectory = NULL,
    outDirectory = NULL,
    prefix = NULL,
    addMitoPerc = TRUE,
    addRiboPerc = TRUE,
    filterSamples = TRUE,
    minCells = 400,
    maxCells = 16000,
    filterCells = TRUE,
    minFeature = 200,       # minimum number of expressed genes to keep a droplet
    maxFeature = 10000,     # maximum number of expressed genes to keep a droplet
    minCount = 400,         # minimum number of UMI per cell
    maxCount = 25000,       # maximum number of UMI per cell
    maxPercentMt = 2,       # minimum percentage of mitochondrial RNA to keep a droplet
    maxPercentRb = 1,       # minimum percentage of ribosomal RNA to keep a droplet
    saveRaw = TRUE
) {
    require(Seurat)
    require(magrittr)

    if (is.null(cellrangerDirectory)) stop("Set the cellranger output directory")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(prefix)) stop("Set prefix for saved objects")

    # find available samples
    samples <- list.dirs(cellrangerDirectory, recursive = FALSE) %>% basename

    seuratObj <- lapply(samples, function(sample) {
        s <- Read10X(data.dir = paste0(cellrangerDirectory, sample, "/outs/filtered_feature_bc_matrix/"))
        s %<>% .$`Gene Expression`
        s %<>% CreateSeuratObject(counts = ., project = sample)
    })

    seuratObj <- merge(seuratObj[[1]], y = seuratObj[2:length(seuratObj)])
    # add mitochondrial and ribosomal percent
    if (addMitoPerc) seuratObj %<>% PercentageFeatureSet("^MT-", , col.name = "percent.mt")
    if (addRiboPerc) seuratObj %<>% PercentageFeatureSet("^RP[SL]", col.name = "percent.rb")
    
    # save raw object
    if (saveRaw) saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_rna_raw.RDS"))

    # filter bad quality patients
    if (filterSamples) {
        badSamples <- which(seuratObj$sample %>% table < minCells | seuratObj$sample %>% table > maxCells) %>% names
        goodSamples <- setdiff(seuratObj$sample %>% unique %>% as.character, badSamples)
        seuratObj %<>% subset(sample %in% goodSamples)
    }

    if (filterCells) {
        # filter bad quality cells
        seuratObj <- seuratObj %>% subset(
            nCount_RNA < maxCount &
            nCount_RNA > minCount &
            nFeature_RNA > minFeature & 
            nFeature_RNA < maxFeature
        )
        if (addMitoPerc) seuratObj <- seuratObj %>% subset(percent.mt < maxPercentMt)
        if (addRiboPerc) seuratObj <- seuratObj %>% subset(percent.rb < maxPercentRb)
    }
    # save object
    saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_rna_filt.RDS"))
    return(seuratObj)
}

mergeATAC <- function(
    cellrangerDirectory = NULL,
    outDirectory = NULL,
    prefix = NULL,
    filterPeaks = TRUE,
    minPeakWidth = 20,
    maxPeakWidth = 10000,
    filterSamples = TRUE,
    minCells = 400,
    maxCells = 16000,
    filterCells = TRUE,
    minCount = 1000,       # minimum number of UMI per cell, 3K used in pbmc vignette, 1K in multiome vignette
    maxCount = 100000,     # maximum number of UMI per cell, 20K used in pbmc vignette, 100K in multiome vignette
    saveRaw = TRUE
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(GenomicRanges)
    require(Signac)
    require(Seurat)
    require(parallel)

    if (is.null(cellrangerDirectory)) stop("Set the cellranger output directory")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(prefix)) stop("Set prefix for saved objects")
    
    # find available samples
    samples <- list.dirs(cellrangerDirectory, recursive = FALSE) %>% basename

    # create a common peak set
    combined.peaks <- mclapply(samples, function(sample) {
        # read in peak sets
        file_to_read <- paste0(cellrangerDirectory, sample, "/outs/atac_peaks.bed")
        peaks <- read.table(file = file_to_read, col.names = c("chr", "start", "end"))
        # convert to genomic ranges
        peaks %<>% makeGRangesFromDataFrame
    }, mc.cores = 16, mc.preschedule = FALSE) %>% GRangesList

    # create a unified set of peaks to quantify in each dataset
    combined.peaks <- reduce(x = c(combined.peaks %>% unlist))
    # filter out bad peaks based on peak width
    if (filterPeaks) {
        peakwidths <- width(combined.peaks)
        combined.peaks <- combined.peaks[peakwidths  < minPeakWidth & peakwidths > minPeakWidth]
    }

    # load metadata
    metadata <- mclapply(samples, function(sample) {
        file_to_read <- paste0(cellrangerDirectory, sample, "/outs/per_barcode_metrics.csv")
        m <- read.table(file = file_to_read, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)
        # filter keeping droplets that are labelled as cells
        m %<>% filter(is_cell == 1)
    }, mc.cores = 16, mc.preschedule = FALSE)
    saveRDS(metadata, paste0(outDirectory, "/", prefix, "_metadata.RDS"))

    # create fragment objects
    frags <- mclapply(samples, function(sample) {
        file_to_read <- paste0(cellrangerDirectory, sample, "/outs/atac_fragments.tsv.gz")
        CreateFragmentObject(path = file_to_read, cells = rownames(metadata[[sample]]), validate.fragments = TRUE)
    }, mc.cores = 16, mc.preschedule = FALSE)

    # quantify peaks in each dataset
    counts <- mclapply(samples, function(sample) {
        FeatureMatrix(frags[[sample]], combined.peaks, rownames(metadata[[sample]]))
    }, mc.cores = 16, mc.preschedule = FALSE)

    # create chromatin assay
    c_assays <- mclapply(samples, function(sample) {
        CreateChromatinAssay(counts[[sample]], fragments = frags[[sample]])
    }, mc.cores = 16, mc.preschedule = FALSE)

    # remove samples with too many of too few cells
    if (filterSamples) {
        badSamples <- which(sapply(c_assays, ncol) < minCells | sapply(c_assays, ncol) > maxCells) %>% names
        goodSamples <- setdiff(names(c_assays), badSamples)
        c_assays <- c_assays[goodSamples]
        metadata <- metadata[goodSamples]
    }

    # create seurat object with metadata 
    seuratObj <- lapply(c_assays %>% names, function(sample) {
        s <- c_assays[[sample]] %>% CreateSeuratObject(assay = "ATAC")
        # add metadata
        meta <- metadata[[sample]]
        s@meta.data <- cbind(s@meta.data, meta)
        return(s)
    })

    # merge all samples
    seuratObj <- merge(seuratObj[[1]], y = seuratObj[2:length(seuratObj)])
    
    # save raw seurat object
    if (saveRaw) saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_atac_raw.RDS"))

    # filter bad quality cells
    if (filterCells) {
        seuratObj <- seuratObj %>% subset(
            nCount_ATAC < maxCount &
            nCount_ATAC > minCount
        )
    }

    # save seurat object
    saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_atac_filt.RDS"))
    return(seuratObj)
}

setMetadata <- function(
    seuratObj = NULL
) {
    require(stringr)
    require(magrittr)

    if (is.null(seuratObj)) stop("Input Seurat object to set metadata")

    seuratObj$sample <- seuratObj$orig.ident %>% factor
    seuratObj$orig.ident <- NULL
    seuratObj$patient <- seuratObj$sample %>% str_sub(1, 5) %>% factor
    seuratObj$timepoint <- ifelse(seuratObj$sample %>% str_sub(6) == "BL", "Baseline", "Surgery") %>% factor
    return(seuratObj)
}

renameCells <- function(
    seuratObj = NULL
) {
    require(Seurat)
    require(stringr)
    require(magrittr)

    if (is.null(seuratObj)) stop("Input Seurat object to rename cells")

    seuratObj %<>% RenameCells(
        new.names = paste0(
            seuratObj %>% colnames %>% str_split("-") %>% sapply(`[`, 1),
            "_",
            seuratObj$sample)
    )
    return(seuratObj)
}

processRNA <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    prefix = NULL,
    scale = TRUE,
    SCT = TRUE,
    splitSCT = TRUE,
    splitSCTby = NULL,
    integrate = TRUE,
    integrationSpit = NULL,
    integrationKweight = 50,
    varsToRegress = NULL, # c("percent.mt", "percent.rb", "nCount_RNA")
    dimsPCA = 20,
    clusterResolution = 0.5
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)

    if (is.null(seuratObj)) stop("Set the Seurat object to process")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(prefix)) stop("Set prefix for saved objects")

    # save a copy of the original seurat object
    sObj <- seuratObj

    if (scale) {
        seuratObj <- sObj
        DefaultAssay(seuratObj) <- "RNA"
        seuratObj %<>% FindVariableFeatures
        seuratObj %<>% ScaleData(vars.to.regress = varsToRegress)
        seuratObj %<>% RunPCA(features = VariableFeatures(object = seuratObj))
        seuratObj %<>% FindNeighbors(dims = 1:dimsPCA)
        seuratObj %<>% FindClusters(resolution = clusterResolution)
        seuratObj %<>% RunUMAP(dims = 1:dimsPCA)
        saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_rna_scale.RDS"))
    }

    if (SCT) {
        seuratObj <- sObj
        DefaultAssay(seuratObj) <- "RNA"
        seuratObj %<>% NormalizeData
        seuratObj %<>% SCTransform(method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = varsToRegress)
        seuratObj %<>% RunPCA
        seuratObj %<>% FindNeighbors(dims = 1:dimsPCA)
        seuratObj %<>% FindClusters(resolution = clusterResolution)
        seuratObj %<>% RunUMAP(dims = 1:dimsPCA)
        saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_rna_SCT.RDS"))
    }

    if (splitSCT) {
        if (is.null(splitSCTby)) stop("Set variable to split object before SCT")
        seuratObj <- sObj
        DefaultAssay(seuratObj) <- "RNA"
        seuratObj %<>% NormalizeData
        seuratObj %<>% SplitObject(split.by = splitSCTby)
        seuratObj %<>% lapply(SCTransform, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = varsToRegress)
        seuratObj <- merge(seuratObj[[1]], y = seuratObj[2:length(seuratObj)])
        DefaultAssay(seuratObj) <- "SCT"
        seuratObj %<>% RunPCA(features = VariableFeatures(object = seuratObj))
        seuratObj %<>% FindNeighbors(dims = 1:dimsPCA)
        seuratObj %<>% FindClusters(resolution = clusterResolution)
        seuratObj %<>% RunUMAP(dims = 1:dimsPCA)
        seuratObj %<>% PrepSCTFindMarkers
        saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_rna_SCT_split.RDS"))
    }

    if (integrate) {
        if (is.null(integrationSplit)) stop("Set variable to split object before integration")
        seuratObj <- sObj
        DefaultAssay(seuratObj) <- "RNA"
        seuratObj %<>% NormalizeData
        seuratObj %<>% SplitObject(split.by = integrationSplit)
        seuratObj %<>% lapply(SCTransform, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = varsToRegress)
        intFeatures <- SelectIntegrationFeatures(object.list = seuratObj, nfeatures = 10^4)
        seuratObj %<>% PrepSCTIntegration(anchor.features = intFeatures)
        seuratObj %<>% lapply(RunPCA, features = intFeatures, verbose = FALSE)
        seuratObjAnchors <- FindIntegrationAnchors(
            object.list = seuratObj, anchor.features = intFeatures,
            normalization.method = "SCT", reduction = "rpca", dims = 1:dimsPCA, verbose = T
        )
        seuratObj <- IntegrateData(
            anchorset = seuratObjAnchors,
            normalization.method = "SCT",
            dims = 1:dimsPCA, verbose = TRUE, k.weight = integrationKweight
        )
        seuratObj %<>% RunPCA
        seuratObj %<>% FindNeighbors(dims = 1:dimsPCA)
        seuratObj %<>% FindClusters(resolution = clusterResolution)
        seuratObj %<>% RunUMAP(dims = 1:dimsPCA)
        saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_rna_intRPCA.RDS"))
    }
    return(seuratObj)
}

processATAC <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    prefix = NULL,
    removeATACassay = TRUE,
    integrate = TRUE,
    integrationSpit = NULL,
    integrationKweight = 50,
    dimsPCA = 30,
    clusterResolution = 1.2
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(Signac)
    require(EnsDb.Hsapiens.v86)
    require(BSgenome.Hsapiens.UCSC.hg38)

    if (is.null(seuratObj)) stop("Set the Seurat object to process")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(prefix)) stop("Set prefix for saved objects")
    
    # get gene annotations for hg38
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevels(annotations) <- paste0("chr", seqlevels(annotations))

    # call peaks via MACS2
    peaks <- CallPeaks(seuratObj)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(seuratObj),
        features = peaks,
        cells = colnames(seuratObj)
    )

    # create a new assay using the MACS2 peak set and add it to the Seurat object
    seuratObj[["peaks"]] <- CreateChromatinAssay(
        counts = macs2_counts,
        fragments = Fragments(seuratObj),
        annotation = annotation
    )

    DefaultAssay(seuratObj) <- "peaks"
    if (removeATACassay) seuratObj[["ATAC"]] <- NULL

    # save raw seurat object
    seuratObjRaw <- seuratObj
    saveRDS(seuratObjRaw, paste0(outDirectory, "/", prefix, "_peaks_raw.RDS"))

    # process assay
    seuratObj %<>% FindTopFeatures(., min.cutoff = 5)
    seuratObj %<>% RunTFIDF
    seuratObj %<>% RunSVD
    seuratObj %<>% RunUMAP(reduction = "lsi", dims = 2:dimsPCA)
    seuratObj %<>% FindNeighbors(reduction = "lsi", dims = 2:dimsPCA)
    seuratObj %<>% FindClusters(algorithm = 3, resolution = clusterResolution, verbose = FALSE)

    # save seurat object
    saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_peaks.RDS"))

    # annotate peaks
    peak_annot <- ClosestFeature(object = seuratObj, regions = granges(seuratObj))
    saveRDS(peak_annot, paste0(outDirectory, "/", prefix, "_peak_annot.RDS"))

    if (integrate) {
        if (is.null(integrationSplit)) stop("Set variable to split object before integration")
        # split by factor
        seuratObj_split <- seuratObjRaw %>% SplitObject(split.by = integrationSplit)
        # compute LSI for each patient
        seuratObj_split %<>% lapply(function(s) {
            s %<>% FindTopFeatures(min.cutoff = 5)
            s %<>% RunTFIDF
            s %<>% RunSVD
            return(s)
        })
        # process the combined object
        seuratObj <- seuratObjRaw
        seuratObj %<>% FindTopFeatures(min.cutoff = 10)
        seuratObj %<>% RunTFIDF
        seuratObj %<>% RunSVD
        seuratObj %<>% RunUMAP(reduction = "lsi", dims = 2:30)
        # find integration anchors
        integrationAnchors <- FindIntegrationAnchors(
            object.list = seuratObj_split,
            anchor.features = rownames(seuratObj),
            reduction = "rlsi",
            dims = 2:dimsPCA
        )
        # integrate LSI embeddings
        integrated <- IntegrateEmbeddings(
            anchorset = integrationAnchors,
            reductions = seuratObj[["lsi"]],
            new.reduction.name = "integrated_lsi",
            dims.to.integrate = 1:dimsPCA,
            k.weight = integrationKweight
        )
        # process assay
        integrated %<>% RunUMAP(reduction = "integrated_lsi", dims = 2:30)
        integrated %<>% FindNeighbors(reduction = "integrated_lsi", dims = 2:30)
        integrated %<>% FindClusters(algorithm = 3, resolution = clusterResolution, verbose = FALSE)

        # save seurat object
        saveRDS(integrated, paste0(outDirectory, "/", prefix, "_peak_int.RDS"))
    }
    return(seuratObj)
}

combineRNAandATAC <- function(
    seuratObjRNA = NULL,
    seuratObjATAC = NULL,
    outDirectory = NULL,
    prefix = NULL,
    dimsPCA = 50,
    dimsLsi = 40,
    RNAisIntegrated = FALSE,
    ATACisIntegrated = FALSE
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(Signac)

    if (is.null(seuratObjRNA)) stop("Set the Seurat RNA object")
    if (is.null(seuratObjATAC)) stop("Set the Seurat ATAC object")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(prefix)) stop("Set prefix for saved objects")

    # find intersecting cells
    seuratObjRNA$cellID <- colnames(seuratObjRNA)
    seuratObjATAC$cellID <- colnames(seuratObjATAC)
    cells_joint <- intersect(seuratObjRNA$cellID, seuratObjATAC$cellID)
    # subset seurat objects
    seuratObjRNA %<>% subset(cellID %in% cells_joint)
    seuratObjATAC %<>% subset(cellID %in% cells_joint)

    seuratObj <- seuratObjRNA
    seuratObj[["peaks"]] <- seuratObjATAC@assays$peaks
    # add metadata
    seuratObj@meta.data <- cbind(seuratObjRNA@meta.data, seuratObjATAC@meta.data)
    # remove umap, add lsi reduction
    seuratObj@reductions$umap <- NULL
    if (!ATACisIntegrated) seuratObj@reductions$lsi <- seuratObjATAC@reductions$lsi
    else seuratObj@reductions$lsi <- seuratObjATAC@reductions$integrated_lsi

    # build a joint neighbor graph using both assays ----
    seuratObj %<>% FindMultiModalNeighbors(
        reduction.list = list("pca", "lsi"),
        dims.list = list(1:dimsPCA, 2:dimsLsi), verbose = TRUE
    )

    # build a joint UMAP visualization
    if (!RNAisIntegrated) seuratObj %<>% RunUMAP(nn.name = "weighted.nn", assay = "RNA", verbose = TRUE)
    else seuratObj %<>% RunUMAP(nn.name = "weighted.nn", assay = "integrated", verbose = TRUE)

    # save seurat object
    saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_joint_scale_lsi.RDS"))
}

chromVar <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    prefix = NULL,
    removePeakAssay = TRUE
) {
    require(magrittr)
    require(Seurat)
    require(Signac)
    require(EnsDb.Hsapiens.v86)
    require(TFBSTools)
    require(JASPAR2022)
    require(BSgenome.Hsapiens.UCSC.hg38)
    require(BiocParallel)

    register(MulticoreParam(8))

    if (is.null(seuratObj)) stop("Set the Seurat object to process")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(prefix)) stop("Set prefix for saved objects")

    # get a list of motif position frequency matrices from the JASPAR database
    pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
    saveRDS(pfm, paste0(outDirectory, "/", prefix, "_pfm.RDS"))
    # add motif information
    seuratObj %<>% AddMotifs(genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
    # compute motif activities
    seuratObj %<>% RunChromVAR(genome = BSgenome.Hsapiens.UCSC.hg38)
    DefaultAssay(seuratObj) <- "chromvar"
    # remove redundant assay
    if (removePeakAssay) seuratObj[["peaks"]] <- NULL
    # save seurat object
    saveRDS(seuratObj, paste0(outDirectory, "/", prefix, "_chromvar.RDS"))
    return(seuratObj)
}

plotQC <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    prefix = NULL,
    caption = NULL,
    RNA = TRUE,
    minFeature = 200,       # minimum number of expressed genes to keep a droplet
    maxFeature = 10000,     # maximum number of expressed genes to keep a droplet
    minCount = 400,         # minimum number of UMI per cell
    maxCount = 25000,       # maximum number of UMI per cell
    maxPercentMt = 2,       # minimum percentage of mitochondrial RNA to keep a droplet
    maxPercentRb = 1,       # minimum percentage of ribosomal RNA to keep a droplet
    ATAC = FALSE,
    minCountATAC = 1000,    # minimum number of UMI per cell, 3K used in pbmc vignette, 1K in multiome vignette
    maxCountATAC = 100000   # maximum number of UMI per cell, 20K used in pbmc vignette, 100K in multiome vignette

) {
    require(magritrr)
    require(reshape2)
    require(Seurat)
    require(ggplot2)
    require(ggpubr)
    require(patchwork)

    if (is.null(seuratObj)) stop("Set the Seurat object to plot")
    if (is.null(outDirectory)) stop("Set directory for saved figures")

    if (!dir.exists(outDirectory)) dir.create(outDirectory, recursive = TRUE)

    plot_cell_numbers <- function(seuratObj, caption) {
        fig <- table(seuratObj$sample) %>% data.frame %>% dplyr::mutate(patient = Var1 %>% str_sub(1, 5)) %>%
            mutate(time = ifelse(Var1 %>% str_sub(6, 7) == "BL", "Baseline", "Surgery")) %>%
            ggplot(., aes(x = patient, y = Freq, fill = time)) +
                geom_bar(stat = "identity", position = position_dodge()) +
                labs(
                    x = "Patient", y = "Number of cells", fill = "Timepoint",
                    subtitle = paste0(caption, " object, ", seuratObj %>% ncol, " cells in total")
                ) +
                theme_classic() +
                theme(
                    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
                    legend.position = "top",
                    plot.subtitle = element_text(hjust = 0.5)
                )
        # add table with cell numbers under the barplot
        tab <- table(seuratObj$patient, seuratObj$timepoint) %>% data.frame %>%
            acast(., Var1 ~ Var2) %>% t %>% ggtexttable
        # combine and return in one plot
        return(ggarrange(fig, tab, ncol = 1, heights = c(1, 0.2)))
    }
    
    f <- plot_cell_numbers(seuratObj, caption)
    ggsave(paste0(outDirectory, "/", prefix, "_cellNums.pdf"), f, width = 23, height = 5)

    if (RNA) {
        Idents(seuratObj) <- "patient"
        fig_v_rna_1 <- VlnPlot(seuratObj, features = c("nCount_RNA"), split.by = "timepoint", split.plot = FALSE, pt.size = 0) +
            scale_y_continuous(trans = "log10", breaks = c(10^2, 10^3, 10^5, minCount, maxCount)) +
            geom_hline(yintercept = minCount, color = "grey30", linetype = "dashed") +
            geom_hline(yintercept = maxCount, color = "grey30", linetype = "dashed") +
            labs(title = "Number of UMI counts per cell", x = "") +
            theme(axis.title.x = element_blank(), axis.text.x = element_blank())

        fig_v_rna_2 <- VlnPlot(seuratObj, features = c("nFeature_RNA"), split.by = "timepoint", split.plot = FALSE, pt.size = 0) +
            scale_y_continuous(trans = "log10", breaks = c(10^2, 10^3, 10^4, minFeature, maxFeature)) +
            geom_hline(yintercept = minFeature, color = "grey30", linetype = "dashed") +
            geom_hline(yintercept = maxFeature, color = "grey30", linetype = "dashed") +
            labs(title = "Number of identified genes per cell", x = "") +
            theme(axis.title.x = element_blank(), axis.text.x = element_blank())

        fig_v_rna_3 <- VlnPlot(seuratObj, features = c("percent.mt"), split.by = "timepoint", split.plot = FALSE, pt.size = 0) +
            scale_y_continuous(breaks = c(0, 50, 100, maxPercentMt)) +
            geom_hline(yintercept = maxPercentMt, color = "grey30", linetype = "dashed") +
            labs(title = "Percentage of mitochondrial genes per cell", x = "Patient") +
            theme(axis.title.x = element_blank(), axis.text.x = element_blank())

        fig_v_rna_4 <- VlnPlot(seuratObj, features = c("percent.rb"), split.by = "timepoint", split.plot = FALSE, pt.size = 0) +
            scale_y_continuous(breaks = c(0, 50, 100, maxPercentRb)) +
            geom_hline(yintercept = maxPercentRb, color = "grey30", linetype = "dashed") +
            labs(title = "Percentage of ribosomal genes per cell", x = "Patient")

        f <- (fig_v_rna_1 / fig_v_rna_2 / fig_v_rna_3 / fig_v_rna_4) +
            plot_layout(guides = "collect") & theme(legend.position = "bottom")

        ggsave(paste0(outDirectory, "/", prefix, "_violins_rna.pdf"), f, width = 14, height = 12)
    }
    
    if (ATAC) {
        f <- VlnPlot(seuratObj, features = c("nCount_ATAC"), split.by = "timepoint", split.plot = FALSE, pt.size = 0) +
            scale_y_continuous(trans = "log10", breaks = c(10^2, 10^3, 10^5, minCountATAC, maxCountATAC)) +
            geom_hline(yintercept = minCountATAC, color = "grey30", linetype = "dashed") +
            geom_hline(yintercept = maxCountATAC, color = "grey30", linetype = "dashed") +
            labs(title = "Number of UMI counts per cell", x = "")
        
        ggsave(paste0(outDirectory, "/", prefix, "_violins_atac.pdf"), f, width = 14, height = 5)
    }
    return(f)
}

getMarkerGenes <- function(
    David = TRUE,
    EORTC = TRUE,
    EORTCxlLocation = "~/Documents/Andreas/Bordet/Resources/MATTIA_single_cell_stuff.xlsx",
    Cibersort = TRUE,
    CibersortMatrixLocation = "~/Documents/Andreas/Bordet/Resources/LM22.txt"
) {
    require(magrittr)
    require(stringr)
    require(readxl)
    require(data.table)

    markerGenes <- list()
    if (David) {
        markerGenes$David <- list(
            `T cell` = c("CD3D", "CD3E"),
            `CAF` = "DCN",
            `epithelial luminal` = c("KRT19", "CD24"),
            `myeloid` = "LYZ",
            `epithelial luminal` = c("KRT19", "CD24"),
            `vascular endothelial` = c("VWF", "TFF3"), # TFF3-
            `pericyte` = "RGS5",
            `B cell` = "MS4A1", 
            `plasma cell` = "JCHAIN",
            `epithelial basal` = "KRT5",
            `MAST` = "CPA3",
            `pDC` = "LILRA4" 
        )

        markerGenes$David2 <- list(
            epithelial = c("CDH1", "EPCAM", "KRT8", "KRT5", "KRT19", "AURKA"),
            stromal = c("HTRA1", "FAP", "COL1A1", "PDGFRA", "PDGFRB", "ACTA2", "PECAM1"),
            immune = c("PTPRC", "CD68", "MS4A1", "JCHAIN", "CD3D", "CD8A", "FOXP3", "CD79A"),
            normal_breast = c("ANKRD30A", "SLPI", "KRT18", "KRT14", "ESR1", "MYLK", "KIT", "SOX10", "FOXA1", "ACTA2", "PROM1", "KRT19")
        )
    }

    if (EORTC) {
        # read EORTC table from Mattia
        markerGenes$EORTC <- read_excel(EORTCxlLocation) %>%
            select(cell_ontology, markers) %>% data.frame %>% filter(!is.na(markers))
            markers_EORTC <- split(markers_EORTC, markers_EORTC$cell_ontology)
            markers_EORTC %<>% lapply(., function(set) {
                set %>% select(markers) %>% as.character %>% str_split(", ") %>% unlist
            })
        # remove symbols from genes
        markerGenes$EORTC %<>% lapply(function(set) {
            sapply(set, function(gene) {
                gene %>% str_remove_all("\"") %>% str_remove_all("c\\(") %>% str_remove_all("\\)")
            }) %>% unname
        })
    }

    if (Cibersort) {
        # get list of CIBERSORT markers
        lm22 <- data.fread(CibersortMatrixLocation)
        cTypes <- colnames(lm22)[-1] %>% setNames(.,  .)
        markerGenes$Cibersort <- lapply(cTypes, function(ctype) {
            vector <- lm22[[ctype]]
            genes <- lm22$`Gene symbol`
            # replace some gene names
            genes[genes == "FAIM3"] <- "FCMR"
            genes[genes == "EMR1"] <- "ADGRE1"
            df <- data.frame(vector, genes)
            # remove "GUSBP11", "IGLL3P", "LOC100130100"
            df$vector[df$genes == "GUSBP11"] <- 0
            df$vector[df$genes == "IGLL3P"] <- 0
            df$vector[df$genes == "LOC100130100"] <- 0
            top10genes <- df %>% arrange(-vector) %>% head(10) %>% select(genes) %>% unlist %>% unname
        })

    }
    
    return(markerGenes)
}

plotUmapMarkers <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    markerGeneList = NULL,
    cutoff = "q0"
) {
    require(Seurat)
    require(ggplot2)
    require(magritrr)

    if (is.null(seuratObj)) stop("Set the Seurat object to plot")
    if (is.null(outDirectory)) stop("Set directory for saved figures")
    if (is.null(markerGeneList)) stop("Set gene marker list")

    if (!dir.exists(outDirectory)) dir.create(outDirectory, recursive = TRUE)

    f <- lapply(markerGeneList %>% names, function(genesetName) {
        message(genesetName, " ------------")
        dir.create(paste0(outDirectory, "/", genesetName))
        genesets <- markerGeneList[[genesetName]]
        # iterate over gene sets
        lapply(genesets %>% names, function(geneSet) {
            # remove genes not in assay
            na_genes <- genesets[[geneSet]] %>% unlist %>% unname
            na_genes <- na_genes[!na_genes %in% rownames(seuratObj)] %>% unique
            gS <- genesets[[geneSet]] %>% lapply(function(set) set[!(set %in% na_genes)])
            gS <- gS[which(gS %>% sapply(length) != 0)]
            # iterate over genes
            lapply(gS, function(gene) {
                fileToSave <- paste0(outDirectory, "/", genesetName, "/", geneSet, " ", gene, ".png")
                if (!file.exists(fileToSave)) {
                    f <- FeaturePlot(seuratObj, gene, reduction = "umap", min.cutoff = cutoff) +
                        NoAxes() + NoLegend() +
                        labs(caption = paste("cell ontology:", geneSet)) +
                        theme(plot.caption = element_text(hjust = 0))
                    ggsave(fileToSave, f, width = 8, height = 6)
                }
            })
        })
    })
    return(f)
}

plotUmapMetadata <- function(
    seuratObj = NULL,
    outDirectory = NULL
) {
    require(Seurat)
    require(ggplot2)
    require(magritrr)
    require(Polychrome)

    if (is.null(seuratObj)) stop("Set the Seurat object to plot")
    if (is.null(outDirectory)) stop("Set directory for saved figures")
    
    if (!dir.exists(outDirectory)) dir.create(outDirectory, recursive = TRUE)
    
    s <- seuratObj
    dir <- outDirectory
    # plot sample
    Idents(s) <- "sample"
    colorsToUse <- createPalette(s$sample %>% unique %>% length, c("#010101")) %>% unname
    f1 <- DimPlot(s, label = FALSE, cols = colorsToUse) + ggtitle("Samples") + theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(dir, "/samples.png"), f, width = 10, height = 6)
    # plot patient
    Idents(s) <- "patient"
    colorsToUse <- createPalette(s$patient %>% unique %>% length, c("#010101")) %>% unname
    f2 <- DimPlot(s, label = FALSE, cols = colorsToUse) + ggtitle("Patients") + theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(dir, "/patients.png"), f, width = 8, height = 6)
    # plot timepoint
    Idents(s) <- "timepoint"
    f3 <- DimPlot(s, label = FALSE) + ggtitle("Timepoint") + theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(dir, "/timepoint.png"), f, width = 7, height = 6)
    # plot clusters
    if ("seurat_clusters" %in% names(s@meta.data)) {
        Idents(s) <- "seurat_clusters"
        colorsToUse <- createPalette(s$seurat_clusters %>% unique %>% length, c("#010101")) %>% unname
        f4 <- DimPlot(s, label = TRUE, cols = colorsToUse) + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0(dir, "/clusters.png"), f, width = 8, height = 6)
    } else {
        s$peaks_snn_res.1.2 %<>% factor(levels = 0:(s$peaks_snn_res.1.2 %>% unique %>% length - 1))
        Idents(s) <- "peaks_snn_res.1.2"
        colorsToUse <- createPalette(s$peaks_snn_res.1.2 %>% unique %>% length, c("#010101")) %>% unname
        f4 <- DimPlot(s, label = TRUE, cols = colorsToUse) + ggtitle("Clusters Peaks") + theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0(dir, "/clusters_peaks.png"), f, width = 8, height = 6)
        if ("RNA_snn_res.0.5" %in% names(s@meta.data)) {
            s$RNA_snn_res.0.5 %<>% factor(levels = 0:(s$RNA_snn_res.0.5 %>% unique %>% length - 1))
            Idents(s) <- "RNA_snn_res.0.5"
            colorsToUse <- createPalette(s$RNA_snn_res.0.5 %>% unique %>% length, c("#010101")) %>% unname
            f <- DimPlot(s, label = TRUE, cols = colorsToUse) + ggtitle("Clusters RNA") + theme(plot.title = element_text(hjust = 0.5))
            ggsave(paste0(dir, "/clusters_rna.png"), f, width = 8, height = 6)
        }
        else if ("SCT_snn_res.0.5" %in% names(s@meta.data)) {
            s$SCT_snn_res.0.5 %<>% factor(levels = 0:(s$SCT_snn_res.0.5 %>% unique %>% length - 1))
            Idents(s) <- "SCT_snn_res.0.5"
            colorsToUse <- Polychrome::createPalette(s$SCT_snn_res.0.5 %>% unique %>% length, c("#010101")) %>% unname
            f <- DimPlot(s, label = TRUE, cols = colorsToUse) + ggtitle("Clusters RNA") + theme(plot.title = element_text(hjust = 0.5))
            ggsave(paste0(dir, "/clusters_rna.png"), f, width = 8, height = 6)
        }
    }
    return(list(sample = f1, patient = f2, timepoint = f3, cluster = f4))
}

plotCelltypes <- function(
    seuratObj,
    outDirectory = NULL,
    prefix = NULL
) {
    require(dplyr)
    require(magrittr)
    require(ggplot2)
    require(Seurat)

    if (is.null(seuratObj)) stop("Set the Seurat object to plot")
    if (is.null(outDirectory)) stop("Set directory for saved figures")
    if (is.null(markerGeneList)) stop("Set gene marker list")

    if (!dir.exists(outDirectory)) dir.create(outDirectory, recursive = TRUE)

    plotSampleCelltypeBarplot <- function(sMeta, celltypeVar) {
        sMeta %>% count(sample, .data[[celltypeVar]]) %>% group_by(sample) %>% mutate(freq = n / sum(n)) %>%
            ggplot(., aes(x = sample, y = freq, fill = celltype)) +
            geom_bar(stat = "identity") +
            labs(x = "Sample", y = "", fill = "Cell type") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "top")
    }

    Idents(seuratObj) <- "celltype"
    f1 <- DimPlot(seuratObj)
    ggsave(paste0(outDirectory, "/", prefix, "_celltypes_umap.png"), f, width = 8, height = 6)
    Idents(s) <- "celltype_3"
    f2 <- DimPlot(seuratObj)
    ggsave(paste0(outDirectory, "/", prefix, "_celltypes_3_umap.png"), f, width = 8, height = 6)
    f3 <- plotSampleCelltype(seuratObj@meta.data, "celltype")
    ggsave(paste0(outDirectory, "/", prefix, "_celltypes_samples_barplot.png"), f, width = 9, height = 3)
    f4 <- plotSampleCelltype(seuratObj@meta.data, "celltype_3")
    ggsave(paste0(outDirectory, "/", prefix, "_celltypes_3_samples_barplot.png"), f, width = 9, height = 3)
    
    return(list(umap1 = f1, umap2 = f2, barplot1 = f3, barplot2 = f4))
}

infercnvCreateReference <- function(
    seuratObj = NULL,
    seuratObjMetadata = NULL,
    celltypeColumnName = NULL,
    targetNumberOfCells = 500,
    markersImmune = c("MS4A1", "ITGAX", "CD2"),
    markersStromal = c("FAP", "PDGFRB", "PECAM1"),
    saveRDS = FALSE,
    outDirectory = NULL
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(infercnv)


    if (is.null(seuratObj)) stop("Set the Seurat object")
    if (is.null(seuratObjMetadata)) stop("Set the Seurat object metadata")
    if (is.null(celltypeColumnName)) stop("Set the celltype metadata column name")
    if (saveRDS & is.null(outDirectory)) stop("Set directory to save infercnv reference")

    # add celltype info
    seuratObj$celltype <- seuratObjMetadata[[celltypeColumnName]] %>% as.character

    # subset object including only immune and stromal baseline cells ----
    seuratObj_BL_imstr <- seuratObj %>% subset(celltype %in% c("Immune", "Stromal") & timepoint == "Baseline")
    # create a table with number of cells in celltype by patient
    cellNums <- table(seuratObj_BL_imstr$sample, seuratObj_BL_imstr$celltype) %>%
        as.data.frame %>% reshape2::dcast(Var1 ~ Var2)
    rownames(cellNums) <- cellNums$Var1
    cellNums$Var1 <- NULL
    # now multiply such that sum equals to 500 (target reference cells)
    cellNums %<>% divide_by(cellNums %>% sum) %>% multiply_by(targetNumberOfCells) %>% round

    # remove samples with zero cells
    rmSamples <- cellNums[cellNums$Immune == 0 & cellNums$Stromal == 0,] %>% rownames
    cellNums <- cellNums[!(rownames(cellNums) %in% rmSamples),]; rm(rmSamples)

    # split by sample
    seuratObj_BL_imstr_splt <- seuratObj_BL_imstr %>% SplitObject(split.by = "sample")
    # remove samples with no immune and stromal cells
    seuratObj_BL_imstr_splt <- seuratObj_BL_imstr_splt[names(seuratObj_BL_imstr_splt) %in% rownames(cellNums)]

    # find best scoring immune and stromal cells for each sample
    cells_imstr_strat <- lapply(seuratObj_BL_imstr_splt %>% names, function(sName){
        message(sName)
        # get assay
        s <- seuratObj_BL_imstr_splt[[sName]]
        # find how many immune and stromal cells we want from each sample
        cellsImmune <- cellNums[rownames(cellNums) == sName, 1]
        cellsStromal <- cellNums[rownames(cellNums) == sName, 2]
        cell_im <- cell_st <- NULL
        # immune cells
        if (cellsImmune != 0) {
            s_im <- s %>% subset(celltype == "Immune")
            lol <- s_im@assays$RNA@data[rownames(s_im) %in% markersImmune,] %>% colSums
            cell_im <- lol %>% sort(decreasing = TRUE) %>% names %>% head(cellsImmune)
        }
        # stromal cells
        if (cellsStromal != 0) {
            s_st <- s %>% subset(celltype == "Stromal")
            lol <- s_st@assays$RNA@data[rownames(s_st) %in% markersStromal,] %>% colSums
            cell_st <- lol %>% sort(decreasing = TRUE) %>% names %>% head(cellsStromal)
        }
        return(list(immune = cell_im, stromal = cell_st))
    }) %>% unlist %>% unname

    # create reference set seurat object
    seuratObj$cellid <- colnames(seuratObj)
    ref <- seuratObj %>% subset(cellid %in% targetNumberOfCells)
    if (saveRDS) saveRDS(ref, paste0(outDirectory, "/infercnv_reference.RDS"))
    return(ref)
}

infercnvCreateObject <- function(
    seuratObj = NULL,
    referenceObj = NULL,
    geneOrderFile = "~/Documents/Andreas/Bordet/Projects/neoRhea/sc_RNA_ATAC/resources/refdata-gex-GRCh38-2020-A_gen_pos.ordered.txt",
    saveRDS = TRUE,
    outDirectory = NULL
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(infercnv)

    if (is.null(seuratObj)) stop("Set the Seurat object")
    if (is.null(referenceObj)) stop("Set the infercnv reference object")
    if (saveRDS & is.null(outDirectory)) stop("Set directory to save infercnv object")
    
    # create observation set (only epithelial cells) seurat object
    seuratObjSplit <- seuratObj %>% subset(celltype == "Epithelial") %>% SplitObject(split.by = "sample")
    # add all cells that are not epithelial from all patients in one object
    seuratObjSplit[["Immune"]] <- seuratObj %>% subset(celltype %in% c("Immune"))
    # split stromal to baseline and surgery, if all together infercnv won't compute :/
    seuratObjSplit[["Stromal_B"]] <- seuratObj %>% subset(celltype %in% c("Stromal") & timepoint == "Baseline")
    seuratObjSplit[["Stromal_S"]] <- seuratObj %>% subset(celltype %in% c("Stromal") & timepoint == "Surgery")
    
    gene_order <- read.table(geneOrderFile, row.names = 1)

    # prepare infercnv input per sample
    preInfercnv <- lapply(seuratObjSplit, function(seuratObjchunk){
        message(seuratObjchunk$sample %>% unique)
        # merge seurat objects
        message("Merging patient and normal reference objects")
        tempS <- merge(seuratObjchunk, referenceObj)
        tempS$origin <- c(rep("Observation", ncol(seuratObjchunk)), rep("Reference", ncol(referenceObj)))
        # make cell type annotations vector
        annot_f <- data.frame(V2 = tempS$origin, row.names = tempS %>% colnames)
        # create infercnv object
        message("Creating infercnv object")
        infercnv_obj <- CreateInfercnvObject(
            raw_counts_matrix = tempS@assays$RNA@counts,
            annotations_file = annot_f,
            gene_order_file = gene_order,
            ref_group_names = c("Reference"))
        return(infercnv_obj)
    })

    if (saveRDS) saveRDS(ref, paste0(outDirectory, "/infercnv_object.RDS"))
    return(preInfercnv)

}

infercnvRun <- function(
    infercnvObject = NULL,
    outDirectory = NULL,
    prefix = NULL
) { 
    if (is.null(infercnvObject)) stop("Set the infercnv object")
    if (is.null(outDirectory)) stop("Set output directory")

    require(magrittr)
    require(infercnv)

    # knn classifier to split cells to CNV abberated and CNV clean
    classifyInfercnv <- function(infercnvDir, thres = 0.05){
        # thres is percentage of cells to include in the instable profile
        require(class)
        require(cola)
        require(scales)
        # read and cluster reference cells
        ref.data <- read.table(paste0(infercnvDir, "/infercnv.references.txt"), header = TRUE, row.names = 1, check.names = FALSE)
        ref.data.scaled <- ref.data %>% data.matrix %>% rescale(to = c(-1, 1))
        mat <- adjust_matrix(ref.data.scaled)  # optional
        message("Running concensus partition on reference")
        res <- consensus_partition(mat, partition_method = "skmeans", max_k = 5) # 7mins
        # find optimal number of clusters
        k <- suggest_best_k(res)
        groups <- res@object_list[[k-1]]
        # extract cluster of each observation cell
        clusters <- groups$class_df$class
        # get profiles for each of the reference clusters
        ref_profiles <- lapply(clusters %>% unique %>% sort, function(cluster){
            cellsInCluster <- groups$class_df %>% subset(class == cluster) %>% rownames
            ref_data <- ref.data.scaled[, colnames(ref.data.scaled) %in% cellsInCluster] %>% rowMeans
        })
        names(ref_profiles) <- clusters %>% unique %>% sort
        # get observation cells
        obs.data <- read.table(paste0(infercnvDir, "/infercnv.observations.txt"), header = TRUE, row.names = 1, check.names = FALSE)
        obs.data.scaled <- obs.data %>% data.matrix %>% scales::rescale(to = c(-1, 1))
        # use the top 5 percent of observation cells with highest CNV score to compute the average CNV profile
        genomic.instability.scores <- apply(obs.data.scaled, 2, function(x) mean(x^2) )
        # get cells above thres
        thres <- 1 - thres
        groupObs <- which(genomic.instability.scores > quantile(genomic.instability.scores, thres)) %>% names
        # get profile of observation cells
        profObs <- obs.data.scaled[, colnames(obs.data.scaled) %in% groupObs] %>% rowMeans
        # concatenate reference and super instable cell count matrices
        trainData <- cbind(ref.data.scaled, obs.data.scaled[, colnames(obs.data.scaled) %in% groupObs])
        # concatenate reference and super instable cell profiles
        annotTrain <- cbind(do.call(cbind, ref_profiles), profObs)
        # get profile of observation cells
        profObs <- obs.data.scaled[, colnames(obs.data.scaled) %in% groupObs] %>% rowMeans
        # concatenate reference and super instable cell count matrices
        trainData <- cbind(ref.data.scaled, obs.data.scaled[, colnames(obs.data.scaled) %in% groupObs])
        # concatenate reference and super instable cell profiles
        annotTrain <- cbind(do.call(cbind, ref_profiles), profObs)
        # give cluster names
        annot <- c(clusters %>% unique %>% sort %>% paste0("Ref_", .), "Tumor")
        # knn implementation
        pr <- knn(annotTrain %>% t, obs.data.scaled %>% t, cl = annot %>% t, k = 1)
        names(pr) <- colnames(obs.data)
        return(pr)
    }

    dir.create(outDirectory, recursive = TRUE)
    # run infercnv for each of the list elements separately
    infer_res <- lapply(infercnvObj, function(iObj) {
        infercnv_obj_default <- infercnv::run(
            iObj, cutoff = 0.1, out_dir = outDirectory,
            cluster_by_groups = FALSE, plot_steps = FALSE, denoise = TRUE, save_rds = FALSE,
            analysis_mode = "subclusters",
            HMM = TRUE, no_prelim_plot = FALSE, png_res = 600, num_threads = 8
        )
        # classify cells as tumor or non-tumor
        pr <- classifyInfercnv(outDirectory, thres = 1)
        saveRDS(pr, paste0(outDirectory, "/knn_predictions.RDS"))
        return(pr)
    })
    return(infer_res)
}

dge <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    suffix = NULL,
    barcodesToInclude = NULL,
    identToCompare = NULL,
    ident1 = NULL,
    logFCthres = 0,
    minPerc = 0
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(future)

    if (is.null(seuratObj)) stop("Set the Seurat object to process")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(suffix)) stop("Set suffix for saved objects")
    if (is.null(barcodesToInclude)) stop("Set the cell barcodes to include")

    plan("multicore", workers = 8)
    options(future.globals.maxSize = 5000 * 1024^2)

    # filter out not included barcodes
    seuratObj$barcode <- seuratObj %>% colnames
    seuratObj <- seuratObj %>% subset(barcode %in% barcodesToInclude)
    # set ident
    Idents(seuratObj) <- identToCompare
    # check if there are enough cells in each ident
    cellNums <- table(seuratObj[[identToCompare]])
    if (any(cellNums <= 3) | length(cellNums) == 1) dge <- "Not enough cells" %>% data.frame
    else {
        # do DGE
        dge <- FindMarkers(
            seuratObj,
            ident.1 = ident1,
            logfc.threshold = logFCthres,
            min.pct = minPerc
        )
    }
    saveRDS(dge, paste0(outDirectory, "/dge_", suffix, ".RDS"))
    return(dge)
}

dap <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    suffix = NULL,
    barcodesToInclude = NULL,
    identToCompare = NULL,
    ident1 = NULL,
    logFCthres = 0,
    minPerc = 0.05,
    latentVars = "atac_peak_region_fragments"
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(Signac)
    require(future)


    if (is.null(seuratObj)) stop("Set the Seurat object to process")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(suffix)) stop("Set suffix for saved objects")
    if (is.null(barcodesToInclude)) stop("Set the cell barcodes to include")

    plan("multicore", workers = 8)
    options(future.globals.maxSize = 5000 * 1024^2)

    # filter out not included barcodes
    seuratObj$barcode <- seuratObj %>% colnames
    seuratObj <- seuratObj %>% subset(barcode %in% barcodesToInclude)
    # set ident
    Idents(seuratObj) <- identToCompare
    # check if there are enough cells in each ident
    cellNums <- table(seuratObj[[identToCompare]])
    if (any(cellNums <= 3) | length(cellNums) == 1) dap <- "Not enough cells" %>% data.frame
    else {
        # do DAP
        dap <- FindMarkers(
            seuratObj,
            ident.1 = ident1,
            test.use = "LR",
            logfc.threshold = logFCthres,
            min.pct = minPerc,
            latent.vars = latentVars
        )
    }
    saveRDS(dap, paste0(outDirectory, "/dap_", suffix, ".RDS"))
    return(dap)
}

dam <- function(
    seuratObj = NULL,
    outDirectory = NULL,
    suffix = NULL,
    barcodesToInclude = NULL,
    identToCompare = NULL,
    ident1 = NULL
) {
    require(dplyr)
    require(stringr)
    require(magrittr)
    require(Seurat)
    require(Signac)

    if (is.null(seuratObj)) stop("Set the Seurat object to process")
    if (is.null(outDirectory)) stop("Set directory for saved objects")
    if (is.null(suffix)) stop("Set suffix for saved objects")
    if (is.null(barcodesToInclude)) stop("Set the cell barcodes to include")

    # filter out not included barcodes
    seuratObj$barcode <- seuratObj %>% colnames
    seuratObj <- seuratObj %>% subset(barcode %in% barcodesToInclude)
    # set ident
    Idents(seuratObj) <- identToCompare
    # check if there are enough cells in each ident
    cellNums <- table(seuratObj[[identToCompare]])
    if (any(cellNums <= 3) | length(cellNums) == 1) dam <- "Not enough cells" %>% data.frame
    else {
        # do DAM
        dam <- FindMarkers(
            seuratObj,
            ident.1 = ident1,
            fc.name = "avg_log2FC",
            mean.fxn = rowMeans
        )
    }
    saveRDS(dam, paste0(outDirectory, "/dam_", suffix, ".RDS"))
    return(dam)
}