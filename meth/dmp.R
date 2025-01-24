library(magrittr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(edgeR)
library(minfi)
library(ELMER)

dir.create("plots")
dir.create("Robj")

rm(list = ls())

# get methylation annotations ----
# first FFPE methylation annotations
patientsAnnot <- readxl::read_xlsx("data/Echantillons_DIAGENODE_METHYLATION_Second_Shipment_ANnotation.xlsx", skip = 1)
# load FFPE annotation file and add basename location
baseDir <- "data/009_0392/raw_data"
targets_FFPE <- read.csv(file.path(baseDir, "1_SampleAnnotation.txt"), stringsAsFactors = F, sep = "\t")
targets_FFPE %<>% arrange(ChipID)
# add basenames in targets_FFPE dataframe
basenames = dir(list.dirs(baseDir)[2:9], pattern = ".idat", recursive = T, full.names = T)
basenames %<>% str_sub(end = -10) %>% unique %>% sort
targets_FFPE$Basename <- basenames; rm(basenames, baseDir)
# subset including only patient samples
targets_FFPE <- targets_FFPE[grepl("NeoDia", targets_FFPE$SampleID),]
targets_FFPE %<>% arrange(SampleID)
# add sample type to targets_FFPE file
patientsAnnot <- reshape2::melt(patientsAnnot, id.vars = c("...1"))
colnames(patientsAnnot) <- c("Patient_ID", "Timepoint", "Sample_ID")
patientsAnnot %<>% arrange(Sample_ID)
targets_FFPE$Timepoint <- ifelse(patientsAnnot$Timepoint == "Biopsy", "Baseline", "Surgery")
targets_FFPE$Patient <- patientsAnnot$Patient_ID; rm(patientsAnnot)
# load OCT methylation annotations
baseDir = "data/009_0278/raw_data"
targets_OCT = read.csv(file.path(baseDir, "sample_annotation.csv"), stringsAsFactors = F)
targets_OCT$Basename = file.path(baseDir, targets_OCT$Sentrix_ID, paste0(targets_OCT$Sentrix_ID, "_", targets_OCT$Sentrix_Position))
targets_OCT %<>% select(ChipID, neoDiaID, Basename, Sample_type, Patient) 
names(targets_OCT) <- names(targets_FFPE)
# add wet lab info and concatenate targets in a dataframe
targets_FFPE$Tissue <- "FFPE"
targets_OCT$Tissue <- "OCT"
annot_meth <- rbind(targets_OCT, targets_FFPE); rm(baseDir, targets_FFPE, targets_OCT)

# load meth data & do quality control (add pval in coldata)----
rgSet <- read.metharray.exp(targets = annot_meth)
# rename samples
colnames(rgSet) <- rgSet$SampleID
# calculate the detection p-values
detP <- detectionP(rgSet)
# add colMean pval as variable in colData
rgSet$detPval <- colMeans(detP); rm(detP)
# minfi qc report
qcReport(rgSet, sampNames = rgSet$SampleID, sampGroups = rgSet$Tissue, pdf = "plots/qcReportMinfi_tissue.pdf")
qcReport(rgSet, sampNames = rgSet$SampleID, sampGroups = rgSet$Timepoint, pdf = "plots/qcReportMinfi_timepoint.pdf")
# remove bad samples
badSamples <- c("NeoDia_91", "NeoDia_96")
goodSamples <- setdiff(rgSet$SampleID, badSamples)
rgSet <- rgSet[, rgSet$SampleID %in% goodSamples]
# convert RGChannel set to genomic metyl and ratio sets mapping probes to genome
mSet <- preprocessFunnorm(rgSet) # functional normalization
# filter probes that measure SNPs (RS) or non-CpG methylation ones (CH)
mSet <- dropMethylationLoci(mSet, dropRS = T, dropCH = T)
# filter reactive probes
keep <- !(featureNames(mSet) %in% unlist(maxprobes::xreactive_probes(array_type = "EPIC")))
mSet <- mSet[keep,]; rm(keep)

# checkpoint!
saveRDS(rgSet, "Robj/rgSet.RDS")
saveRDS(mSet, "Robj/mSet.RDS")
saveRDS(mSet@colData, "Robj/mSet_meta.RDS")

# plot probes pie by annotations
ann <- getAnnotation(mSet)

# DMP
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

dmps <- getResTable(
    betaMeth = mSet@assays@data$Beta,
    annot = mSet@colData,
    individual = "Patient",
    factorOfInterest = "Timepoint",
    contrast = "Surgery-Baseline"
    )

saveRDS(dmps, "Robj/dmps.RDS")

# split responders non responders
mSet$patID <- mSet$Patient %>% str_sub(-5, -1)
clinD <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/clinicoPathological/data/clinDf.RDS")
clinD %<>% dplyr::rename(`Ultrasound response` =  `DR Result`)
clinD$`Ultrasound response` <- ifelse(clinD$`Ultrasound response` %in% c("SD", "PD"), "No", "Yes")
clinD$`Ultrasound response` %>% table
clinD <- clinD[match(mSet$patID, clinD$`Patient ID`),]
all.equal(clinD$`Patient ID`, mSet$patID)
mSet$`Ultrasound response` <- clinD$`Ultrasound response`

mSetSplit <- lapply(mSet$`Ultrasound response` %>% unique %>% setNames(., .), function(us) {
    mSet[, mSet$`Ultrasound response` == us]
})

dmps <- lapply(mSetSplit, function(mS) {
    getResTable(
        betaMeth = mS@assays@data$Beta,
        annot = mS@colData,
        individual = "Patient",
        factorOfInterest = "Timepoint",
        contrast = "Surgery-Baseline"
        )
})

saveRDS(dmps, "Robj/dmps_DR.RDS")

# get distal probes
distal.probes <- ELMER::get.feature.probe(genome = "hg19", met.platform = "EPIC", promoter = F)
distal.probes %>% names %>% head

saveRDS(distal.probes %>% names, "Robj/distal_probes.RDS")

# only baseline, responders vs resistors
mSet_f <- mSet[,mSet$Timepoint == "Baseline"]

dmps <- getResTable(
    betaMeth = mSet_f@assays@data$Beta,
    annot = mSet_f@colData,
    individual = "Patient",
    factorOfInterest = "Ultrasound response",
    contrast = "Yes-No"
    )

saveRDS(dmps, "Robj/dmps_DR_BL.RDS")
