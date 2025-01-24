library(readxl)
library(magrittr)
library(stringr)
library(dplyr)

# read data
data <- list(
    B = data.table::fread("data/BiopsyTXT/Merge_cell_seg_data.txt"),
    S = data.table::fread("data/SurgeryTXT/Merge_cell_seg_data.txt")
)
patDf <- read_xlsx("data/batchAug24/Patients.xlsx", skip = 2) %>%
    select(Slide, `Patient ID`) %>% mutate(Slide = Slide %>% as.character) %>%
    rename(sample = Slide, pat = `Patient ID`)

data_og <- data
# process files
data <- lapply(names(data), function(set) {
    df <- data[[set]]
    # remove NAs
    df <- df[which(df$`Phenotype-KI67+` != ""),]
    # format`
    df %<>% select(
        `Sample Name`, `Tissue Category`,
        `Phenotype-KI67+`, `Phenotype-CD103+`, `Phenotype-CD39+`,
        `Phenotype-CD4+`, `Phenotype-CD8+`, `Phenotype-CK+`) %>%
        rename(
            sample = `Sample Name`, tissue = `Tissue Category`,
            ki67 = `Phenotype-KI67+`, cd103 = `Phenotype-CD103+`,
            cd39 = `Phenotype-CD39+`, cd4 = `Phenotype-CD4+`,
            cd8 = `Phenotype-CD8+`, ck = `Phenotype-CK+`
        ) %>%
        filter(tissue != "NO TISSUE") %>%
        mutate(
            sample = sample %>% str_sub(1, 1),
            tissue = tissue %>% str_to_title %>% factor,
            ki67 = if_else(ki67 == "OTHER+", FALSE, TRUE),
            cd103 = if_else(cd103 == "OTHERS", FALSE, TRUE),
            cd39 = if_else(cd39 == "OTHERS", FALSE, TRUE),
            cd4 = if_else(cd4 == "OTHERS", FALSE, TRUE),
            cd8 = if_else(cd8 == "OTHERS", FALSE, TRUE),
            ck = if_else(ck == "OTHERS", FALSE, TRUE)
        ) %>%
        # add neorhea ids
        left_join(patDf) %>%
        mutate(sample = paste0(pat, "_", set)) %>%
        select(-pat)
}) %>% do.call(rbind, .)

data_splt <- split(data, data$sample)

data_splt$`17010_B` %>% head

dfAll <- data_splt %>% sapply(function(df) {
    return(
        list(
            ki_ck = sum(df$ki67 & df$ck)/sum(df$ck),
            ki_cd4_str = sum(df$ki67 & df$cd4 & df$tissue == "Stroma")/sum(df$cd4 & df$tissue == "Stroma"),
            ki_cd4_tum = sum(df$ki67 & df$cd4 & df$tissue == "Tumor")/sum(df$cd4 & df$tissue == "Tumor"),
            ki_cd8_str = sum(df$ki67 & df$cd8 & df$tissue == "Stroma")/sum(df$cd8 & df$tissue == "Stroma"),
            ki_cd8_tum = sum(df$ki67 & df$cd8 & df$tissue == "Tumor")/sum(df$cd8 & df$tissue == "Tumor"),
            ki_t_str = sum(df$ki67 & (df$cd8 | df$cd4) & df$tissue == "Stroma")/sum((df$cd8 | df$cd4) & df$tissue == "Stroma"),
            ki_t_tum = sum(df$ki67 & (df$cd8 | df$cd4) & df$tissue == "Tumor")/sum((df$cd8 | df$cd4) & df$tissue == "Tumor"),
            cd103_cd4_str = sum(df$cd103 & df$cd4 & df$tissue == "Stroma")/sum(df$cd4 & df$tissue == "Stroma"),
            cd103_cd4_tum = sum(df$cd103 & df$cd4 & df$tissue == "Tumor")/sum(df$cd4 & df$tissue == "Tumor"),
            cd103_cd8_str = sum(df$cd103 & df$cd8 & df$tissue == "Stroma")/sum(df$cd8 & df$tissue == "Stroma"),
            cd103_cd8_tum = sum(df$cd103 & df$cd8 & df$tissue == "Tumor")/sum(df$cd8 & df$tissue == "Tumor"),
            cd103_t_str = sum(df$cd103 & (df$cd4 | df$cd8) & df$tissue == "Stroma")/sum((df$cd4 | df$cd8) & df$tissue == "Stroma"),
            cd103_t_tum = sum(df$cd103 & (df$cd4 | df$cd8) & df$tissue == "Tumor")/sum((df$cd4 | df$cd8) & df$tissue == "Tumor"),
            cd39_cd4_str = sum(df$cd39 & df$cd4 & df$tissue == "Stroma")/sum(df$cd4 & df$tissue == "Stroma"),
            cd39_cd4_tum = sum(df$cd39 & df$cd4 & df$tissue == "Tumor")/sum(df$cd4 & df$tissue == "Tumor"),
            cd39_cd8_str = sum(df$cd39 & df$cd8 & df$tissue == "Stroma")/sum(df$cd8 & df$tissue == "Stroma"),
            cd39_cd8_tum = sum(df$cd39 & df$cd8 & df$tissue == "Tumor")/sum(df$cd8 & df$tissue == "Tumor"),
            cd39_t_str = sum(df$cd39 & (df$cd4 | df$cd8) & df$tissue == "Stroma")/sum((df$cd4 | df$cd8) & df$tissue == "Stroma"),
            cd39_t_tum = sum(df$cd39 & (df$cd4 | df$cd8) & df$tissue == "Tumor")/sum((df$cd4 | df$cd8) & df$tissue == "Tumor")
        )
    )
}) %>% t %>% data.frame %>% apply(2, unlist %>% unname) %>% data.frame %>% tibble::rownames_to_column("sample") %>%
    mutate(
        pat = sample %>% str_sub(1, 5),
        time = sample %>% str_sub(7, 7) %>% recode("B" = "Baseline", "S" = "Surgery")
    )

# add CCCA Yes vs No
clinD <- readRDS("~/Documents/Andreas/Bordet/Projects/neoRhea/shinyApp/data/clinD.RDS")
clinD %<>% select(`Patient ID`, `Response Ki67`)
clinD$`Response Ki67` %<>% recode("Responder" = "Yes", "Non-responder" = "No")
clinD %<>% rename(pat = `Patient ID`, CCCA = `Response Ki67`)
dfAll <- left_join(dfAll, clinD)
saveRDS(dfAll, "informTabs.RDS")
