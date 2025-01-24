library(data.table)
library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)
library(tibble)

# convert to dataframe by sample
dataFiles <- data.frame(file = list.files("qtiffs/Results_qupath_20240715", pattern = ".txt$", full.names = TRUE)) %>%
    mutate(
        sample = file %>% basename %>% str_sub(1, 6),
        patient = sample %>% str_sub(1, 5),
        timepoint = sample %>% str_sub(6, 6) %>% recode("B" = "Baseline", "S" = "Surgery") 
    )
# filter out bad staining samples
dataFiles %<>% filter(!sample %in% c("17024S", "25032S")) 

# for each sample make a df
tabs <- lapply(dataFiles$sample %>% setNames(., .), function(sample) {
    message(sample)
    df <- fread(dataFiles$file[dataFiles$sample == sample]) %>%
        select(Classification) %>%
        mutate(
            CK = if_else(grepl("CK", Classification), TRUE, FALSE),
            CD45 = if_else(grepl("CD45", Classification), TRUE, FALSE),
            Ki67 = if_else(grepl("Ki67", Classification), TRUE, FALSE),
            PDL1 = if_else(grepl("PDL1", Classification), TRUE, FALSE),
            PDL2 = if_else(grepl("PDL2", Classification), TRUE, FALSE)
        ) %>% drop_na()
    df$sample <- sample
    return(df)
})
tabs <- do.call(rbind, tabs)
tabs$sample %<>% factor
saveRDS(tabs, "qtiffs/qupathTabs.RDS")