library(stringr)

# output of: ll sc_RNA_ATAC/fastqs/*/*
fileContent <- readLines("~/Downloads/fastqFiles.txt")

# get indexes of lines containing directories
dirIdxs <- grep(":$", fileContent)
# extract directory names from the fileContent
dirNames <- fileContent[dirIdxs] %>% str_remove(":")
# add prefix where ll was invoced
dirNames %<>% paste0("/scratch/ulb/BCTL/apapagia/", .)
# add end of file to dir indeces
dirIdxs <- c(dirIdxs, length(fileContent))

# Initialize an empty list to store the file_list
files <- sapply(dirNames, function(dir) {
    # find which of the directories is it
    dirIdx <- which(dir == dirNames)
    # get ranges
    ranges <- dirIdxs[dirIdx]:dirIdxs[dirIdx+1 ]
    # get all files of that dir
    files <- fileContent[ranges]
    # remove non fastq file entries
    files <- files[-grep(":$|^total|.html$", files)]
    files <- files[!(files == "")]
    # get last "word", the filename without ls details
    files %<>% sapply(word, -1) %>% unname
    # add the parent directory
    files %<>% paste0(dir, "/", .)
    return(files)
})

files %<>% unlist %>% unname

# assign a sample name to each of the files
names(files) <- files %>%
    basename %>%
    str_remove(".fastq.gz") %>%
    # wizzardry 
    str_remove_all("-") %>%
    gsub("(\\d{4,})", "___\\1___", .) %>%
    str_split("___") %>%
    sapply(function(vec) paste0(vec[2], vec[3] %>% str_sub(1, 2))) %>%
    # fix typos
    str_replace("0025012SG", "25012SG") %>%
    str_replace("150015SG", "15001SG") %>%
    str_replace("170016", "17016") %>%
    str_replace("170032", "17032") %>%
    str_replace("210019SG", "25019SG") %>%
    str_replace("250019BL", "25019BL") %>%
    str_replace("250002SG", "25002SG") %>%
    str_replace("250015", "25015") %>%
    str_replace("250028", "25028") %>%
    str_replace("250031", "25031") %>%
    str_replace("3009SG", "30009SG")

# split files by sample
files <- split(files %>% unname, names(files))

# now make a df for each sample with 3 columns: fastqs,sample,library_type
# as in https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count
filesDf <- lapply(files %>% names %>% setNames(., .), function(sample) {
    filesOfSample <- files[[sample]]
    df <- data.frame(
        fastqs = filesOfSample %>% dirname,
        sample = filesOfSample %>% basename %>% str_split("_") %>% sapply(`[`, 1)
    ) %>% unique
    df$library_type = ifelse(grepl("GEX", df$fastqs), "Gene Expression", "Chromatin Accessibility")
    return(df)
})

# write in csvs
lapply(names(filesDf), function(sample) {
    df <- filesDf[[sample]]
    write.csv(df, paste0("cellrangerLibraryCSVs/", sample, ".csv"), row.names = FALSE)
})
 