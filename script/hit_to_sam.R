# .libPaths(c("/usr/local/lib/R/site-library"))

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(Biostrings))

script_dir <- "script"

hit_file <- snakemake@input[[1]]
read_file <- snakemake@input[[2]]
pre_file <- snakemake@input[[3]]
sam_file <- snakemake@output[[1]]

# read_id   pre_id    start   cigar   seq   read_num    identity
hits <- read_tsv(hit_file, col_types = "cciccid") %>%
    arrange(pre_id, desc(read_num))

reads <- read_tsv(read_file, col_names = FALSE, show_col_types = FALSE)
names(reads) <- c("id", "num", "seq")
total_read_size <- sum(nchar(reads$seq) * reads$num) / 1000000

hits <- mutate(hits,
    tpm = read_num / (nchar(seq) / 1000) / total_read_size)

pres <- readDNAStringSet(pre_file, format = "fasta")
pre_lens <- width(pres)
pre_ids <- names(pres)
sq <- map2_chr(pre_ids, pre_lens, ~ str_c("@SQ\tSN:", .x, "\tLN:", .y))

fout <- file(sam_file, "w")
writeLines("@HD\tVN:1.6\tSO:coordinate", fout)
writeLines(sq, fout)
close(fout)

sam <- tibble(
    QNAME = str_c("read", hits$read_id),
    FLAG = 0,
    RNAME = hits$pre_id,
    POS = hits$start,
    MAPQ = 30,
    CIGAR = hits$cigar,
    RNEXT = "*",
    PNEXT = 0,
    TLEN = 0,
    SEQ = hits$seq,
    QUAL = "*",
    READ_NUM = str_c("RN:i:", hits$read_num),
    TPM = str_c("TM:f:", round(hits$tpm, 2)),
    ID = str_c("ID:i:", round(hits$identity, 3))
)

write_tsv(sam, file = sam_file, col_names = FALSE, append = TRUE)
