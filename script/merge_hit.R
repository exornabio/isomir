suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

args <- commandArgs(TRUE)

pre_file <- args[1]
hit_files <- args[2]
read_file <- args[3]
out_file <- args[4]

hit_files <- str_split(hit_files, ",")[[1]]

# read_id   pre_id    start   cigar   seq   read_num    identity
hits <- map_dfr(hit_files, read_tsv, col_types = "cciccid") %>%
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

fout <- file(out_file, "w")
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

write_tsv(sam, file = out_file, col_names = FALSE, append = TRUE)
