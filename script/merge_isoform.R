suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

args <- commandArgs(TRUE)

script_dir <- args[1]
mibase_file <- args[2]
pre_file <- args[3]
isoform_files <- args[4]
read_file <- args[5]
out_file <- args[6]

source(file.path(script_dir, "util.R"))

# pre_id	pre_seq	mirna_id	mirna_seq	start	end
mibase <- read_tsv(mibase_file, show_col_types = FALSE)

isoform_files <- str_split(isoform_files, ",")[[1]]
# mirna_id  read_id  read_seq    read_num    dist
isoforms <- map_dfr(isoform_files, read_tsv,
    col_types = "cccii") %>%
    left_join(mibase, by = "mirna_id") %>%
    arrange(pre_id, desc(read_num))

queries <- DNAStringSet(isoforms$read_seq)
subjects <- DNAStringSet(isoforms$mirna_seq)

query_num <- length(queries)
identities <- vector("double", query_num)
cigars <- vector("character", query_num)
query_indexes <- vector("character", query_num)
offsets <- vector("integer", query_num)

sub_mat <- nucleotideSubstitutionMatrix(
    match = 1,
    mismatch = 0
)

for (i in seq(query_num)) {
    query <- queries[i]
    subject <- subjects[i]
    alns <- pairwiseAlignment(
        pattern = query,
        subject = subject,
        substitutionMatrix = sub_mat,
        gapOpening = 0,
        gapExtension = 1,
        type = "global"
    )
    aln <- alns[1]
    aln_query <- pattern(aln)
    aln_subject <- subject(aln)
    q_start <- start(aln_query)
    q_end <- end(aln_query)
    s_start <- start(aln_subject)
    s_end <- end(aln_subject)
    offsets[i] <- q_start - s_start
    pad_head <- ifelse(q_start > 1, q_start - 1, 0)
    pad_tail <- ifelse(width(query) > q_end, width(query) - q_end, 0)
    str1 <- str_c(
        str_dup("X", pad_head), as.character(aln_query),
        str_dup("X", pad_tail)
    )
    str2 <- str_c(
        str_dup("Y", pad_head), as.character(aln_subject),
        str_dup("Y", pad_tail)
    )
    cigars[i] <- get_cigar(str1, str2)
    query_indexes[i] <- i
    identities[i] <- score(aln) / width(subject)
}


reads <- read_tsv(read_file, col_names = FALSE, show_col_types = FALSE)
names(reads) <- c("id", "num", "seq")
total_read_size <- sum(nchar(reads$seq) * reads$num) / 1000000

isoforms <- mutate(isoforms,
    tpm = read_num / (nchar(read_seq) / 1000) / total_read_size
)

sam <- tibble(
    QNAME = str_c("read", isoforms$read_id),
    FLAG = 0,
    RNAME = isoforms$pre_id,
    POS = isoforms$start - offsets,
    MAPQ = 30,
    CIGAR = cigars,
    RNEXT = "*",
    PNEXT = 0,
    TLEN = 0,
    SEQ = isoforms$read_seq,
    QUAL = "*",
    READ_NUM = str_c("RN:i:", isoforms$read_num),
    TPM = str_c("TM:f:", round(isoforms$tpm, 2)),
    DIST = str_c("DT:i:", isoforms$dist)
)

pres <- readDNAStringSet(pre_file, format = "fasta")
pre_lens <- width(pres)
pre_ids <- names(pres)
sq <- map2_chr(pre_ids, pre_lens, ~ str_c("@SQ\tSN:", .x, "\tLN:", .y))

fout <- file(out_file, "w")
writeLines("@HD\tVN:1.6\tSO:coordinate", fout)
writeLines(sq, fout)
close(fout)

write_tsv(sam, file = out_file, col_names = FALSE, append = TRUE)
