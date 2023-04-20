.libPaths(c("/usr/local/lib/R/site-library"))

suppressMessages(library(Biostrings))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))

script_dir <- "script"

read_file <- snakemake@input[[1]]
pre_file <- snakemake@input[[2]]
output_file <- snakemake@output[[1]]
min_identity <- as.numeric(snakemake@config[["min_identity"]])

source(file.path(script_dir, "util.R"))

pres <- readDNAStringSet(pre_file, format = "fasta")
reads <- read_tsv(read_file, col_names = FALSE, show_col_types = FALSE)
read_ids <- reads[[1]]
read_nums <- reads[[2]]
seqs <- DNAStringSet(reads[[3]])

sub_mat <- nucleotideSubstitutionMatrix(
  match = 1,
  mismatch = 0
)

align <- function(pre) {
  # pre_id <- names(pre)
  # pre <- DNAString(as.character(pre))
  alns <- pairwiseAlignment(
    pattern = seqs,
    subject = pre,
    substitutionMatrix = sub_mat,
    gapOpening = 0,
    gapExtension = 1,
    type = "overlap"
  )
  identities <- score(alns) / width(seqs)
  hit_indexes <- identities >= min_identity
  alns <- alns[hit_indexes]
  aln_read_ids <- read_ids[hit_indexes]
  aln_read_nums <- read_nums[hit_indexes]
  aln_seqs <- seqs[hit_indexes]

  aln_list <- vector("list", length(alns))

  for (i in seq_along(alns)) {
    aln <- alns[i]
    aln_query <- pattern(aln)
    aln_subject <- subject(aln)
    q_start <- start(aln_query)
    q_end <- end(aln_query)
    s_start <- start(aln_subject)

    pad_head <- ifelse(q_start > 1, q_start - 1, 0)
    pad_tail <- ifelse(width(aln_seqs[i]) > q_end, width(aln_seqs[i]) - q_end,
      0)
    str1 <- str_c(
      str_dup("X", pad_head), as.character(aln_query),
      str_dup("X", pad_tail)
    )
    str2 <- str_c(
      str_dup("Y", pad_head), as.character(aln_subject),
      str_dup("Y", pad_tail)
    )
    cigar <- get_cigar(str1, str2)
    start <- s_start - pad_head
    aln_list[[i]] <- tibble(
      read_id = aln_read_ids[i],
      pre_id = names(pre),
      start = start,
      cigar = cigar,
      seq = as.character(aln_seqs[i]),
      read_num = aln_read_nums[i],
      identity = identities[hit_indexes][i]
    )
  }

  aln_list
}

out_list <- vector("list", length(pres))

for (i in seq_len(length(pres))) {
  pre <- pres[i]
  out_list[[i]] <- align(pre)
}

tib <- bind_rows(out_list)

write_tsv(tib, file = output_file)
