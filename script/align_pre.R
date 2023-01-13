suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)
script_dir <- args[1]
read_file <- args[2]
pre_file <- args[3]
min_identity <- as.numeric(args[4])
output_file <- args[5]

source(file.path(script_dir, "util.R"))

pres <- readDNAStringSet(pre_file, format = "fasta")
reads <- read_tsv(read_file, col_names = FALSE, show_col_types = FALSE)

sub_mat <- nucleotideSubstitutionMatrix(
  match = 1,
  mismatch = 0,
  baseOnly = TRUE
)

align <- function(read) {
  read_id <- read[[1]]
  read_num <- read[[2]]
  seq <- read[[3]]
  seq <- DNAStringSet(seq)
  alns <- pairwiseAlignment(
    pattern = pres,
    subject = seq,
    substitutionMatrix = sub_mat,
    gapOpening = 0,
    gapExtension = 1,
    type = "overlap"
  )
  scores <- score(alns)
  max_indexes <- which(scores == max(scores))
  alns <- alns[max_indexes]
  pre_ids <- names(pres[max_indexes])
  identities <- scores[max_indexes] / width(seq)

  aln_list <- vector("list", length(alns))

  for (i in seq_along(alns)) {
    aln <- alns[i]
    aln_query <- pattern(aln)
    aln_subject <- subject(aln)
    q_start <- start(aln_query)
    q_end <- end(aln_query)
    s_start <- start(aln_subject)
    s_end <- end(aln_subject)

    pad_head <- ifelse(s_start > 1, s_start - 1, 0)
    pad_tail <- ifelse(width(seq) > s_end, width(seq) - s_end, 0)
    str1 <- str_c(
      str_dup("X", pad_head), as.character(aln_subject),
      str_dup("X", pad_tail)
    )
    str2 <- str_c(
      str_dup("Y", pad_head), as.character(aln_query),
      str_dup("Y", pad_tail)
    )
    cigar <- get_cigar(str1, str2)
    start <- q_start - pad_head
    aln_list[[i]] <- tibble(
      read_id = read_id,
      pre_id = pre_ids[i],
      start = start,
      cigar = cigar,
      seq = as.character(seq),
      read_num = read_num,
      identity = identities[i]
    )
  }

  aln_list
}

out_list <- vector("list", nrow(reads))

for (i in seq(nrow(reads))) {
  read <- reads[i, ]
  out_list[[i]] <- align(read)
}

tib <- bind_rows(out_list)

# tib <- map2_dfr(seqs, read_nums, align)

tib <- filter(tib, identity >= min_identity)

write_tsv(tib, file = output_file)
