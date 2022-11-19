suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

#' get all sub strings of a seq of certain length
#' @param seq
#' @param seq_id
#' @param kmer_len the length of sub string
#'  -default = 13
get_kmers <- function(seq, seq_id, kmer_len = 13) {
  last_i <- nchar(seq) - kmer_len + 1
  kmers <- vector("character", last_i)
  ids <- vector("character", last_i)
  for (i in 1:last_i) {
    kmer <- str_sub(seq, i, i + kmer_len - 1)
    ids[i] <- str_c(seq_id, ":", i)
    kmers[i] <- kmer
  }
  names(kmers) <- ids
  kmers
}

args <- commandArgs(TRUE)
mibase_file <- args[1]
species <- args[2] # species (e.g. hsa)
kmer_len <- as.integer(args[3])
out_dir <- args[4]

if (!dir.exists(out_dir)) {
  suppressWarnings(dir.create(out_dir))
}

out_mibase_file <- file.path(out_dir, str_c(species, "_mibase.tsv"))
out_mirna_file <- file.path(out_dir, str_c(species, "_mirna.tsv"))
out_pre_file <- file.path(out_dir, str_c(species, "_pre.fa"))
out_sam_file <- file.path(out_dir, str_c(species, "_mirna.sam"))

mibase <- read_tsv(mibase_file, show_col_types = FALSE)
mibase <- select(mibase, Accession:Mature2_Seq)
mibase <- filter(mibase, str_detect(ID, str_c("^", species)))

mibase_5p <- select(mibase,
  pre_id = ID, pre_seq = Sequence,
  mirna_id = Mature1_ID, mirna_seq = Mature1_Seq
)
mibase_3p <- select(mibase,
  pre_id = ID, pre_seq = Sequence,
  mirna_id = Mature2_ID, mirna_seq = Mature2_Seq
)

mibase <- bind_rows(mibase_5p, mibase_3p)
mibase <- filter(mibase, !is.na(mirna_id) | !is.na(mirna_seq))

mibase$pre_seq <- str_replace_all(mibase$pre_seq, "U", "T")
mibase$mirna_seq <- str_replace_all(mibase$mirna_seq, "U", "T")

pos <- str_locate(mibase$pre_seq, mibase$mirna_seq)
mibase <- mutate(mibase, start = pos[, 1], end = pos[, 2])
mibase <- filter(mibase, !is.na(start))
mibase <- unique(mibase)
write_tsv(mibase, file = out_mibase_file)

pre_mibase <- unique(select(mibase, pre_id, pre_seq))
pre <- DNAStringSet(pre_mibase$pre_seq)
names(pre) <- pre_mibase$pre_id
writeXStringSet(pre, format = "fasta", file = out_pre_file)

mature_mibase <- unique(select(mibase, mirna_id, mirna_seq))

kmers <- map2(
  mature_mibase$mirna_seq, mature_mibase$mirna_id, get_kmers,
  kmer_len
)
kmers <- unlist(kmers)

# should check if the length of mirna_seq is shorted than kmer_len?
start <- floor((nchar(mature_mibase$mirna_seq) - kmer_len) / 2) + 1
end <- start + kmer_len - 1
mature_mibase <- mutate(mature_mibase, motif = str_sub(mirna_seq, start, end))

kmer_count_tib <- as.data.frame(table(kmers))
names(kmer_count_tib) <- c("kmer", "count")

mature_mibase <- left_join(mature_mibase, kmer_count_tib,
  by = c("motif" = "kmer")
)
write_tsv(mature_mibase, file = out_mirna_file)

# mibase
# pre_id	pre_seq	mirna_id	mirna_seq	start	end
sam <- tibble(
  QNAME = mibase$mirna_id,
  FLAG = 0,
  RNAME = mibase$pre_id,
  POS = mibase$start,
  MAPQ = 30,
  CIGAR = str_c(nchar(mibase$mirna_seq), "M"),
  RNEXT = "*",
  PNEXT = 0,
  TLEN = 0,
  SEQ = mibase$mirna_seq,
  QUAL = "*"
)

sq <- map2_chr(
  pre_mibase$pre_id, nchar(pre_mibase$pre_seq),
  ~ str_c("@SQ\tSN:", .x, "\tLN:", .y)
)

fout <- file(out_sam_file, "w")
writeLines("@HD\tVN:1.6\tSO:coordinate", fout)
writeLines(sq, fout)
close(fout)

write_tsv(sam, file = out_sam_file, col_names = FALSE, append = TRUE)
