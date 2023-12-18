suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

mibase_file <- snakemake@input[[1]]
custom_mibase_file <- snakemake@input[[2]]

out_mibase_file <- snakemake@output[[1]]
out_mirna_file <- snakemake@output[[2]]
out_pre_file <- snakemake@output[[3]]
out_sam_file <- snakemake@output[[4]]

spe_name <- snakemake@config[["spe"]]
kmer_len <- as.integer(snakemake@config[["kmer_len"]])

mibase <- read_tsv(mibase_file, show_col_types = FALSE) %>% 
  filter(spe == spe_name) %>% 
  filter(!is.na(mirna_id) | !is.na(mirna_seq)) 

custom_mibase <- read_tsv(custom_mibase_file, show_col_types = FALSE) %>% 
  filter(spe == spe_name) %>% 
  filter(!is.na(mirna_id) | !is.na(mirna_seq))

mibase <- rbind(mibase, custom_mibase) %>%
  select(-spe) %>%
  unique()
write_tsv(mibase, file = out_mibase_file)

pre_mibase <- unique(select(mibase, pre_id, pre_seq))
pre <- DNAStringSet(pre_mibase$pre_seq)
names(pre) <- pre_mibase$pre_id
writeXStringSet(pre, format = "fasta", file = out_pre_file)

mature_mibase <- unique(select(mibase, mirna_id, mirna_seq))

# should check if the length of mirna_seq is shorted than kmer_len?
start <- floor((nchar(mature_mibase$mirna_seq) - kmer_len) / 2) + 1
end <- start + kmer_len - 1
mature_mibase <- mutate(mature_mibase, motif = str_sub(mirna_seq, start, end))
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
