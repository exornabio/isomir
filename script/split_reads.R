suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)
read_file <- args[1]
split_line_num <- as.integer(args[2])
read_dir <- args[3]

reads <- read_tsv(read_file, col_names = FALSE, show_col_types = FALSE)
reads <- split(reads, ceiling(seq_len(nrow(reads)) / split_line_num))

for (chunk in names(reads)) {
  write_tsv(reads[[chunk]],
    file = str_c(read_dir, "/", chunk, ".tsv"),
    col_names = FALSE
  )
}
