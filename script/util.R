suppressMessages(library(stringr))

get_cigar <- function(query, ref) {
  query <- str_split(query, "")[[1]]
  ref <- str_split(ref, "")[[1]]
  cigar <- rep("M", length(query))
  cigar[query == "-"] <- "D"
  cigar[ref == "-"] <- "I"
  r <- rle(cigar)
  str_c(as.character(r$lengths), r$values, collapse = "")
}
