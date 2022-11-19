
# command_line_programs <- c("", "")
# programs_present <- sapply(command_line_programs, function(app) {
#   system2("which", app, stderr = NULL, stdout = NULL)
# }) == 0

# if (any(!programs_present)) {
#   stop(paste(command_line_programs[!programs_present]), " is not available")
# }

# R packages
r_packages <-
  c(
    "ShortRead",
    "Biostrings",
    "tidyverse"
  )

r_packages_present <-
  is.element(r_packages, installed.packages()[, 1])
if (any(!r_packages_present)) {
  stop(paste(r_packages[!r_packages_present]), " is not available")
}
