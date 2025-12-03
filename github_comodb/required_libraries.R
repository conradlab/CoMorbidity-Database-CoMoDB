#"R version 4.2.1 (2022-06-23 ucrt)"
{r}
# CRAN packages
cran_pkgs <- c(
  "httr",
  "jsonlite",
  "gprofiler2"
)

# Bioconductor packages
bioc_pkgs <- c(
  "GenomeInfoDb",
  "GenomicRanges",
  "rtracklayer",
  "DiffBind"
)

# Install CRAN packages if missing
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg, " ...")
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all packages
all_pkgs <- c(cran_pkgs)
for (pkg in all_pkgs) {
  library(pkg, character.only = TRUE)
}