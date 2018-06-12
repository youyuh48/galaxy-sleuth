#!/usr/bin/env Rscript

# A command-line interface to sleuth for use with Galaxy This script modified
# from https://github.com/pachterlab/bears_analyses/blob/master/sleuth.R
# https://github.com/nturaga/bioc-galaxy-integration/blob/master/README.md

## Command to run tool:
## Rscript sleuth.R --indir test-rscript --metadata test-rscript/metadata.txt
## --full_model '~condition' --reduced_model '~1'
## --gene_anno_name 'hsapiens_gene_ensembl'

# setup R error handling to go to stderr
options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
library("tools")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
spec <- matrix(c(
  "quiet",          "q", 0, "logical",
  "help",           "h", 0, "logical",
  "indir",          "i", 1, "character",
  "metadata",       "m", 1, "character",
  "full_model",     "f", 1, "character",
  "reduced_model",  "r", 1, "character",
  "gene_anno_name", "a", 2, "character"),
  byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message and exit with a non-zero error
# code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

# enforce the following required arguments
if (is.null(opt$indir)) {
  cat("'indir' is required\n")
  q(status = 1)
}
if (is.null(opt$metadata)) {
  cat("'metadata' is required\n")
  q(status = 1)
}
if (is.null(opt$full_model)) {
  cat("'full_model' is required\n")
  q(status = 1)
}
if (is.null(opt$reduced_model)) {
  cat("'reduced_model' is required\n")
  q(status = 1)
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
} else {
  FALSE
}

suppressPackageStartupMessages({
  library("sleuth")
  library("biomaRt")
})

s2c <- read.table(file.path(opt$metadata), header = TRUE, stringsAsFactors = FALSE)
run_dirs <- s2c$sample
kal_dirs <- c()

for (dir in run_dirs) {
  kal_dirs <- c(kal_dirs, file.path(opt$indir, dir, "kallisto"))
}

s2c <- dplyr::mutate(s2c, path = kal_dirs)

if (!is.null(opt$gene_anno_name)) {
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = opt$gene_anno_name)
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id,
    ext_gene = external_gene_name)
  so <- sleuth_prep(s2c, as.formula(opt$full_model), target_mapping = t2g, read_bootstrap_tpm = TRUE,
    extra_bootstrap_summary = TRUE)
} else {
  so <- sleuth_prep(s2c, as.formula(opt$full_model), read_bootstrap_tpm = TRUE,
    extra_bootstrap_summary = TRUE)
}
so <- sleuth_fit(so, as.formula(opt$full_model), "full")
so <- sleuth_fit(so, as.formula(opt$reduced_model), "reduced")
so <- sleuth_lrt(so, "reduced", "full")
sleuth_deploy(so, opt$indir)

cat("Successfully finished script.\n")
