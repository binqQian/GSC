#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: gsc2gds.R <input.vcf[.gz]> <output.gds>\n", file = stderr())
  quit(status = 2)
}

vcf_file <- args[[1]]
gds_file <- args[[2]]

suppressPackageStartupMessages({
  library(SeqArray)
})

seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT")
