suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
  library(data.table)
  library(dplyr)
  library(R.utils)
})

RAW_DIR <- "data/raw"
META_DIR <- "data/metadata"
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(META_DIR, recursive = TRUE, showWarnings = FALSE)

meta_url <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"
meta_dest <- file.path(RAW_DIR, "TcgaTargetGTEX_phenotype.txt.gz")

if (!file.exists(meta_dest)) {
  message("Downloading metadata...")
  download.file(meta_url, destfile = meta_dest, method = "auto", quiet = FALSE)
}

pheno <- fread(meta_dest, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(pheno) <- make.names(colnames(pheno))

kirc_tumor <- pheno %>%
  filter(X_study == "TCGA",
         grepl("Kidney Clear Cell", detailed_category, ignore.case = TRUE),
         X_sample_type == "Primary Tumor") %>%
  mutate(group = "TCGA_Tumor")

kirc_nat <- pheno %>%
  filter(X_study == "TCGA",
         grepl("Kidney Clear Cell", detailed_category, ignore.case = TRUE),
         X_sample_type == "Solid Tissue Normal") %>%
  mutate(group = "TCGA_NAT")

gtex_kidney <- pheno %>%
  filter(X_study == "GTEX",
         grepl("Kidney", detailed_category, ignore.case = TRUE),
         X_sample_type == "Normal Tissue") %>%
  mutate(group = "GTEx_Healthy")

manifest <- bind_rows(kirc_tumor, kirc_nat, gtex_kidney)
fwrite(manifest, file.path(META_DIR, "kirc_gtex_manifest.tsv"), sep = "\t")

print(table(manifest$group))
