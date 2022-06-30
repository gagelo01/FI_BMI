#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()


df_index_eqtl <- df_index[pmid == 32999275, ]
df_index_eqtl <- df_index_eqtl[,.(id, note)]
df_index_eqtl[, ID_exposure_file := paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id, "/", id, ".vcf.gz")]
df_index_eqtl[,id:=NULL]
setnames(df_index_eqtl, "note", "chrompos")
id_to_include <- c("trait-2-2", "trait-17-1")
dt_vcf <- data.table(ID_outcome_file = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id_to_include, "/", id_to_include,".vcf.gz"))
df_index_eqtl <- tidyr::crossing(df_index_eqtl, dt_vcf)
list_arg <- split(df_index_eqtl, seq(nrow(df_index_eqtl)))

run_mapping <- function(comb, method_list = list("get_coloc")) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  res_all<-  GagnonMR::run_all_pqtl_analyses(comb$ID_exposure_file, comb$ID_outcome_file, chrompos = comb$chrompos, method_list = method_list)
  return(res_all)
}

options(future.globals.maxSize= 1e10)
plan(multisession, workers = 20)

run_mapping_safely <- safely(run_mapping)
res <- future_map(list_arg, function(x) {run_mapping_safely(x, method_list = list("get_coloc"))},
                  .options = furrr_options(seed = TRUE)) 

saveRDS(res, "Data/Modified/vinuela_res_mapping")
# fwrite(res, "Data/Modified/res_mapping.txt")
message("This script finished without errors")
