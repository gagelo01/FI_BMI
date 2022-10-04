#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)#plan(sequential)
#################
inst_map_hypo <-fread( "Data/Modified/inst_map_hypo.txt")
instmapfi <- fread("Data/Modified/instmapfi.txt")

############What you need to change ###########
ID_mrbase_exp <- c("ukb-b-9405", "ukb-b-19953")
exp_mrbase <- paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_exp, "/", ID_mrbase_exp, ".vcf.gz")
ID_server_exp <- c("trait-10-1", "trait-10-3", "trait-1-1", "trait-2-2", paste0("trait-25-", 1:9))
exp_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_exp, "/", ID_server_exp, ".vcf.gz")
arguments <- data.table(id = c(ID_mrbase_exp,ID_server_exp), path = c(exp_mrbase,exp_server), pval = 5e-8, r2 =0.001, kb = 10000)
#################################################
inst <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  GagnonMR::get_inst(x$path, r2 = x$r2, pval = x$pval, kb = x$kb)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)
setDT(inst)
inst[,inst_sel_strat:="statistically_driven"]
inst <- rbindlist(list(inst, inst_map_hypo, instmapfi), fill = TRUE)
inst[inst_sel_strat == "biologically_driven", exposure:=paste0(exposure, "_bio")]
inst[,id.exposure := exposure]

out <- future_map(as.list(c(exp_mrbase,exp_server)), function(x, rsiid = unique(inst$SNP)) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr<-GagnonMR::extract_outcome_variant(snps = rsiid, outcomes = x)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

instmvmr <- TwoSampleMR::convert_outcome_to_exposure(out) %>% as.data.table(.)

fwrite(instmvmr, "Data/Modified/all_inst_mvmr.txt")
fwrite(out, "Data/Modified/all_outcome_mvmr.txt" )
fwrite(inst, "Data/Modified/inst_all_sign_clump.txt")

message("This script finished without errors")
