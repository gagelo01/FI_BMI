#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index<- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

options(future.globals.maxSize= 5e9)
plan(multisession, workers = 6, gc = TRUE) #I should try using multicore
forward_exposure<-df_index[pmid == 24699409, trait]
forward_outcome <- df_index[id == "trait-1-1", trait]

get_inst_safely <- safely(GagnonMR::get_inst)
arguments <- rbind(data.table(id = df_index[pmid == 24699409,]$id, pval = 5e-06), 
                              data.table(id ="trait-1-1", pval = 5e-8))

inst <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  instr <- get_inst(paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", x$id, "/", x$id, ".vcf.gz"), pval = x$pval)
  return(instr)
}, .options = furrr_options(seed = TRUE)) %>% 
  rbindlist(.,fill=TRUE)


out <- future_map(arguments$id, function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr <- gwasvcf::query_gwas(vcf = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", x, "/", x, ".vcf.gz"), 
                              rsid = inst[,unique(SNP)], proxies = "yes", bfile = ldref)
  outr <- outr %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


arguments_forward <- cbind(data.table(exposure = forward_exposure),  data.table(outcome = forward_outcome))
arguments_reverse <- cbind(data.table(exposure = forward_outcome),  data.table(outcome = forward_exposure))
argument <- rbind(arguments_forward, arguments_reverse)
argument_split <- split(argument, 1:nrow(argument))

list_harm_univariate <- map(argument_split, function(x)
TwoSampleMR::harmonise_data(inst[exposure == x$exposure,], out[outcome == x$outcome], action = 1))
lapply(list_harm_univariate, setDT)
lapply(list_harm_univariate, function(x) x[,id.exposure:=exposure])
lapply(list_harm_univariate, function(x) x[,id.outcome:=outcome])

all_mr_methods_safely <- safely(GagnonMR::all_mr_methods)
res_forward <- map(list_harm_univariate[1:9], function(x) {all_mr_methods(x)}) %>%
  rbindlist(., fill = TRUE)
res_reverse <- future_map(list_harm_univariate[10:18], function(x) {all_mr_methods(x)}, .options = furrr_options(seed = TRUE)) %>%
  rbindlist(., fill = TRUE)

res <- rbind(res_forward, res_reverse, fill = TRUE)

fwrite(res, "Data/Modified/res_propokenko.txt")
