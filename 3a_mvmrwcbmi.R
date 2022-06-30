#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(ggrepel)
library(mrclust)
library(TwoSampleMR)
library("xlsx")

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

#wc_giant = ieu-a-79 / wc_ukb = ukb-b-9405
vcffile_exp <- c("/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", "/mnt/sda/gagelo01/Vcffile/MRBase_vcf/ukb-b-9405/ukb-b-9405.vcf.gz")


inst<- lapply(as.list(vcffile_exp), function(x) GagnonMR::get_inst(x)) %>% rbindlist(., fill = TRUE)


d1 <- lapply(as.list(vcffile_exp), function(x){ ok <- gwasvcf::query_gwas(vcf = x, rsid = inst$SNP) %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(.,  "exposure") %>% data.table::as.data.table(.)
return(ok)}) %>% rbindlist(., fill = TRUE)


inst[, id.exposure := exposure]
d1[, id.exposure := exposure]
inst_mvmr <- prepare_for_mvmr(exposure_dat = inst, d1 = d1, harmonise_strictness = 1, clump_exp = NULL)


outcome_mvmr <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz", rsid = inst$SNP) %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(.,  "outcome") %>% data.table::as.data.table(.)

exposure_outcome_harmonized <- mv_harmonise_data(exposure_dat = inst_mvmr,
                                                 outcome_dat = outcome_mvmr,
                                                 harmonise_strictness = 1)



mvmr_results <- mv_multiple_MendelianRandomization(exposure_outcome_harmonized)
# mvmr_results <- mv_multiple(exposure_outcome_harmonized)$result
setDT(mvmr_results)
# mvmr_results<-mvmr_results[,.(method, exposure, outcome, b,se,CILower,CIUpper,pval, nsnp)]
# mvmr_results %>% setnames(., colnames(.), c("method", "Exposure", "Outcome",  "Estimate", "StdError", "CILower", "CIUpper",  "Pvalue", "nsnp"))
mvmr_results

#univariate
shortres<- function(inst, vcffile_out) {
  out_vcf <- gwasvcf::query_gwas(vcf = vcffile_out,
                                 rsid = inst$SNP, proxies = "yes", bfile = ldref)
  
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)
  
  harm <- TwoSampleMR::harmonise_data(inst, out_tsmr, action = 1)
  harm <- TwoSampleMR::steiger_filtering(harm) %>% as.data.table()
  harm[, .(rsq.exposure, rsq.outcome, steiger_dir, steiger_pval)]
  harm <- harm[steiger_dir == TRUE,]
  res <- harm %>% GagnonMR::all_mr_methods(., short = FALSE)
  return(res)
}

univariate_results <- lapply(as.list(inst$exposure %>% unique), function(x) {
res <- shortres(inst[exposure == x,], "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz")
return(res)}) %>% rbindlist(., fill = TRUE)


# univariate_results <- lapply(as.list(vcffile_exp), function(x) GagnonMR::get_pan(x, 
#                                                  "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz")) %>%
#   rbindlist(., fill = TRUE)
# colnames(univariate_results)<-c("outcome", "exposure", "b", "se", "pval", "nsnp")
# univariate_results[, lci := b-1.96*se]
# univariate_results[, uci := b+1.96*se]
# univariate_results[, method := "univariate IVW"]
mvmr_results <- rbind(univariate_results, mvmr_results, fill = TRUE)
mvmr_results[, exposure := sub("ieu-a-79", "wc_giant", exposure) %>% gsub("UKB-b-9405", "wc_ukb", . )]

fwrite(mvmr_results, "Results/mvmrbmiwcfi.txt")
