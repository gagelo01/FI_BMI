#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
library(gwasvcf)
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
exposures_gtex <- fread("Data/Modified/exposures_gtex_hypothalamus.txt")

exposures_gtex_split <- split(exposures_gtex, exposures_gtex$exposure)
####
from_harm_to_coloc_results <- function(harm) {
  dat_window <- harm
  type1 <- ifelse(dat_window[1,]$units.exposure != "log odds", "quant", "cc")
  type2 <- ifelse(dat_window[1,]$units.outcome != "log odds", "quant", "cc")
  
  D1<- list(exposure = dat_window$exposure[1],
            beta=dat_window$beta.exposure,
            varbeta = dat_window$se.exposure^2,
            N = dat_window$samplesize.exposure,
            MAF = dat_window$eaf.exposure %>% ifelse(. < 0.5, ., 1-.),
            type=type1,
            snp = dat_window$SNP,
            position = dat_window$pos.exposure)
  D2 <- list(outcome = dat_window$outcome[1],
             beta=dat_window$beta.outcome,
             varbeta=dat_window$se.outcome^2,
             N = dat_window$samplesize.outcome,
             type=type2,
             MAF = dat_window$eaf.outcome %>% ifelse(. < 0.5, ., 1-.),
             snp = dat_window$SNP,
             position = dat_window$pos.exposure)
  
  if (type1 == "cc") {
    D1$s <- dat_window[, mean(ncase.exposure/ncontrol.exposure, na.rm = TRUE)] 
  }
  if (type2 == "cc") {
    D2$s <- dat_window[, mean(ncase.outcome/ncontrol.outcome, na.rm = TRUE)] 
  }
  
  vres <- coloc::coloc.abf(dataset1=D1,   dataset2=D2)
  
  vresres <- vres$results %>% as.data.table(.)
  df_res <- data.frame(exposure = D1$exposure[1], outcome = D2$outcome[1], 
                       posprob_coloc.mr = vres$summary[6], posprob_coloc.SNP = vresres[which.max(SNP.PP.H4), 
                                                                                       snp], posprob_coloc.SNPexplained_var = vresres[which.max(SNP.PP.H4), 
                                                                                                                                      SNP.PP.H4])
  return(df_res)
  
}
########
dummy_function <- function(exposure) {
exposure  <- exposure[!(is.na(se.exposure) | is.na(eaf.exposure)),]
tsmr <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", 
                    chrompos = exposure[,paste0(chr.exposure, ":", min(pos.exposure), "-", max(pos.exposure))]) %>%
  gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% 
  as.data.table(.)

harm <- TwoSampleMR::harmonise_data(exposure, tsmr, action = 1)
setDT(harm)
harm[,units.outcome := "SD"]
harm[,units.exposure := "SD"]
res_coloc <- from_harm_to_coloc_results(distinct(harm))
return(res_coloc)
}

dummy_function_safely <- safely(dummy_function)
options(future.globals.maxSize= 1e10)
plan(multisession, workers = 6)
res <-  map(exposures_gtex_split, dummy_function_safely) 

saveRDS(res, "Data/Modified/res_mapping_hypo.rds")
message("This script finished without errors")
