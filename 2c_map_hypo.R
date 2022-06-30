#!/usr/bin/env Rscript
library("xlsx")
library(data.table)
library(data.table)
library(tidyverse)
library(GagnonMR)
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()

res <- readRDS("Data/Modified/res_mapping_hypo.rds")
res<-unname(res)
index <- sapply(res, function(x) !is.null(x$error))
res <- rbindlist(lapply(res[!index], function(x)  x$result), fill = TRUE)
res <- res[posprob_coloc.mr > 0.8,]

ldmat <- ieugwasr::ld_matrix_local(res$posprob_coloc.SNP %>% unique, plink_bin = genetics.binaRies::get_plink_binary(), 
                                   bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")
diag(ldmat)<- 0
max(ldmat^2) #yes
ldmat

shortharm<- function(inst, vcffile_out,ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") {
  out_vcf <- gwasvcf::query_gwas(vcf = vcffile_out,
                                 rsid = inst$SNP, proxies = "yes", bfile = ldref, tag_kb = 10000, tag_r2 = 0.8)
  
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)
  
  harm <- TwoSampleMR::harmonise_data(inst, out_tsmr, action = 1)
  return(harm) 
}


dat_vcf <- gwasvcf::query_gwas("/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", 
                               rsid = res[!is.na(exposure), posprob_coloc.SNP],
                               proxies = "no")
inst_map <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(. , "exposure") %>% data.table::as.data.table(.)
inst_map <- GagnonMR::clump_data_local(inst_map)
harm <- shortharm(inst = inst_map,
                    vcffile_out = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz") 
  
harm <- TwoSampleMR::add_rsq(harm)
harm <- TwoSampleMR::steiger_filtering(harm)
setDT(harm)
resmap <- harm %>% TwoSampleMR::mr(., method = "mr_ivw")
resmap <- harm %>% GagnonMR::all_mr_methods(.)

harmsteiger <- harm[!(steiger_dir == FALSE & steiger_pval < 0.05), ]
resmapsteiger <- harmsteiger %>% GagnonMR::all_mr_methods(.)

fwrite(res, "Data/Modified/eQTLcoloc_hypo.txt")
fwrite(harm, "Data/Modified/harm_hypo.txt")
fwrite(resmapsteiger,  "Data/Modified/resmapsteiger_hypo.txt")

