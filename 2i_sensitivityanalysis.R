#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index<- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao<-fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/")]

harm <- fread("Data/Modified/harm.txt")
harm <- harm[SNP!= "rs1121980",] #FTO SNP
###



harm[, exposure %>% unique]
harm <- harm[exposure == "Fasting_Insulin" & outcome == "bmi_ukbgiant", ] 
inst_fi <- harm[,.SD, .SDcols = c("SNP",colnames(harm)[grepl("exposure",colnames(harm))])]
harm[,id.exposure := exposure];harm[,id.outcome := outcome]
thesnp <- harm[exposure == "Fasting_Insulin" & outcome == "bmi_ukbgiant", SNP] 
id <- c("ieu-b-109", "ieu-b-111")
idpath <- c(paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/",id, "/", id, ".vcf.gz"),
            "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-15-5/trait-15-5.vcf.gz")
tsmr_list <- vector(mode = "list", length = length(id))
for(i in 1:length(idpath)) {
tsmr_list[[i]] <- gwasvcf::query_gwas(idpath[i], rsid = thesnp,
                                      proxies = "yes", bfile = ldref, tag_r2 = 0.8) %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
}
tsmr <- rbindlist(tsmr_list, fill = TRUE)
tsmr[,id.outcome:=outcome]
####il faut que j'haromise
col_toselect <- colnames(harm)[!(grepl("outcome", colnames(harm)))]
inst <- harm[outcome == "bmi_ukbgiant",.SD, .SDcols = col_toselect]
inst[, action := NULL]
k <- TwoSampleMR::harmonise_data(exposure = inst, outcome = tsmr, action = 1) %>% as.data.table(.)
snphdl <- k[outcome == "ieu-b-109" & ((beta.outcome*beta.exposure)<0) & pval.outcome < 5e-8, SNP]
snptri <- k[exposure == "ieu-b-111" & ((beta.outcome*beta.exposure)>0) & pval.exposure < 5e-8, SNP]
snp_ir<-unique(snphdl, snptri) #These SNPs show association with high triglycerides or low HDL, hallmarks of common insulin resistance
snp_is_prokopenko <- k[outcome == "Ins30" & pval.outcome < 0.05, ]$SNP

harm_is <- harm[!(SNP %in% snp_ir), ]
harm_is_split <- split(harm_is, harm_is$outcome)
res_is <- map(harm_is_split, function(x) GagnonMR::all_mr_methods(x)) %>%
  rbindlist(.,fill = TRUE)

harm_ir <- harm[(SNP %in% snp_ir), ]
harm_ir_split <- split(harm_ir, harm_ir$outcome)
res_ir <- map(harm_ir_split, function(x) GagnonMR::all_mr_methods(x)) %>%
  rbindlist(.,fill = TRUE)

#the prokopenko GWAS is not enough precise so I remove the sensitivity analysis.
res_ir<-res_ir[outcome!="Ins30"]
harm_ir<-harm_ir[outcome!="Ins30"]
res_is<-res_is[outcome!="Ins30"]
harm_is<-harm_is[outcome!="Ins30"]
fwrite(harm_is,"Data/Modified/harm_is.txt")
fwrite(res_is, "Data/Modified/res_is.txt")

fwrite(harm_ir, "Data/Modified/harm_ir.txt")
fwrite(res_ir, "Data/Modified/res_ir.txt")
############
dtcir<- data.table(SNP = c("rs933360", "rs10830963", "rs7923866", "rs7756992", "rs11671664", "rs4502156", "rs3757840", "rs12549902"),
           effect_allele.outcome = c("A", "G", "C", "G", "A", "T", "T", "A"),
           other_allele.outcome = c("G", "C", "T", "A", "G", "C", "G", "G"),
           eaf.outcome = c(0.62, 0.69, 0.62, 0.3, 0.11, 0.52, 0.45, 0.58),
           beta.outcome = c(-0.051, -0.17, -0.12, -0.11, -0.17, -0.092, -0.090, -0.60),
           se = c(0.0086, 0.016, 0.015, 0.015, 0.025, 0.014, 0.015, 0.010))
dtcir$SNP

