#!/usr/bin/env Rscript
library(biomaRt)
library(data.table)
library(tidyverse)
library(GagnonMR)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao<-fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
datatable <- rbind(ao[id == "ukb-b-9405", ],
      df_index[id %in% c("trait-10-1", "trait-1-1", "trait-2-2","eqtl-2-1", "trait-25-1"), ], fill = TRUE)

datatable[, c("id","build","category","subcategory","ontology","mr","priority","sd","ncase","ncontrol") := NULL]

eQTLcoloc <- fread( "Data/Modified/eQTLcoloc.txt")
eQTLcoloc <- eQTLcoloc[outcome != "Fasting_Insulin_correctedBMI", ]
eQTLcoloc <- eQTLcoloc[order(-posprob_coloc.mr)]
mapharm <- fread("Data/Modified/mapharm.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao<-fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
res_steiger <- fread( "Data/Modified/res_steiger.txt")
mvmr_results <- fread("Results/resmvmr.txt")
FandQ <- fread( "Data/Modified/FandQ_univariate.txt")
FandQ[method == "Inverse variance weighted", N:=  Q_df + 1]
FandQ <- FandQ[method == "Inverse variance weighted",]
egger_intercept <- fread("Data/Modified/egger_intercept.txt")
tissuespec <- fread("Data/Modified/snptissuespecificity.txt")
res_simple <- fread( "Data/Modified/res_univariate.txt")


######
return_format_data<-function(data) {
  return(data[, paste0(round(exp(b), digits = 2), " 95% CI=", round(exp(lci), digits = 2), "-",  round(exp(uci), digits = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_data_noexp <-function(data) {
  return(data[, paste0(round(b, digits = 2), " SD (95% CI=", round(lci, digits = 2), "-",  round(uci, digits = 2), ", p=",pval %>% formatC(., format = "e", digits = 1), ")")])
}
return_format_fstat <-function(data) {
  return(data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")])
}
#Abstract
datatable[,.(trait, sample_size)]
exposures_gtex_hypo <- fread( "Data/Modified/exposures_gtex_hypothalamus.txt", nrows = 1)
exposures_gtex_hypo$samplesize.exposure
res_steiger[method == "Inverse variance weighted", ][exposure == "bmi_ukbgiant" & outcome == "Fasting_Insulin", ] %>%
  return_format_data_noexp(.)
res_steiger[method == "Inverse variance weighted", ][exposure == "whradjbmi" & outcome == "Fasting_Insulin", ] %>%
  return_format_data_noexp(.)
res_steiger[method == "Inverse variance weighted", ][exposure == "Fasting_Insulin" & outcome == "bmi_ukbgiant", ] %>%
  return_format_data_noexp(.)
{
oui<-mvmr_results[outcome == "Fasting_Insulin" & exposure == "bmi_ukbgiant" & method %in% "Multivariable IVW", b ]
non<-mvmr_results[outcome == "Fasting_Insulin" & exposure == "wc_ukb" & method %in% "Multivariable IVW", b ]
non/oui
  }

###########METHOD##################    
eQTLcoloc[hgnc_symbol %in% c("ADCY5", "TCF7L2", "TRIM73", "PMS2P3")]$posprob_colocH4.SNP
eQTLcoloc[hgnc_symbol %in% c("ADCY5", "TCF7L2", "TRIM73",  "PMS2P3")]$hgnc_symbol
ldmat <- ieugwasr::ld_matrix_local(eQTLcoloc[hgnc_symbol%in%c("PMS2P3","TRIM73"), posprob_colocH4.SNP], plink_bin = genetics.binaRies::get_plink_binary(), 
                                   bfile = ldref)
ldmat^2
eQTLcoloc[hgnc_symbol %in% c("ADCY5", "TCF7L2", "TRIM73"), min(posprob_colocH4.SNPexplained_var)]



###############Results##################
#para 1
return_format_fstat(FandQ[exposure == "bmi_ukbgiant" & outcome == "Fasting_Insulin",])
return_format_data_noexp(res_steiger[exposure == "bmi_ukbgiant" & outcome == "Fasting_Insulin" & method == "Inverse variance weighted"])
return_format_fstat(FandQ[exposure == "bmi_ukbgiant_bio" & outcome == "Fasting_Insulin",])
return_format_data_noexp(res_steiger[exposure == "bmi_ukbgiant_bio" & outcome == "Fasting_Insulin" & method == "Inverse variance weighted"])

#para 2
return_format_fstat(FandQ[exposure == "whradjbmi" & outcome == "Fasting_Insulin",])
return_format_data_noexp(res_steiger[exposure == "whradjbmi" & outcome == "Fasting_Insulin" & method == "Inverse variance weighted"])
return_format_data_noexp(res_steiger[exposure == "mri_gfat" & outcome == "Fasting_Insulin" & method == "Inverse variance weighted"])
return_format_data_noexp(res_steiger[exposure == "mri_vat" & outcome == "Fasting_Insulin" & method == "Inverse variance weighted"])
return_format_data_noexp(res_steiger[exposure == "mri_asat" & outcome == "Fasting_Insulin" & method == "Inverse variance weighted"])

mvmr_results[method == "Multivariable IVW" & exposure == "UKB-b-9405", ] %>% return_format_data_noexp(.)
mvmr_results[method == "Multivariable IVW" & exposure == "bmi_ukbgiant" & otherexposure == "UKB-b-9405", ] %>% return_format_data_noexp(.)
mvmr_results[method == "Multivariable IVW" & exposure == "bmi_ukbgiant" & otherexposure == "UKB-b-9405", ][, .(cochranQ,cochranQpval)] %>% distinct
mvmr_results[method == "Multivariable Egger Intercept"& exposure == "bmi_ukbgiant" & otherexposure == "UKB-b-9405"][1,.(b,pval)]
# {mediation <- mvmr_results[outcome == "Fasting_Insulin" & exposure == "bmi_ukbgiant",]
#   total<- mediation[method %in% c("Inverse variance weighted"), b] 
#   direct<-mediation[method %in% "Multivariable IVW", b]
#   (total-direct)/total}
{
  oui<-mvmr_results[outcome == "Fasting_Insulin" & exposure == "bmi_ukbgiant" & otherexposure == "UKB-b-9405" & method %in% "Multivariable IVW", b ]
  non<-mvmr_results[outcome == "Fasting_Insulin" & exposure == "UKB-b-9405" & method %in% "Multivariable IVW", b ]
  non/oui
}

#para 3
return_format_fstat(FandQ[exposure == "Fasting_Insulin" & outcome == "whradjbmi",])
return_format_data_noexp(res_simple[exposure == "Fasting_Insulin" & outcome ==  "bmi_ukbgiant" & method == "Inverse variance weighted"])
return_format_data_noexp(res_simple[exposure == "Fasting_Insulin" & outcome ==  "whradjbmi" & method == "Inverse variance weighted"])
k <- GagnonMR::from_genecard_to_generegion("FTO")
return_format_data_noexp(res_steiger[exposure == "Fasting_Insulin" & outcome ==  "bmi_ukbgiant" & method == "Inverse variance weighted"])
#para 4
return_format_fstat(FandQ[exposure == "Fasting_Insulin_bio" & outcome ==  "bmi_ukbgiant",])
return_format_data_noexp(res_steiger[exposure == "Fasting_Insulin_bio" & outcome ==  "bmi_ukbgiant" & method == "Inverse variance weighted"])
return_format_data_noexp(res_steiger[exposure == "Fasting_Insulin_bio" & outcome ==  "mri_asat" & method == "Inverse variance weighted"])
return_format_data_noexp(res_steiger[exposure == "Fasting_Insulin_bio" & outcome ==  "mri_gfat" & method == "Inverse variance weighted"])
return_format_data_noexp(res_steiger[exposure == "Fasting_Insulin_bio" & outcome ==  "mri_vat" & method == "Inverse variance weighted"])


####################Discussion######################
{
  oui<-mvmr_results[outcome == "Fasting_Insulin" & exposure == "bmi_ukbgiant" & otherexposure == "UKB-b-9405" & method %in% "Multivariable IVW", b ]
  non<-mvmr_results[outcome == "Fasting_Insulin" & exposure == "UKB-b-9405" & method %in% "Multivariable IVW", b ]
  non/oui
}
