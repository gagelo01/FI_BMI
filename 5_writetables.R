#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library("xlsx")
library(writexl)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao<-fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
res_steiger <- fread( "Data/Modified/res_steiger.txt")
mvmr_results <- fread("Results/resmvmr.txt")
FandQ <- fread( "Data/Modified/FandQ_univariate.txt")
egger_intercept <- fread("Data/Modified/egger_intercept.txt")
tissuespec <- fread("Data/Modified/snptissuespecificity.txt")

####
cleanify<-function(data) {
    if("exposure"%in%names(data)) {
    data[, exposure := gsub("UKB-b-9405",  "WC_UKB", exposure)]
    data <- data[exposure != "whr", ]}
  if("outcome"%in%names(data)) {
    data[, outcome := gsub("UKB-b-9405",  "WC_UKB", outcome)]
    data<-data[outcome != "whr", ]}
  data[,c("id.exposure","id.outcome"):=NULL]
  return(data)
}
#############

datatable <- rbind(ao[id == "ukb-b-9405", ], 
                   df_index[id %in% c("trait-10-1", "trait-1-1", "trait-2-2","eqtl-2-1","trait-25-1"), ], fill = TRUE)
datatable[, c("id","build","category","subcategory","ontology","mr","priority","sd","ncase","ncontrol") := NULL]
exposures_gtex_hypo <- fread( "Data/Modified/exposures_gtex_hypothalamus.txt", nrows = 1)
newrow <- data.table(trait = "hypothalamus_gene_expression", group_name = "publix", year = 2020, author = "The GTEx Consortium",
             consortium = "GTEX", sex = "Males and Females", population = "Mixed", unit = "SD", 
             sample_size = exposures_gtex_hypo$samplesize.exposure, pmid = "PMID: 32913098", initial_build = "HG38/GRCh38",
             adjustments = "5PC WGS sequencing platform (HiSeq 2000 or HiSeq X), WGS library construction protocol (PCR-based or PCR-free) and donor sex")
datatable[trait == "ENSG00000000419", trait := "Langerhans Islet's gene expression"]  
datatable <- rbind(datatable, newrow, fill = TRUE)
datatable[trait=="Waist circumference",adjustments:=c("age, sex, genetic array, and 10 principal component of ancestry")]
#
inst_all_sign_clump[inst_sel_strat == "biologically_driven", exposure:=paste0(exposure, "_bio")]
inst_all_sign_clump[,id.exposure := exposure]
inst_all_sign_clump[,exposure := gsub("UKB-b-9405","WC_UKB",exposure)]

#
res_steiger[,exposure := gsub("UKB-b-9405","WC_UKB",exposure)]
res_steiger[,outcome := gsub("UKB-b-9405","WC_UKB",outcome)]
tochange<- c("res_steiger")
for(i in 1:length(tochange)) {
  assign(tochange[i], get(tochange[i])[!(method %in% c("Robust adjusted profile score (RAPS)", "Weighted mode")), ])
}

#

map(list(mvmr_results, egger_intercept,res_steiger, mvmr_results, inst_all_sign_clump), function(x) assign(deparse(substitute(x)), cleanify(x)))

#table caption
dt_title <- data.table(title = paste0("Supplementary Table ", 1:7),
                       caption = c("Description of the datasets used.",
                                   "Instruments and relevant statistics",
                                   "Genetic instrument effect on gene level accross different tissues",
                                   "Instrument strength and heterogeneity statistics for bivariable MR",
                                   "Bivariable Mendelian Randomization results",
                                   "Bivariable Mendelian Randomization Egger's intercept",
                                   "Multivariable Mendelian randomization results."))

#                       
list_supdat <- list("Tables captions and titles" = dt_title,
                    "Supplementary Table 1" = datatable,
                    "Supplementary Table 2" = inst_all_sign_clump,
                    "Supplementary Table 3" = tissuespec,
                    "Supplementary Table 4" = FandQ,
                    "Supplementary Table 5" = res_steiger, 
                    "Supplementary Table 6" = egger_intercept, 
                    "Supplementary Table 7" = mvmr_results)
for(i in 1:length(list_supdat)) {
  writexl::write_xlsx(x = list_supdat[[i]],
                      path = paste0("Results/", gsub(" ", "_", names(list_supdat)[[i]]), ".xlsx"))
}
writexl::write_xlsx(x = list_supdat,
                    path = "Results/supplementary_tables_clean.xlsx")

message("This script finished without errors")

