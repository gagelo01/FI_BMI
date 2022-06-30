#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library("xlsx")
library(writexl)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao<-fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
datatable <- rbind(ao[id == "ukb-b-9405", ], 
                   df_index[id %in% c("trait-10-1", "trait-1-1", "trait-2-2","eqtl-2-1"), ], fill = TRUE)
datatable[, c("id","build","category","subcategory","ontology","mr","priority","sd","ncase","ncontrol") := NULL]
exposures_gtex_hypo <- fread( "Data/Modified/exposures_gtex_hypothalamus.txt", nrows = 1)
newrow <- data.table(trait = "hypothalamus_gene_expression", group_name = "publix", year = 2020, author = "The GTEx Consortium",
             consortium = "GTEX", sex = "Males and Females", population = "Mixed", unit = "SD", 
             sample_size = exposures_gtex_hypo$samplesize.exposure, pmid = "PMID: 32913098", initial_build = "HG38/GRCh38")
datatable[trait == "ENSG00000000419", trait := "Langerhans Islet's gene expression"]  
datatable <- rbind(datatable, newrow, fill = TRUE)
datatable


res_steiger <- fread( "Data/Modified/res_steiger.txt")
resexcludecluster <- fread( "Data/Modified/resexcludecluster.txt") #"Results exclusing cluster 1and3" = resexcludecluster, I would not include that
dt_clust <-fread("Data/Modified/dt_clust.txt") #           "Supplementary Table 2" = "Results of mr-clust on the effect of fasting insulin on BMI",
resmap <- fread( "Data/Modified/resmap.txt")
mvmr_results <- fread( "Results/mvmrbmiwcfi.txt")
tissuespec <- fread("Data/Modified/snptissuespecificity.txt")
tochange<- c("res_steiger", "resmap")
for(i in 1:length(tochange)) {
  assign(tochange[i], get(tochange[i])[!(method %in% c("Robust adjusted profile score (RAPS)", "Weighted mode")), ])
}

map(list(res_steiger, resexcludecluster,resmap,mvmr_results), function(x) x[,c("id.exposure","id.outcome"):=NULL])

#table caption
table_caption <- data.table(table = paste0("Supplementary Table ",1:5),
                            caption = c("Description of the datasets used",
                              "Results bivariable with steiger filtering",
                              "Multivariable results to evaluate direct effect of general and abdominal adiposity on Fasting insulin",
                              "Genetic instrument effect on gene level accross different tissues",
           "Effect of insulin secretion on adiposity indices using a biological rational for genetic instrument selection"))

writexl::write_xlsx(x = list("Tables caption" = table_caption,
                             "Supplementary Table 1" = datatable,
                             "Supplementary Table 2" = res_steiger,
                             "Supplementary Table 3" = mvmr_results,
                             "Supplementary Table 4" = tissuespec,
                             "Supplementary Table 5" = resmap),
                    path = "Results/supplementary_tables.xlsx")

