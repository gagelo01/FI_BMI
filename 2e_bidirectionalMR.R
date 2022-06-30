#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index<- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
#####
shortharm<- function(inst, vcffile_out) {
  out_vcf <- gwasvcf::query_gwas(vcf = vcffile_out,
                                 rsid = inst$SNP, proxies = "yes", bfile = ldref, tag_kb = 10000, tag_r2 = 0.8)
  
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)
  
  harm <- TwoSampleMR::harmonise_data(inst, out_tsmr, action = 1)
  return(harm) 
}

######
idserver1<-c("trait-10-1", "trait-1-1") #whradjbmi, bmi, 
serverdata1<- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/",idserver1, "/", idserver1, ".vcf.gz")
mrbasedata<-"/mnt/sda/gagelo01/Vcffile/MRBase_vcf/ukb-b-9405/ukb-b-9405.vcf.gz"
serverdata2<- "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz"

argumentsforward<-data.table(exposure = serverdata2, outcome = c(mrbasedata,serverdata1))
argumentsreverse<-argumentsforward
colnames(argumentsreverse)<-ifelse(colnames(argumentsreverse) == "exposure", "outcome", "exposure") 
arguments <- rbind(argumentsforward, argumentsreverse)

#####harmonise
inst <- map(as.list(unique(arguments$exposure)), function(x) GagnonMR::get_inst(x)) %>% 
  rbindlist(., fill = TRUE)

outcome<-map(as.list(unique(arguments$outcome)), function(vcffile_out) {
  out_vcf <- gwasvcf::query_gwas(vcf = vcffile_out,
                                 rsid = unique(inst$SNP), proxies = "yes", bfile = ldref, tag_kb = 10000, tag_r2 = 0.8)
  out_tsmr <- out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% data.table::as.data.table(.)
  return(out_tsmr) })%>% rbindlist(., fill = TRUE)

k<-c(rep("Fasting_Insulin", 3), "UKB-b-9405", "whradjbmi", "bmi_ukbgiant")
arguments <- data.table(exposure = k, outcome = rev(k))

arguments_split <- split(arguments, 1:nrow(arguments))
harm <- map(arguments_split, function(x) {
  harm <- TwoSampleMR::harmonise_data(exposure_dat =  inst[exposure == x[,exposure], ],
                              outcome_dat = outcome[outcome == x[,outcome], ],
                              action = 1)
  return(harm)}) %>% rbindlist(.,fill = TRUE)


harm <- TwoSampleMR::steiger_filtering(harm)
setDT(harm)
harm[, id.exposure:= exposure]
harm[, id.outcome := outcome]
harm_split <- split(harm, 1:nrow(harm))
harm$fstat <- GagnonMR::fstat_fromdat(harm_split)
harm[,exposure_outcome := paste0(exposure, "_", outcome)]
fwrite(harm, "Data/Modified/harm.txt")
harm_split <- split(harm, harm$exposure_outcome)
res_simple <- map(harm_split, function(x) GagnonMR::all_mr_methods(x)) %>%
  rbindlist(., fill= TRUE)

harm[,exposure_outcome := paste0(exposure, "_", outcome)]
harm <- harm[!(steiger_dir == FALSE & steiger_pval < 0.05),]
harm_split <- split(harm, harm$exposure_outcome)
res_steiger <- map(harm_split, function(x) GagnonMR::all_mr_methods(x)) %>%
  rbindlist(., fill= TRUE)

fwrite(res_steiger, "Data/Modified/res_steiger.txt")
fwrite(res_simple, "Data/Modified/res_simple.txt")

###
res <- readRDS( "Data/Modified/vinuela_res_mapping")
index <- sapply(res, function(x) !is.null(x$error))
k <- rbindlist(lapply(res[!index], function(x)  x$result), fill = TRUE)
k <- k[posprob_coloc_PPH4 > 0.8, ]
ensembl_human <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 
ensemblhugo <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                                          filters = "ensembl_gene_id", values = k$exposure,
                                          mart = ensembl_human)
k <- merge(k, ensemblhugo, by.x = "exposure", by.y = "ensembl_gene_id", all.x = TRUE)
fwrite(k,"Data/Modified/eQTLcoloc.txt")
message("This script finished without errors")
