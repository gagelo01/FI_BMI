#!/usr/bin/env Rscript
# https://zenodo.org/record/3408356#.YjCY7XrMKUk
# InsPIRE
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
library(tictoc)
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()

gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")
gencode[,gene_id2 := gsub("\\..*","",gene_id)]
gencode[,length(unique(gene_id2)) , gene_id][V1>1,]
traduction <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction <- traduction[,!c("maf")]
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

test <- fread("Data/Raw/InsPIRE_Islets_Gene_eQTLs_Nominal_Pvalues.txt")
colnom<- colnames(test)[2:length(colnames(test))]
test[,FreqALT:=NULL]
colnames(test)<-colnom
test[, c( "GeneChr","TSS", "NumSNPtest","DistanceToGene", "FreqREF", "Lead", "T-stat" ) := NULL]
test[, GeneID := gsub("\\..*","",GeneID)]
test[!(GeneID %in% gencode$gene_id2), unique(GeneID)] #I believe those are pseudogene,
test <- merge(test, distinct(gencode[,.(gene_name,gene_id2)]),by.x = "GeneID",by.y="gene_id2", all.x = TRUE)
all_out <- merge(test, traduction, by.x = c("SNPchr", "SNPposition"), by.y = c("chr", "position"), all = FALSE)
all_out <- all_out[(ALT == a0 | ALT == a1) & (REF == a0 | REF == a1) & a0 != a1  & ALT != REF, ] #because low number removed, coded on the forward strand
all_out <- all_out[SNPchr %in% 1:22, ]
all_out[, SNPchr := as.integer(SNPchr)]
all_out[ALT == a0, Slope := Slope*-1]
all_out[ALT == a0, FreqALT := 1-FreqALT] 
all_out[ALT == a0, ALT := a1]
all_out[REF == a1, REF := a0] #less than mrbase, possibly because I do not use the same traduction file


df_traitnote <- all_out[,paste0(SNPchr, ":", min(SNPposition), "-", max(SNPposition)) %>% unique, by = "GeneID"]
df_traitnote <- merge(df_traitnote, all_out[, .N,by="GeneID"], by = "GeneID")

newrow <- data.frame(id = paste0("eqtl-2-",1:nrow(df_traitnote)), trait = df_traitnote$GeneID, group_name = "InsPIRE",
                     year = 2020,author = "Vinuela, Ana",consortium = "Fadista_et_al, Oxford_dataset, Varshney_et_al",sex = "Males and Females",
                     population = "European", unit = "SD", nsnp = df_traitnote$N, sample_size = 420,
                     initial_build = "HG19/GRCh37", category = "eqtl", pmid = 32999275, ncase = NA,
                     sd = 1, note = df_traitnote$V1, ncontrol = NA)
df_index <- rbind(df_index, newrow)
setDT(df_index)
list_eqtl <- split(all_out, all_out$GeneID)
rm(list = "all_out")

formateqtl_wrapper <- function(x, df_index) {
formattovcf_createindex(all_out = x,
                        snp_col = "rsID",
                        outcome_name = unique(x$GeneID),
                        beta_col = "Slope",
                        se_col = "SE",
                        pval_col = "Pvalue",
                        eaf_col = "FreqALT",
                        effect_allele_col = "ALT",
                        other_allele_col =  "REF",
                        ncase_col = NULL,
                        ncontrol_col = NULL,
                        samplesize_col = 420,
                        chr_col = "SNPchr",
                        pos_col = "SNPposition",
                        units = "SD",
                        traduction = NULL,
                        out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                        df_index = df_index,
                        group_name = "InsPIRE",
                        year = 2020,
                        author = "Vinuela, Ana",
                        consortium = "Fadista_et_al, Oxford_dataset, Varshney_et_al",
                        sex = "Males and Females",
                        population = "European",
                        initial_build = "HG19/GRCh37",
                        category = "eqtl",
                        pmid = 32999275,
                        note = note,
                        should_create_id = FALSE,
                        ID = df_index[pmid == 32999275 & trait == unique(x$Gene), id])}

df_index_copy <- df_index
tic()
options(future.globals.maxSize= 1e11)
plan(multisession, workers = 6)

future_map(list_eqtl, function(eqtl) {formateqtl_wrapper(x = eqtl, df_index = df_index_copy)},
           .options = furrr_options(seed = TRUE))

fwrite(df_index, "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

toc()

message("This script finished without errors")

