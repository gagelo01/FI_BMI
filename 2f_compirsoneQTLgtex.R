#!/usr/bin/env Rscript
library(biomaRt)
library(data.table)
library(tidyverse)
library(GagnonMR)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"


tissue_dir <- list.files("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/GTEx_EUR_Analysis_v8_eQTL_expression_matrices")
tissue_loop <- sub(".v8.normalized_expression_EUR_chr1_22.bed","",tissue_dir)
vec_tissue <- tissue_loop[c(31,1,2,23,24,34)]
eQTLcoloc <- fread( "Data/Modified/eQTLcoloc.txt")
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")
gencode[, gene_id2 := gsub("\\..*", "", gene_id)]
eQTLcoloc <- merge(eQTLcoloc, distinct(gencode[, .(gene_id2, gene_name)]), by.x = "exposure", by.y = "gene_id2", all.x = TRUE)
eQTLcoloc[,hgnc_symbol:=NULL]
# eQTLcoloc_hypo <- fread( "Data/Modified/eQTLcoloc_hypo.txt")
# eQTLcoloc_hypo[, gene_name := gsub("Brain_Hypothalamus-", "", exposure)]
# eQTLcoloc <- rbind(eQTLcoloc, eQTLcoloc_hypo, fill = TRUE)
# mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl_hugo<- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = eQTLcoloc$exposure, mart = mart)     
# ensembl_hugo
# 
# 
# gene_hugo_mart <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id" ), 
#                                    filters = "ensembl_gene_id ", values = eQTLcoloc$exposure, mart = ensembl) %>% 
#   as.data.table(.)

arguments <- tidyr::crossing(tissue = vec_tissue, gene = eQTLcoloc[!is.na(gene_name), unique(gene_name)])
arguments_split <- split(arguments, 1:nrow(arguments))


dfeqtl <- map(arguments_split, function(x) get_eQTL(tissue = x[names(x)=="tissue"], 
                                                    gene = x[names(x)=="gene"])) %>% 
  rbindlist(.,fill = TRUE)

#
translation <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt")
traduction <-  fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")

exposures_formated <- format_gtex_data(exposures = dfeqtl, translation = translation, gencode = gencode)

exposures_formated<-merge(exposures_formated, traduction[,.(rsid, chr, position)], by.x = "SNP", by.y = "rsid")
exposures_formated[, chr.exposure := chr]
exposures_formated[, pos.exposure := position]
exposures_formated[, chr := NULL]
exposures_formated[,position := NULL]
exposures_formated[, gene.exposure := NULL]
exposures_formated<- separate(exposures_formated, col = "id.exposure", into = c("id.exposure", "gene.exposure"), sep = "-")
exposures_gtex <- exposures_formated
fwrite(exposures_gtex, "Data/Modified/exposures_gtex_hyprcoloc.txt")

message("This script finished without errors")

