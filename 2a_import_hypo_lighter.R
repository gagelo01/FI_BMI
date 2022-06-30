#!/usr/bin/env Rscript
library(biomaRt)
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"


tissue_dir <- list.files("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/GTEx_EUR_Analysis_v8_eQTL_expression_matrices")
tissue_loop <- sub(".v8.normalized_expression_EUR_chr1_22.bed","",tissue_dir)
vec_tissue <- tissue_loop[15]
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")

translation_chr <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt", select = "chr")
colnoms <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt", nrows = 1)
translation_colnoms <- colnames(colnoms)

traduction_chr <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt", select = "chr")
colnoms <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt", nrows = 1)
traduction_colnoms <- colnames(colnoms)

gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")

#####
fread_wrapper_extract <- function(data_path, index, dt_chr, colnoms) {
  salsa <- which(dt_chr == index)
  valsa <- (1:dt_chr[,.N] %in% salsa)
  seq <- rle(valsa)
  idx  <- c(0, cumsum(seq$lengths))[which(seq$values)] + 1
  indx <- data.frame(start=idx, length=seq$length[which(seq$values)])
  result <- do.call(rbind,apply(indx,1, function(x) return(fread(data_path,nrows=x[2],skip=x[1]))))
  colnames(result)<-colnoms
  return(result)
}
######


options(future.globals.maxSize= 1e10)
plan(multisession, workers = 6)

errors <- vector(mode = "list", length = 22)
exposures_gtex <- vector(mode = "list", length = 22)
for(i in 1:22) {
arguments <- tidyr::crossing(tissue = vec_tissue, 
                             gene = gencode[gene_type =="protein_coding" & chr %in% i, unique(gene_name)])
arguments_split <- split(arguments, 1:nrow(arguments))

get_eQTL_safely <- safely(get_eQTL)

dfeqtl <- future_map(arguments_split, function(x) {get_eQTL_safely(tissue = x[names(x)=="tissue"], 
                                                                   gene = x[names(x)=="gene"])},
                     .options = furrr_options(seed = TRUE))# %>% 
#rbindlist(.,fill = TRUE)

index <- sapply(dfeqtl, function(x) is.null(x$error))
hypo <- lapply(dfeqtl[index], function(x) x$result) %>% rbindlist(., fill = TRUE)
hypo <- hypo[!is.na(probe),]
hypo[, c("chr", "start", "end", "strand", "n_cis_variant_tested", "distance_probe_variant", "end_top_variant", "top_variant", "zscore", "tissue.tissue", "hgncgene.gene") := NULL]

errors[[i]] <- dfeqtl[!index]

translation <- fread_wrapper_extract(data_path = "/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt",
                                     index = i,
                                     dt_chr = translation_chr,
                                     colnoms = translation_colnoms)

traduction <- fread_wrapper_extract(data_path = "/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt",
                                    index = i,
                                    dt_chr = traduction_chr,
                                    colnoms = traduction_colnoms)

exposures_formated <- format_gtex_data(exposures = hypo, translation = translation, gencode = gencode)

exposures_formated<-merge(exposures_formated, traduction[,.(rsid, chr, position)], by.x = "SNP", by.y = "rsid")
exposures_formated[, chr.exposure := chr]
exposures_formated[, pos.exposure := position]
exposures_formated[, chr := NULL]
exposures_formated[,position := NULL]
exposures_formated[, gene.exposure := NULL]
exposures_formated<- separate(exposures_formated, col = "id.exposure", into = c("id.exposure", "gene.exposure"), sep = "-")
exposures_gtex[[i]] <- exposures_formated

}

exposures_gtex <- rbindlist(exposures_gtex, fill = TRUE)
fwrite(exposures_gtex, "Data/Modified/exposures_gtex_hypothalamus.txt")
saveRDS(errors, "Data/Modified/hypo_errors.rds")

message("This script finished without errors")
