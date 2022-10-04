#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

eQTLcoloc <- fread( "Data/Modified/eQTLcoloc.txt")
eQTLcoloc <- eQTLcoloc[outcome != "Fasting_Insulin_correctedBMI",]
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")
gencode[, gene_id2 := gsub("\\..*", "", gene_id)]
eQTLcoloc <- merge(eQTLcoloc, distinct(gencode[, .(gene_id2, gene_name)]), by.x = "exposure", by.y = "gene_id2", all.x = TRUE)
vep_input <- eQTLcoloc[!is.na(gene_name), posprob_colocH4.SNP]
all_out_vcf <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", rsid = vep_input)

VariantAnnotation::writeVcf(all_out_vcf, file = "Data/Modified/VEPinput.vcf")
fwrite(list(vep_input), "Data/Modified/VEPinput.txt")

system2(command = "sh", args = "Analysis/2d_VEP.sh")
# then go see this website https://www.ensembl.org/info/docs/tools/vep/index.html
#follow everything with default settings except change for a 1Mb window
vepoutput <-fread("Data/Modified/VEPoutput.txt", skip = 39)
{test <- all_out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table
  setnames(vepoutput, "#Uploaded_variation", "Uploaded_variation")
  merge(test[,.(SNP, chr.exposure, pos.exposure)], vepoutput[,.(Uploaded_variation, Location),], 
        by.x = "SNP", by.y = "Uploaded_variation") %>% distinct(.)}
# vepoutput<- vepoutput[Gene %in% eQTLcoloc[!is.na(gene_name), exposure], ]
vepoutput <- separate(vepoutput, "Consequence", sep = ",", into = paste0("consequence", 1:4))

vepoutput <- vepoutput[, unique(c(consequence1, consequence2, consequence3, consequence4)), by =c("Uploaded_variation", "Gene")]
setnames(vepoutput, "V1", "consequence")
vepoutput<-vepoutput[!is.na(consequence)]
vepoutput[, consequences := paste(consequence, collapse = ","), by = c("Uploaded_variation", "Gene")]
setnames(vepoutput, "Uploaded_variation", "SNP")
df_VEP <- distinct(vepoutput[, c("SNP", "Gene", "consequences")])
df_VEP <- df_VEP[SNP %in% vep_input, ]
toMatch<- c("missense_variant", "stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift")
df_VEP[ , is_altering_variant := grepl(paste(toMatch,collapse="|"), consequences)]
fwrite(df_VEP, "Data/Modified/df_VEP.txt")

#
message(paste0("This script finished without errors at ", Sys.time()))
id_to_include <- df_index[pmid == 32999275, ][trait %in% c("ENSG00000127946", "ENSG00000127957","ENSG00000178809") ,]$id
vcf_path <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id_to_include, "/", id_to_include,".vcf.gz")
k <- map(as.list(vcf_path), ~ gwasvcf::query_gwas(vcf = .x, chrompos = df_index[id %in% id_to_include, note]) %>% 
           gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table(.)) %>% rbindlist(., fill = TRUE)

k
mk <- merge(eQTLcoloc[!is.na(gene_name), .(posprob_colocH4.SNP, exposure, hgnc_symbol)], df_VEP, by.x = "posprob_colocH4.SNP", by.y = "SNP")

k[SNP %in% mk$posprob_colocH4.SNP,]
