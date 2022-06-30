#!/usr/bin/env Rscript
library(biomaRt)
library(data.table)
library(tidyverse)
library(GagnonMR)
library(gassocplot)
library(ggplotify)
library(cowplot)
library(ggrepel)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

eQTLcoloc <- fread( "Data/Modified/eQTLcoloc.txt")
eQTLcoloc <-  eQTLcoloc[outcome != "Fasting_Insulin_correctedBMI",]
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")
gencode[, gene_id2 := gsub("\\..*", "", gene_id)]
eQTLcoloc <- merge(eQTLcoloc, distinct(gencode[, .(gene_id2, gene_name)]), by.x = "exposure", by.y = "gene_id2", all.x = TRUE)

exposures_gtex <- fread( "Data/Modified/exposures_gtex_hyprcoloc.txt")
df_index<-fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

df_index_eqtl <- df_index[pmid == 34644572, ][trait %in% eQTLcoloc[!is.na(gene_name), exposure]]
df_index_eqtl <- df_index_eqtl[,.(id, note)]
df_index_eqtl[, ID_exposure_file := paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id, "/", id, ".vcf.gz")]
df_index_eqtl[,id:=NULL]
setnames(df_index_eqtl, "note", "chrompos")
df_index_eqtl_split <- split(df_index_eqtl, 1:nrow(df_index_eqtl))

inst_islet<- map(df_index_eqtl_split, function(x) {
inst_tsmr <- gwasvcf::query_gwas(vcf = x$ID_exposure_file, chrompos = x$chrompos) %>%
  gwasglue::gwasvcf_to_TwoSampleMR(.) %>%
  as.data.table(.)
return(inst_tsmr)}) %>% rbindlist(., fill = TRUE)

inst_islet_small <- merge(inst_islet, eQTLcoloc[,.(exposure, posprob_colocH4.SNP, gene_name)], 
      by.x = c("exposure", "SNP"), by.y = c("exposure", "posprob_colocH4.SNP"))


exposures_gtex_small <- merge(exposures_gtex, eQTLcoloc[,.(posprob_colocH4.SNP, gene_name)], 
      by.x = c("gene.exposure", "SNP"), by.y = c("gene_name", "posprob_colocH4.SNP"))

nametochange<-colnames(exposures_gtex_small)[grepl("exposure", colnames(exposures_gtex_small))]
setnames(exposures_gtex_small, nametochange, gsub("exposure", "outcome", nametochange))
exposures_gtex_small[,id.outcome := outcome]
harm <- TwoSampleMR::harmonise_data(inst_islet_small[!is.na(gene_name)], exposures_gtex_small, 1)
setDT(harm)
harm<-harm[,.(SNP, effect_allele.exposure, other_allele.exposure, outcome, beta.exposure, beta.outcome)]
harm <- separate(harm, col = outcome, into = c("tissue", "gene"), sep = "-")
setnames(harm, "beta.exposure", "beta.pancreatic_islet")
harm[,tissue := paste0("beta.", tissue)]

theres <- dcast(harm,  SNP+ effect_allele.exposure + other_allele.exposure+gene+beta.pancreatic_islet ~ tissue, value.var = "beta.outcome")

#PMS2 te rsid is not in GTEX, and apparently there are no proxies
# {ldmat <-  exposures_gtex[gene.exposure == "PMS2",c(unique(SNP),"rs7798471") ] %>%
# ieugwasr::ld_matrix_local(., plink_bin = genetics.binaRies::get_plink_binary(), bfile = ldref)
# test <- as.data.frame(ldmat)
# test$rowname <- rownames(test)
# nom <- colnames(test)[grepl("rs7798471",colnames(test))]
# test <- test[,c("rowname", nom)]
# setDT(test)
# test[(rs7798471_C_T^2)>0.1,] #No good proxies, hence
# }
#faire mr
ldmat <- ieugwasr::ld_matrix_local(eQTLcoloc[!is.na(hgnc_symbol), posprob_colocH4.SNP], plink_bin = genetics.binaRies::get_plink_binary(), 
                                   bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")
ldmat
diag(ldmat)<- 0
max(ldmat^2) #yes
ldmat
theres #rs1167827 is notspecific therefore remove ?

dat_vcf <- gwasvcf::query_gwas("/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz", 
                               rsid = eQTLcoloc[!is.na(gene_name) & gene_name != "PMS2P3", unique(posprob_colocH4.SNP)],
                               proxies = "no")
inst_map <- dat_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(. , "exposure") %>% data.table::as.data.table(.)


id_out<-c("trait-14-8","trait-14-6", "trait-14-7","trait-14-10", "trait-1-1", "trait-10-1", "trait-16-4", "dis-6-1")
path_out<-c(paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id_out, "/", id_out, ".vcf.gz"),
            "/mnt/sda/gagelo01/Vcffile/MRBase_vcf/ukb-b-9405/ukb-b-9405.vcf.gz")
out <- map(path_out, function(x) GagnonMR::extract_outcome_variant(snps = inst_map$SNP, outcomes = x)) %>%
  rbindlist(., fill = TRUE)
resharm <- TwoSampleMR::harmonise_data(exposure_dat = inst_map, outcome_dat = out, action = 1) %>% 
  as.data.table(.)

resharm <- TwoSampleMR::add_rsq(resharm)
resharm <- TwoSampleMR::steiger_filtering(resharm)
fwrite(resharm, "Data/Modified/mapharm.txt")

resmap <- resharm %>%
  TwoSampleMR::mr(., method = "mr_ivw")

##Faire hyprcoloc
#import fi
inst_fi<- map(as.list(df_index_eqtl$chrompos), function(x) {
  inst_tsmr <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-2-2/trait-2-2.vcf.gz",
                                   chrompos = x) %>%
    gwasglue::gwasvcf_to_TwoSampleMR(.) %>%
    as.data.table(.)
  return(inst_tsmr)}) %>% rbindlist(., fill = TRUE)
inst_fi[,id.exposure:=exposure]
inst_fi <- distinct(inst_fi)
inst_bmi<- map(as.list(df_index_eqtl$chrompos), function(x) {
  inst_tsmr <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz",
                                   chrompos = x) %>%
    gwasglue::gwasvcf_to_TwoSampleMR(.) %>%
    as.data.table(.)
  return(inst_tsmr)}) %>% rbindlist(., fill = TRUE)
inst_bmi[,id.exposure := exposure]
inst_bmi <- distinct(inst_bmi)

exposures_gtex[SNP %in% eQTLcoloc[!is.na(gene_name), posprob_colocH4.SNP], ][,unique(SNP) %>% length,by = "id.exposure"]
inst_bmi[SNP %in% eQTLcoloc[!is.na(gene_name), posprob_colocH4.SNP], ][,unique(SNP) %>% length,by = "id.exposure"]

#hyprcoloc
stack_assoc_plot_wrapper <- function(df_aligned, res_hypr, annotate_snp=NULL) {
  stopifnot(is.factor(df_aligned$exposure))
  df_reshaped <- reshape(df_aligned, idvar = c("SNP", "chr.exposure", "pos.exposure"), timevar = "exposure", direction = "wide")
  
  ldmat <- ieugwasr::ld_matrix_local(
    df_reshaped$SNP,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = ldref
  )
  snpname <- do.call(rbind, strsplit(rownames(ldmat), split = "_"))[,1]
  df_reshaped <- df_reshaped[SNP %in% snpname, ]
  
  markers <- df_reshaped[, .(SNP, chr.exposure, pos.exposure)]
  setnames(markers, colnames(markers), c("marker", "chr", "pos"))
  
  df_aligned[, z := beta.exposure/se.exposure]
  z<-reshape(df_aligned[SNP %in% snpname,.(SNP, exposure, z)], idvar = "SNP", timevar = "exposure", direction = "wide")  
  z[,SNP:=NULL]  
  setnames(z, colnames(z), gsub("z.", "", colnames(z)))
  setcolorder(z, levels(df_aligned$exposure))
  zscores<-as.matrix(z)  
  
  top_snp <- res_hypr[traits != "None", ][1, candidate_snp]
  if(is.na(top_snp)){top_snp<-annotate_snp}
  
  res <- gassocplot::stack_assoc_plot(markers = markers,
                                      z = zscores,
                                      corr = ldmat, 
                                      traits= colnames(zscores),
                                      top.marker= top_snp)
  
  return(res)
  
}

sensitivity.plot_wrapper <- function(df_aligned) {
  df_reshaped <- reshape(df_aligned, idvar = c("SNP", "chr.exposure", "pos.exposure"), timevar = "exposure", direction = "wide")
  effect_est <- df_reshaped[, .SD, .SDcols = names(df_reshaped)[grepl("beta.exposure", names(df_reshaped))]] %>% as.matrix
  effect_se <- df_reshaped[, .SD, .SDcols = names(df_reshaped)[grepl("^se.exposure", names(df_reshaped))]] %>% as.matrix
  res<- hyprcoloc::sensitivity.plot(effect.est = effect_est,
                                    effect.se = effect_se,
                                    trait.names = gsub("beta.exposure.", "", colnames(effect_est), fixed = TRUE),
                                    snp.id = df_reshaped$SNP,
                                    similarity.matrix = TRUE)
  return(res)
}

drawheatmap <- function(heat) {
  levels <- colnames(heat)
  heat <- as.data.frame(heat)
  heat$row <- rownames(heat)
  rownames(heat)<-NULL  
  setDT(heat)
  heat <- melt(heat, id.vars = "row")
  heat[,row:=factor(row, levels = levels)] 
  heat[,variable:=factor(variable, levels = levels)]  
  
  
  
  
  g <- ggplot(heat, aes(x = variable, tissue, y = row, fill = value))  +
    geom_tile() +
    # scale_fill_gradient(low = "lightblue", high = "blue3",limits=c(0,1)) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF",limits=c(0,1)) +
    labs(fill = "") +
    theme(
      panel.background = element_blank(),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
      legend.position = "top",
      legend.title = element_text(
        color = "gray20",
        size = 12
      ),
      legend.text = element_text(
        color = "gray20",
        size = 10
      ),
      legend.title.align = 0.5,
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      axis.title = element_blank(),
      axis.line = element_line(size = 1, colour = "gray20"),
      axis.ticks = element_line(size = 1, colour = "gray20"),
      # axis.text.y = element_text(
      #   size = 10,
      #   colour = "gray20"
      # ),
      axis.text.y=element_blank(),
      axis.text.x = element_text(
        angle = 60,
        size = 10,
        hjust = 1,
        colour = "gray20"
      ),
      axis.ticks.length = unit(.25, "cm"))
  
  g
}

inst_liver <- exposures_gtex[id.exposure == "Liver",]
inst_liver[,id.exposure := exposure]
inst_liver <- distinct(inst_liver)

inst_islet <- merge(inst_islet, distinct(eQTLcoloc[,.(exposure, gene_name)]), by.x = "exposure", by.y = "exposure")
inst_islet[, exposure := paste0("Pancreatic_islet-", gene_name)]
inst_islet[,id.exposure := exposure]
setnames(inst_islet, "gene_name", "gene.exposure")

genename <- eQTLcoloc[!is.na(gene_name), gene_name]
for(i in 1:length(genename)) {
dt <- rbindlist(list(inst_islet[gene.exposure == genename[i], ], inst_liver[gene.exposure == genename[i]], inst_fi), fill = TRUE)
aligned <- prepare_for_mvmr(exposure_dat = dt, d1 = dt, harmonise_strictness = 1, should_clump = FALSE)
aligned[, chr.exposure := chr.exposure %>% as.character(.) %>% as.numeric(.)]
hyprres <- run_hypr_on_aligned(aligned)

k<-aligned$exposure %>% unique
levels <- c(k[grep("Pancreatic_islet", k)], k[grep("Liver", k)], "Fasting_Insulin")
aligned[,exposure := factor(exposure, levels = rev(levels))]

if(all(is.na(hyprres$candidate_snp))){
  annotate_snp<-eQTLcoloc[gene_name == genename[i], posprob_colocH4.SNP]
}else{annotate_snp<-NULL}

A <- stack_assoc_plot_wrapper(df_aligned = aligned, res_hypr = hyprres, annotate_snp=annotate_snp)
res <- sensitivity.plot_wrapper(df_aligned = aligned)
B<-drawheatmap(res[[2]])

twopanel <-  ggdraw() +
  draw_plot(ggplotify::as.ggplot(A) + theme(text = element_text(size = 0.4)), x = 0.08, y =0, width = .6, height = 1) +
  draw_plot(B, x = .65, y =0.1, width = .35, height = 0.7) +
  draw_plot_label(label = c("", ""), size = 25,
                  x = c(0, 0.62), y = c(0.9, 0.9))
saveRDS(ggplotify::as.ggplot(A), paste0("Results/stackassoc_plot_", genename[i], ".rds"))
saveRDS(object = twopanel, file = paste0("Results/twopanel_hypr_plot_", genename[i], ".rds"))
ggsave(plot = twopanel, filename = paste0("Results/", "twopanel_hypr_plot_", genename[i], ".png"),
       width = 590/72,height = 583/72,units="in",scale=1, device = "png")

}

####
source("Analysis/my_mr_scatter_plot.R")
genename<-c("TCF7L2", "ADCY5",  "TRIM73") #remove "PMS2P3" because in LD and none specific

resharm<-fread( "Data/Modified/mapharm.txt")
resharm <- merge(resharm, distinct(eQTLcoloc[gene_name %in% genename,.(posprob_colocH4.SNP, hgnc_symbol)]), by.x = "SNP", by.y = "posprob_colocH4.SNP")
resharm[, Locus := hgnc_symbol]
resharm[,id.outcome := outcome]
resharm[,id.exposure := exposure]
resharm[, exposure_outcome := paste0(exposure, "_", outcome)]
resharm_split <- split(resharm, resharm$exposure_outcome)
resmap <-  map(resharm_split, GagnonMR::all_mr_methods) %>%
  rbindlist(.,fill=TRUE)
q_int_f <- function(harm) {
  q <- TwoSampleMR::mr_heterogeneity(harm) %>% as.data.table
  int <- TwoSampleMR::mr_pleiotropy_test(harm) %>% as.data.table
  q<-q[method != "MR Egger",]
  q[,c("method", "id.exposure","id.outcome"):=NULL];
  setnames(int, c("se", "pval"),c("egger_intercept_se","egger_intercept_pval"))
  int[,c("id.exposure", "id.outcome") := NULL]
  q_int<- merge(int, q, by = c("exposure", "outcome"))
  q_int[,nsnp := harm[,.N]]
  q_int[,fstat := GagnonMR::fstat_fromdat(list(harm))]
  return(q_int)
}
q_int <-  map(resharm_split,   q_int_f) %>%
  rbindlist(.,fill=TRUE)

setDT(resmap)
k <- resmap[outcome == "bmi_ukbgiant" &  grepl("Inverse variance weighted", method),]
m <-resharm[outcome == "bmi_ukbgiant"]
m[, align:=""]
scatter <- my_mr_scatter_plot(mr_results = k, dat = m, equation_facet_grid = "") +
  ylab("SNP effect on body mass index") +
  xlab("SNP effect on fasting insulin") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# annotation <- data.table(x = 0.01, y = 0.015, label = paste0("b = ", round(k$b, digits =2),"; p = ", formatC(k$pval, format = "e", digits = 1)))
# scatter <- scatter +
#   geom_text(data=annotation, aes( x=x, y=y, label=label),                 , 
#               color="orange", size=7 , angle=45, fontface="bold" )

for(i in 1:length(genename)){
assign(paste0("X", i), 
       readRDS(file = paste0("Results/twopanel_hypr_plot_", genename[i], ".rds")))
  assign(paste0("Y", i), 
         readRDS(file = paste0("Results/stackassoc_plot_", genename[i], ".rds")))
}

fourpanel <- cowplot::plot_grid(scatter, X3, X1,X2, labels=c("A)", "B)", "C)", "D)"))

ggsave(plot = fourpanel, filename = "Results/fig2_fourpanel.png",
       width = 1000/72,height = 1000/72,units="in",scale=1, device = "png")

threepanel <- cowplot::plot_grid(Y1, Y2,Y3, labels=c("A)", "B)", "C)"), nrow = 1, ncol = 3)
threepanel

ggsave(plot = threepanel, filename = "Results/fig4_threepanel.png",
       width = 1000/72,height = 600/72,units="in",scale=1, device = "png")

fwrite(resmap, "Data/Modified/resmap.txt")
fwrite(resharm, "Data/Modified/harmmap_vinuela.txt")
fwrite(theres, "Data/Modified/snptissuespecificity.txt")



