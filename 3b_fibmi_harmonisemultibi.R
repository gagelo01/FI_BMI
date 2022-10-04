#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(writexl)
library(furrr)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/FI_BMI")
all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr.txt")
all_inst_mvmr[,id.exposure := exposure]
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvmr.txt" )
all_outcome_mvmr[,id.outcome := outcome]
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")

##############
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
ao_small[, tomirge := tolower(id)]
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)#plan(multicore, workers = 20, gc = TRUE)#plan(sequential)
#######################
#change only this section 
#univariable
idserver_exposure <- c("trait-1-1", "trait-10-1", "trait-10-3", paste0("trait-25-", 1:9))
idserver_outcome <- c("trait-2-2")
exposure_vec<-c("bmi_ukbgiant_bio")
outcome_vec<-c("Fasting_Insulin_bio")
idmrbase_exposure <- "ukb-b-9405" #c("ukb-b-9405", "ukb-b-19953")
idmrbase_outcome <- NULL
#multivariable
k<-"bmi_ukbgiant"
exposurefor <- c("whr", "UKB-b-9405") #"UKB-b-9405" #The phenotype you want to correct for , but will be added as "exposure" which means there instruments will be selected
correctfor <- NULL #The phenotype you want to correct for , but will be added as "correctfor" which means there instruments won't selected
if(!is.null(exposurefor)){k <- purrr::cross2(k,exposurefor)}
mvmr_exposure<-lapply(k, unlist)
mvmr_outcome <- c("Fasting_Insulin" )
pval <- c(5e-8)
if(is.null(correctfor)){correctfor<-"NULL"}
arguments_mvmr <- purrr::cross(list(mvmr_exposure, mvmr_outcome, correctfor, pval))
#######
u<-c("logOR", "log odds", "log orr")
k <- c(df_index[unit %in% u,trait],  ao_small[unit %in% u, id])
exp_vec<-c(all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_exposure), unique(exposure)],
           df_index[id %in% idserver_exposure, trait])
out_vec<-c(all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_outcome), unique(exposure)],
           df_index[id %in% idserver_outcome, trait])
arguments_uni <- rbind(tidyr::crossing(exposure = c(exp_vec[!(exp_vec%in%k)], exposure_vec), outcome = out_vec),
                       tidyr::crossing(exposure = c(out_vec[!(out_vec%in%k)], outcome_vec), outcome = exp_vec))

arguments_uni <- distinct(arguments_uni)
setDT(arguments_uni)
arguments_uni <- arguments_uni[!(exposure == outcome),]
harm_univariate <- map(split(arguments_uni, 1:nrow(arguments_uni)), function(x) {
  harm <- TwoSampleMR::harmonise_data(exposure_dat = inst_all_sign_clump[exposure == x[,exposure], ],
                                      outcome_dat = all_outcome_mvmr[outcome == x[,outcome], ],
                                      action = 1)
  return(harm)}) %>% rbindlist(.,fill = TRUE)
harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate)
setDT(harm_univariate)
fwrite(harm_univariate, "Data/Modified/harm_univariate.txt")
harm_univariate[, exposure_outcome :=paste0(exposure, "_", outcome)]
harm_steiger <- harm_univariate[!(steiger_dir == FALSE & steiger_pval < 0.05),]

egger_intercept <- mr_pleiotropy_test(harm_steiger)

list_harm<- list(harm_univariate= harm_univariate, harm_steiger=harm_steiger)
list_res<-vector(mode="list", length = length(list_harm)) 
names(list_res)<-names(list_harm)
for(i in 1:length(list_harm)) {
list_harm_univariate <- split(list_harm[[i]], list_harm[[i]]$exposure_outcome)
list_res_univariate <- future_map(list_harm_univariate,
                                  function(x) {
                                    GagnonMR::all_mr_methods(x,
                                                             skip_presso = FALSE)},
                                  .options = furrr_options(seed = TRUE))

res_univariate <- rbindlist(list_res_univariate, fill = TRUE)
res_univariate[, c("id.exposure", "id.outcome") := NULL]
list_res[[i]] <- res_univariate
}

names(list_res)<- gsub("harm", "res", names(list_res))
list_harm_steiger <- split(harm_steiger,by = "exposure_outcome")
FandQ <- lapply(list_harm_steiger, function(x) {
  res <- TwoSampleMR::mr_heterogeneity(x)
  if(nrow(res)==0){res<-data.frame(exposure = x$exposure, outcome = x$outcome)}
  x <- TwoSampleMR::add_rsq(x)
  res$fstat<-GagnonMR::fstat_fromdat(list(x))
  res$rsq <- sum(x$rsq.exposure)
  return(res)
}) %>% rbindlist(.,fill = TRUE)

FandQ[,c("id.exposure", "id.outcome") := NULL]
fwrite(FandQ, "Data/Modified/FandQ_univariate.txt")
fwrite(list_res$res_univariate, "Data/Modified/res_univariate.txt")
fwrite(list_res$res_steiger, "Data/Modified/res_steiger.txt")
setDT(egger_intercept)
fwrite(egger_intercept, "Data/Modified/egger_intercept.txt")
#mvmr
performmvmr <- function(exposure_vec, outcome_vec, correctfor=NULL,pval_threshold = 1,
                        clump_exp_arg = "none", 
                        clump_r2 = 0.001, clump_kb = 10000) { #clump_exp_arg either none, first, or second
  message(paste0("MVMR for the effect of ", paste(exposure_vec, collapse = " + "), " on ", outcome_vec, " while correcting for ",  paste(correctfor, collapse = " + ")))
  exposure_dat <- inst_all_sign_clump[inst_all_sign_clump$exposure %in% exposure_vec,]
 if(!is.null(correctfor)){k<-all_inst_mvmr[(all_inst_mvmr$exposure %in% correctfor) & (all_inst_mvmr$SNP %in% unique(exposure_dat$SNP)),]} else {k<-NULL }
  
  exposure_dat<- rbind(exposure_dat,k, fill = TRUE)
  d1 <- all_inst_mvmr[(all_inst_mvmr$exposure %in% c(exposure_vec,correctfor)) & (all_inst_mvmr$SNP %in% unique(exposure_dat$SNP)),]
  
  if(clump_exp_arg == "none") {clump_exp<-NULL} else if(clump_exp_arg == "first"){clump_exp<-exposure_vec[1]} else if(clump_exp_arg == "second"){clump_exp<-exposure_vec[2]}
  inst_mvmr <- prepare_for_mvmr(exposure_dat = exposure_dat, d1 =d1,clump_r2 = clump_r2, clump_kb = clump_kb, pval_threshold = pval_threshold, clump_exp = clump_exp, harmonise_strictness = 1)
  
  exposure_outcome_harmonized <- TwoSampleMR::mv_harmonise_data(exposure_dat = inst_mvmr,
                                                                outcome_dat = all_outcome_mvmr[all_outcome_mvmr$outcome == outcome_vec,],
                                                                harmonise_strictness = 1)
  mvmr_results <- GagnonMR::mv_multiple_MendelianRandomization(exposure_outcome_harmonized = exposure_outcome_harmonized)
  mvmr_results$clump_exposure <- mvmr_results$clump_exp %>% ifelse(is.null(.), "none", .)
  mvmr_results <- mvmr_results[!(mvmr_results$exposure %in% correctfor),]
  mvmr_results$correctfor <- paste(correctfor, collapse = " + ")
  return(mvmr_results)
}


resmvmr <- future_map(arguments_mvmr, function(x) {
  resnone<-performmvmr(exposure_vec = x[[1]], 
                       outcome_vec = x[[2]], 
                       correctfor=if(x[[3]]=="NULL"){NULL}else{x[[3]]}, 
                       clump_exp_arg = "none")
  return(resnone)
}, .options = furrr_options(seed = TRUE))

names(resmvmr) <- sapply(arguments_mvmr, function(x) paste0(paste(x[[1]],collapse = " + "), " ~ ", x[[2]], " correctfor = ", x[[3]]," (pval=", x[[4]],")"))
mvmr_results <- lapply(resmvmr, function(x)
  x[, otherexposure := apply(.SD, 1,function(x) paste(setdiff(unique(exposure), x), collapse = "+")), .SDcols = "exposure"])
mvmr_results<-rbindlist(mvmr_results)
fwrite(mvmr_results, "Results/resmvmr.txt")
saveRDS(resmvmr, "Data/Modified/res_mvmr.rds")

message("This script finished without errors")


