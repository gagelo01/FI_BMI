library(data.table)
library(TwoSampleMR)
library(xlsx)
library(data.table)
library(ckbplotr)
library(ggforestplot)
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

res_steiger <- fread("Data/Modified/res_steiger.txt")
res_simple <- fread("Data/Modified/res_simple.txt")
resmapsteiger_hypo <- fread("Data/Modified/resmapsteiger_hypo.txt")
harm <- fread( "Data/Modified/harm.txt")
harm_is <- fread("Data/Modified/harm_is.txt")
res_is <- fread( "Data/Modified/res_is.txt")
harm_ir <- fread("Data/Modified/harm_ir.txt")
res_ir <- fread("Data/Modified/res_ir.txt")
mvmr_results <- fread("Results/mvmrbmiwcfi.txt")
resmap <- fread( "Data/Modified/resmap.txt")
res_at <- fread("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdistribution_FI/Data/Modified/res_univariate.txt")
res_at <- res_at[!(method %in% c("Weighted mode", "MR Egger", "Robust adjusted profile score (RAPS)")), ]
res_steiger <- rbind(res_steiger, res_at[exposure == "Fasting_Insulin" & outcome %in% c("Volume_VAT", "Volume_ASAT"),], fill = TRUE)

tochange<- c("res_ir", "res_is", "res_simple", "res_steiger", "resmap", "resmapsteiger_hypo")
for(i in 1:length(tochange)) {
assign(tochange[i], get(tochange[i])[!(method %in% c("Robust adjusted profile score (RAPS)", "Weighted mode")), ])
}
     
##Figure 1 and 4

res_univariate <- res_steiger

#Figure 1 
#Causal effet of BMI on FI using different genetic instruments selection strategy. 
#Using biologically informed genetic instrument and Using statistically informed genetic instruments
i<-1
out_name = c("Fasting_Insulin")
exp_name = c("bmi_ukbgiant")

doA <- res_univariate[exposure %in% exp_name[[i]] & outcome %in% out_name[i]]
doA[,instrument_selection := "Statistically informed"]
doB<-resmapsteiger_hypo[exposure %in% exp_name[[i]] & outcome %in% out_name[i]]
doB[,instrument_selection := "Biologically informed"]
doA <- rbindlist(list(doA, doB),fill = TRUE)
doA[,exposure := gsub("Fasting_Insulin", "Fasting Insulin", exposure)]
doA[,exposure := gsub("whradhbmi", "WHRadjBMI", exposure)]
doA[,exposure := gsub("bmi_ukbgiant", "BMI", exposure)]

resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))


mylabels <- data.frame(heading1 = doA$instrument_selection,
                       heading2 = doA$method,
                       heading3 = as.character(NA),
                       variable = as.character(1:nrow(doA)))


make_forest_plot(panels = list(resultsA),
                 col.key = "variable",
                 row.labels = mylabels,
                 exponentiate = FALSE,
                 pointsize = 2, 
                 rows = unique(mylabels$heading1),
                 col.stderr = NULL,
                 col.lci = "lci",
                 col.uci = "uci",
                 col.left         = c("n"),
                 col.left.heading = c("n SNP"),
                 col.right = "P_value",
                 col.right.heading = c("Effect (95% CI)", "P-value"),
                 xlab = c("Effect of 1-SD increase in BMI on fasting insulin (SD)"),
                 blankrows = c(0,1,0,0),
                 col.right.hjust = 1,
                 nullval = 0,
                 panel.headings = NULL,
                 xlim = c(0, 1))

ggsave(paste0("Results/", "Figure1", ".png"),
       width=550/72,height=328/72, units="in", scale=1,
       device = "png")
###plotting Figure 2 and supplementary figure 1
exp_name = c("whradjbmi")
out_name = c("Fasting_Insulin")
file_name = c("Figure2")
# panel_heading <- c("Effect of  WHRadjBMI on Fasting Insulin", "Effect of Fasting Insulin on WHRadjBMI")
panel_heading <- ""
i<-1
doA <- res_univariate[exposure %in% exp_name[[i]] & outcome %in% out_name[i]]
doA<- doA[method != "MR Egger",]
doB<-res_univariate[ outcome %in% exp_name[[i]] & exposure %in% out_name[i]]
doB<-doB[method != "MR Egger"]
doA[,exposure := gsub("Fasting_Insulin", "Fasting Insulin", exposure)]
doA[,exposure := gsub("whradhbmi", "WHRadjBMI", exposure)]
doA[,exposure := gsub("bmi_ukbgiant", "BMI", exposure)]

resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))

resultsB <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doB$b, digits =2),
                       lci =  round(doB$lci, digits = 2),
                       uci =  round(doB$uci, digits = 2),
                       n = doB$nsnp,
                       P_value = formatC(doB$pval, format = "e", digits = 1))

mylabels <- data.frame(heading1 = doA$method,
                       heading2 = as.character(NA),
                       heading3 = as.character(NA),
                       variable = as.character(1:nrow(doA)))


make_forest_plot(panels = list(resultsA),
                 col.key = "variable",
                 row.labels = mylabels,
                 exponentiate = FALSE,
                 pointsize = 2, 
                 rows = unique(mylabels$heading1),
                 panel.names = panel_heading,
                 col.stderr = NULL,
                 col.lci = "lci",
                 col.uci = "uci",
                 col.left         = c("n"),
                 col.left.heading = c("n SNP"),
                 col.right = "P_value",
                 col.right.heading = c("Effect (95% CI)", "P-value"),
                 xlab = c("Effect of 1-SD Increase of WHRadjBMI on Fasting Insulin (SD)"),
                 blankrows = c(0,1,0,0),
                 col.right.hjust = 1,
                 nullval = 0,
                 xlim = c(0,0.4))

ggsave(paste0("Results/", file_name[i], ".png"),
       width=541/72,height=200/72, units="in", scale=1,
       device = "png")
#Figure 3
unimeth<-"Inverse variance weighted" 
multimeth<- c("Multivariable IVW", "Multivariable Median",
              "Multivariable Lasso", "Multivariable Egger")             
data <- mvmr_results[method %in% c(unimeth, multimeth),]
data[, Category_other := ifelse(method %in% unimeth, "Univariable", "Multivariable")]
data[, Category_other := factor(Category_other, levels = c("Univariable", "Multivariable"))]
data[,name := gsub("bmi_ukbgiant", "Body mass index", exposure) %>% gsub("wc_ukb", "Waist circumference", .)]
data[Category_other == "Multivariable", name := name %>% ifelse(. == "Body mass index", "BMI adjusted for WC", .) %>% ifelse(. == "Waist circumference", "WC adjusted for BMI", .)]

forestplot(
  df = data,
  name = name,
  se = se,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Effect of 1 SD increase in WC/BMI on FI (SD)",
  ci = 0.95,
  colour = method,
  xlim = data[method != "MR Egger", round(c(min(lci),max(uci)), digits = 1)]
) + 
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = c(-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)) +
  theme(legend.position="right") +
  ggforce::facet_col(
    # facets = ~outcome_category,
    facets = ~Category_other,
    scales = "free_y",
    space = "free"
  )

ggsave(paste0("Results/", "Figure3colour", ".png"),
       width=350/72,height=250/72, units="in", scale=1,
       device = "png")

#clean

doA <- data
resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))


mylabels <- data.frame(heading1 = doA$name,
                       heading2 = doA$method,
                       heading3 = as.character(NA),
                       variable = as.character(1:nrow(doA)))

ckbplotr::make_forest_plot(panels = list(resultsA),
                           col.key = "variable",
                           row.labels = mylabels,
                           exponentiate = TRUE,
                           pointsize = 2, 
                           rows = unique(mylabels$heading1),
                           col.stderr = NULL,
                           col.lci = "lci",
                           col.uci = "uci",
                           col.left         = c("n"),
                           col.left.heading = c("n SNP"),
                           col.right = "P_value",
                           col.right.heading = c("Effect (95% CI)", "P-value"),
                           xlab = c("Effect of 1 SD increase in WC/BMI on FI (SD)"),
                           blankrows = c(0,0,0,0),
                           col.right.hjust = 1,
                           panel.headings = NULL,
                           scalepoints = FALSE)

ggsave(paste0("Results/", "Figure3Clean", ".png"),
       width=483/72,height=250/72, units="in", scale=1,
       device = "png")
#Figure 4
# derived in the script 2g_hyprcolocmap.

#Figure 5-6

exp_name = list(c("Fasting_Insulin"), c("Fasting_Insulin"))
out_name = list(c("bmi_ukbgiant", "whradjbmi"), c("Volume_VAT", "Volume_ASAT"))
list_name <- c("Figure5", "Figure6")
xlab = list(paste0("Effect of 1-SD increase in fasting insulin on ", c("BMI", "WHRadjBMI")),
            paste0("Effect of 1-SD increase in fasting insulin on ", c("Volume_VAT", "Volume_ASAT")))
for(i in 1:length(list_name)) {
doA<-res_steiger[exposure %in% exp_name[[i]] & outcome %in% out_name[[i]]]
doA[,instrument_selection := "Statistically informed"]
doB<-resmap[exposure %in% exp_name[[i]] & outcome %in% out_name[[i]]]
doB[,instrument_selection := "Biologically informed"]
doA <- rbindlist(list(doA, doB),fill = TRUE)
doA <- doA[method != "MR Egger",]
doA[,exposure := gsub("Fasting_Insulin", "Fasting Insulin", exposure)]
doA[,exposure := gsub("whradhbmi", "WHRadjBMI", exposure)]
doA[,exposure := gsub("bmi_ukbgiant", "BMI", exposure)]
doA[, name := gsub("_", " ", outcome)]
doB<- doA[outcome==out_name[[i]][2],]
doA <- doA[outcome==out_name[[i]][1],]
resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))

resultsB <- data.frame(variable = as.character(1:nrow(doB)),
                       estimate = round(doB$b, digits =2),
                       lci =  round(doB$lci, digits = 2),
                       uci =  round(doB$uci, digits = 2),
                       n = doB$nsnp,
                       P_value = formatC(doB$pval, format = "e", digits = 1))

mylabels <- data.frame(heading1 = doA$instrument_selection,
                       heading2 = doA$method,
                       heading3 = as.character(NA),
                       variable = as.character(1:nrow(doA)))


a<- make_forest_plot(panels = list(resultsA, resultsB),
                 col.key = "variable",
                 row.labels = mylabels,
                 exponentiate = FALSE,
                 pointsize = 2, 
                 rows = unique(mylabels$heading1),
                 col.stderr = NULL,
                 col.lci = "lci",
                 col.uci = "uci",
                 col.left         = c("n"),
                 col.left.heading = c("n SNP"),
                 col.right = "P_value",
                 col.right.heading = c("Effect (95% CI)", "P-value"),
                 xlab = xlab[[i]],
                 blankrows = c(0,1,0,0),
                 col.right.hjust = 1,
                 nullval = 0,
                 panel.headings = NULL)
a
ggsave(paste0("Results/", list_name[i], ".png"),
       width=864/72,height=250/72, units="in", scale=1,
       device = "png")
}

#Figure 5-6 combined

exp_name = list(c("Fasting_Insulin"), c("Fasting_Insulin"))
out_name = list(c("bmi_ukbgiant", "whradjbmi"), c("Volume_VAT", "Volume_ASAT"))
list_name <- c("Figure4Combined")
header <- c("Statistically informed", "Biologically informed")
xlab = paste0("Effect of 1-SD increase in fasting insulin on adiposity indices")
  doA<-res_steiger[exposure %in% unlist(exp_name) & outcome %in% unlist(out_name)]
  doA[,instrument_selection := "Statistically informed"]
  doB<-resmap[exposure %in% unlist(exp_name) & outcome %in% unlist(out_name)]
  doB[,instrument_selection := "Biologically informed"]
  doA <- rbindlist(list(doA, doB),fill = TRUE)
  doA <- doA[method != "MR Egger",]
  doA[,name := outcome %>% gsub("whradhbmi", "WHRadjBMI", .) %>% gsub("bmi_ukbgiant", "BMI", .) %>% gsub("_", " ", .) ]
doA[,outcome := factor(outcome, levels = unlist(out_name))]
doA[,method := factor(method, levels = c("Inverse variance weighted", "Weighted median","MR-PRESSO (outlier-corrected)" ,"IVW radial", "Contamination mixture"))]

  doB<- doA[instrument_selection ==header[2] ,]
  doA<- doA[instrument_selection ==header[1] ,]
  tobind <- doA[method %in% c("MR-PRESSO (outlier-corrected)"), ]
  tobind[,c("nsnp", "b", "se", "pval", "lci", "uci") := NA]
  tobind[,instrument_selection := "Biologically informed"]
  doB<- rbind(doB,tobind)
  doB <- doB[method != "IVW radial",]
  doA <- doA[order(outcome, method)]
  doB <- doB[order(outcome, method)]
  resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                         estimate = round(doA$b, digits =2),
                         lci =  round(doA$lci, digits = 2),
                         uci =  round(doA$uci, digits = 2),
                         n = doA$nsnp,
                         P_value = formatC(doA$pval, format = "e", digits = 1))
  
  resultsB <- data.frame(variable = as.character(1:nrow(doB)),
                         estimate = round(doB$b, digits =2),
                         lci =  round(doB$lci, digits = 2),
                         uci =  round(doB$uci, digits = 2),
                         n = doB$nsnp,
                         P_value = formatC(doB$pval, format = "e", digits = 1))
  
  mylabels <- data.frame(heading1 = as.character(doA$name),
                         heading2 = as.character(doA$method),
                         heading3 = as.character(NA),
                         variable = as.character(1:nrow(doA)))
  
  
  a<- make_forest_plot(panels = list(resultsA, resultsB),
                       col.key = "variable",
                       row.labels = mylabels,
                       exponentiate = FALSE,
                       pointsize = 2, 
                       rows = unique(mylabels$heading1),
                       col.stderr = NULL,
                       col.lci = "lci",
                       col.uci = "uci",
                       col.left         = c("n"),
                       col.left.heading = c("n SNP"),
                       col.right = "P_value",
                       col.right.heading = c("Effect (95% CI)", "P-value"),
                       xlab = xlab,
                       blankrows = c(0,1,0,0),
                       col.right.hjust = 1,
                       nullval = 0,
                       panel.headings = header)
  a
  ggsave(paste0("Results/", list_name, ".png"),
         width=864/72,height=408/72, units="in", scale=1,
         device = "png")
  
#Figure 6
res_ir[,align:="insuline_resistance_SNPs"]
res_is[,align:="other_SNPs"]
harm_ir[,align:="insuline_resistance_SNPs"]
harm_is[,align:="other_SNPs"]
harm_tot<- rbindlist(list(harm_ir, harm_is), fill = TRUE)
harm_tot[,Locus:=SNP]
res_tot <- rbindlist(list(res_is, res_ir), fill = TRUE)

mr_results <- res_tot
dat<-harm_tot
source("Analysis/my_mr_scatter_plot.R")

k <- my_mr_scatter_plot( dat = harm_tot, mr_results = res_tot, equation_facet_grid = " ~ align",legend.position = "top")
k + xlab("SNP effect on fasting insulin") +ylab("SNP effect on body mass index")

ggsave(paste0("Results/", "Figure7", ".png"),
       width=557/72,height=385/72, units="in", scale=1,
       device = "png")


#Figure 6 
resmap <- fread( "Data/Modified/resmap.txt")
resharm <- fread( "Data/Modified/harmmap_vinuela.txt")

at <- resmap[outcome %in% c("Volume_ASAT", "Volume_VAT"), ]
at <- at[!(method %in% c("Weighted mode","Robust adjusted profile score (RAPS)", "MR Egger")), ] 
at[, name := gsub("_", " ", outcome)]
doA <- at
resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))


mylabels <- data.frame(heading1 = doA$name,
                       heading2 = doA$method,
                       heading3 = as.character(NA),
                       variable = as.character(1:nrow(doA)))

ckbplotr::make_forest_plot(panels = list(resultsA),
                           col.key = "variable",
                           row.labels = mylabels,
                           exponentiate = TRUE,
                           pointsize = 2, 
                           rows = unique(mylabels$heading1),
                           col.stderr = NULL,
                           col.lci = "lci",
                           col.uci = "uci",
                           col.left         = c("n"),
                           col.left.heading = c("n SNP"),
                           col.right = "P_value",
                           col.right.heading = c("Effect (95% CI)", "P-value"),
                           xlab = c("Effect of 1-SD increase in FI on VAT/ASAT (SD)"),
                           blankrows = c(0,0,0,0),
                           col.right.hjust = 1,
                           panel.headings = NULL,
                           scalepoints = FALSE)

ggsave(paste0("Results/", "Figure3Clean", ".png"),
       width=483/72,height=250/72, units="in", scale=1,
       device = "png")

