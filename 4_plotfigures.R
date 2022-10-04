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
mvmr_results <- fread("Results/resmvmr.txt")
harm_univariate <- fread( "Data/Modified/harm_univariate.txt")
harm_univariate[, exposure_outcome :=paste0(exposure, "_", outcome)]
harm_steiger <- harm_univariate[!(steiger_dir == FALSE & steiger_pval < 0.05),]
##########
fileformat <- "tiff"
dpi <- 300
########
tochange<- c("res_steiger")
for(i in 1:length(tochange)) {
  assign(tochange[i], get(tochange[i])[!(method %in% c("Robust adjusted profile score (RAPS)", "Weighted mode")), ])
}
##Figure 1 and 4

res_univariate <- res_steiger
mr_results<-res_univariate
mr_results[,align:=ifelse(grepl("_bio$", exposure), "Biologically driven", "Statistically driven")]
mr_results[,exposure:=gsub("_bio$","", exposure)]
mr_results[,name := outcome %>% gsub("whradhbmi", "WHRadjBMI", .) %>% gsub("bmi_ukbgiant", "BMI", .) %>% gsub("_", " ", .) ]
mr_results[,method := factor(method, levels = c("Inverse variance weighted", "Weighted median", "MR Egger", "MR-PRESSO (outlier-corrected)" ,"IVW radial", "Contamination mixture"))]
mr_results[,align := factor(align, levels = c("Statistically driven", "Biologically driven"))]
mr_results[method ==  "Contamination mixture", se := mean(c(b-lci, uci-b))/1.96]
mr_results <- mr_results[order(align, method),]
#Figure 1
i<-1
out_name = c("Fasting_Insulin")
exp_name = c("bmi_ukbgiant")
doA <- mr_results[exposure == exp_name[i] & outcome == out_name[1], ]

resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))


mylabels <- data.frame(heading1 = as.character(doA$align),
                       heading2 = as.character(doA$method),
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
                 xlab = c("Effect of one SD increase in BMI on fasting insulin (SD)"),
                 blankrows = c(0,1,0,0),
                 col.right.hjust = 1,
                 nullval = 0,
                 panel.headings = NULL,
                 xlim = c(0, 1))

ggsave(paste0("Results/", "Figure2", ".", fileformat),
       width=550/72,height=328/72, units="in", scale=1,
       device = fileformat, dpi = dpi)
###plotting Figure 2 and supplementary figure 1
exp_name = c("whradjbmi")
out_name = c("Fasting_Insulin")
file_name = c("Figure3")
panel_heading <- ""
i<-1
doA <- mr_results[exposure %in% exp_name[[i]] & outcome %in% out_name[i]]

resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                       estimate = round(doA$b, digits =2),
                       lci =  round(doA$lci, digits = 2),
                       uci =  round(doA$uci, digits = 2),
                       n = doA$nsnp,
                       P_value = formatC(doA$pval, format = "e", digits = 1))


mylabels <- data.frame(heading1 = as.character(doA$method),
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
                 xlab = c("Effect of one SD Increase of WHRadjBMI on Fasting Insulin (SD)"),
                 blankrows = c(0,1,0,0),
                 col.right.hjust = 1,
                 nullval = 0,
                 xlim = c(0,0.4))

ggsave(paste0("Results/", file_name[i], ".", fileformat ),
       width=541/72,height=200/72, units="in", scale=1,
       device = fileformat, dpi = dpi)

#Figure 3
unimeth<-"Inverse variance weighted" 
multimeth<- c("Multivariable IVW", "Multivariable Median",
              "Multivariable Lasso", "Multivariable Egger")             
data <- mvmr_results[otherexposure =="UKB-b-9405" | exposure == "UKB-b-9405" & outcome == "Fasting_Insulin",]
data<-rbind(data, mr_results[exposure %in% c("UKB-b-9405", "bmi_ukbgiant") & outcome == "Fasting_Insulin" & align == "Statistically driven",], fill = TRUE)
data<-data[method %in%c(unimeth,multimeth),]
data[, Category_other := ifelse(method %in% unimeth, "Univariable", "Multivariable")]
data[, Category_other := factor(Category_other, levels = c("Univariable", "Multivariable"))]
data[,name := gsub("bmi_ukbgiant", "Body mass index", exposure) %>% gsub("UKB-b-9405", "Waist circumference", .)]

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

ggsave(paste0("Results/", "Figure4", ".", fileformat),
       width=350/72,height=250/72, units="in", scale=1,
       device = fileformat, dpi = dpi)

#Figure 4
exp_name = list(c("Fasting_Insulin"), c("Fasting_Insulin"))
out_name = c("mri_vat", "mri_asat", "mri_gfat", "mri_vatAsatRatio",  "mri_vatGfatRatio", "whradjbmi", "bmi_ukbgiant")
xlab = paste0("Effect of one SD increase in fasting insulin on adiposity indices")

forfig4<-mr_results[method != "MR Egger" & outcome %in% out_name, ]
forfig4[,outcome := factor(outcome, levels = out_name)]

k <- forestplot(
  df = forfig4,
  name = name,
  se = se,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = xlab,
  logodds = FALSE,
  ci = 0.95,
  colour = method,
  # xlim = data[, round(c(min(exp(b))-0.1,max(exp(b))+0.1), digits = 1)]
  # xlim = data[, round(c(min(exp(lci)),max(exp(uci))), digits = 1)]
)

k<-k+ facet_wrap(facets =  ~ align)
k
ggsave(paste0("Results/", "Figure5", ".", fileformat), plot = k,
       width=864/72,height=408/72, units="in", scale=1,
       device = fileformat, dpi = dpi)

#Supplementary_Figure 2
source("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI/Analysis/my_mr_scatter_plot.R")
dat<- harm_steiger
forsupfig2<-mr_results[method != "MR Egger" & outcome %in% out_name, ]
forsupfig2[,outcome := factor(outcome, levels = out_name)]

dat[,align := ifelse(grepl("_bio$", exposure), "Biologically driven", "Statistically driven")]
dat[,exposure:=gsub("_bio$","", exposure)]
dat[, Locus := ""]
dat[, align := factor(align, levels = c("Statistically driven","Biologically driven"))]
dat<- dat[exposure == "Fasting_Insulin" & outcome %in% out_name, ]
dat[, outcome := factor(outcome, levels = out_name)]
dat <- dat[order(align,outcome),]

k <- my_mr_scatter_plot( dat = dat, mr_results = forsupfig2, equation_facet_grid = "outcome ~  align", legend.position = "top")
k + xlab("SNP effect on FI") +ylab("SNP effect on adiposity traits")

ggsave(paste0("Results/", "Supplementary_Figure2", ".", fileformat),
       width=435/72,height=758/72, units="in", scale=1, dpi = dpi,
       device = fileformat)
