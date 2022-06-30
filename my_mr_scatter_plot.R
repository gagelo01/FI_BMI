my_mr_scatter_plot <- function (mr_results, dat, 
                                equation_facet_grid, legend.position = c(.32,.89) ) { #equation to parse into facet grid example includes "outcome ~  align", " ~  align", and "", if we do not want to have a grid
  
  d <- dat 
  
  d <- plyr::mutate(d)
  if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
    return(blank_plot("Insufficient number of SNPs"))
  }
  d <- subset(d, mr_keep)
  index <- d$beta.exposure < 0
  d$beta.exposure[index] <- d$beta.exposure[index] * 
    -1
  d$beta.outcome[index] <- d$beta.outcome[index] * 
    -1
  mrres <- mr_results
  mrres$a <- 0
  if ("MR Egger" %in% mrres$method) {
    d_split <- split(d, d$align)
    k_intercept <- lapply(d_split, function(data){
      temp <- TwoSampleMR::mr_egger_regression(data$beta.exposure, 
                                               data$beta.outcome, data$se.exposure, data$se.outcome, 
                                               default_parameters())
      return(data.table(align = unique(data$align), intercept = temp$b_i))}) %>%
      rbindlist(., fill = TRUE)
    mrres <- merge(mrres, k_intercept, by = "align")
    mrres[method == "MR Egger",a:=intercept];mrres[,intercept := NULL]
  }
  if ("MR Egger (bootstrap)" %in% mrres$method) {
    temp <- TwoSampleMR::mr_egger_regression_bootstrap(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
    mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
  }
  
  
  ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                         y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - 
                                                                                                    se.outcome, ymax = beta.outcome + se.outcome), 
                                                                                     colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - 
                                                                                                                                                          se.exposure, xmax = beta.exposure + se.exposure), 
                                                                                                                                           colour = "grey", height = 0) + ggplot2::geom_point() + 
    ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                    slope = b, colour = method), show.legend = TRUE) + 
    # ggplot2::scale_colour_manual(values = c("#a6cee3", 
    #                                         "#1f78b4", "#b2df8a", "#6a3d9a","#33a02c", "#fb9a99", 
    #                                         "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
    #                 "#b15928")) +
    ggplot2::scale_colour_manual(values = c(rgb(112, 54, 153, maxColorValue = 255),
                                            rgb(66, 94, 191, maxColorValue = 255), 
                                            rgb(84,201,237, maxColorValue = 255), 
                                            rgb(59,166,18, maxColorValue = 255), 
                                            rgb(255,110,26, maxColorValue = 255),
                                            rgb(149,199,71, maxColorValue = 255), 
                                            rgb(161,15,125, maxColorValue = 255), 
                                            rgb(249, 106, 27, maxColorValue = 255), 
                                            rgb(214,15,102, maxColorValue = 255),
                                            rgb(8, 161, 217, maxColorValue = 255),
                                            rgb(255,186,51, maxColorValue = 255), 
                                            rgb(54, 150, 214, maxColorValue = 255))) +
   ggplot2::labs(colour = "MR Test", 
                                                                        x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", 
                                                                                                                             d$outcome[1])) + ggplot2::theme(legend.position = "top", 
                                                                                                                                                             legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
    expand_limits(x = 0, y = 0) +
    ggrepel::geom_text_repel(label=d$Locus, check_overlap = T) + 
    theme(
      legend.background = element_rect(fill = "white", size = 4, colour = "white"),
      axis.ticks = element_line(colour = "black", size = 0.2),
      # panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.key=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      # panel.border = element_rect(colour = "grey70", fill=NA, size=1),
        axis.line.x.bottom = element_line(color = "black", size = 0.2),
        axis.line.y.left   = element_line(color = "black", size = 0.2),
        axis.line.y.right  = element_blank(),
        axis.text.y.right  = element_blank(),
        panel.border       = element_blank(),
      legend.position=legend.position) + 
    facet_grid(eval(parse(text = equation_facet_grid))) +
    expand_limits(x = 0, y = 0)
}
