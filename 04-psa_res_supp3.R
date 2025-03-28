# Appendix 9 plots 

CURRENT_RESDIR = "~/IncipienTB/PSA_RNAgam"

# =================================================================================================
# Input data preparation
# =================================================================================================

prod_by_coo_all_array = readRDS(paste0(CURRENT_RESDIR,"/prod_by_coo_all_array.RDS"))
eff_by_coo_all_array = readRDS(paste0(CURRENT_RESDIR,"/eff_by_coo_all.RDS"))


## subgroup
pop_size = cohort2019_by_coo2$total_pop[order(cohort2019_by_coo2$place_of_birth)]
v_risk1 = cohort2019_by_coo2$place_of_birth %in% riskcat_I
prod_by_strat_mean_risk1 = apply(prod_by_coo_all_array[v_risk1,,], c(2:3), function(x) weighted.mean(x, pop_size[v_risk1]/sum(pop_size[v_risk1]))) 
eff_by_strat_mean_risk1 = apply(eff_by_coo_all_array[v_risk1,,], c(2:3), function(x) weighted.mean(x, pop_size[v_risk1]/sum(pop_size[v_risk1])))

v_risk4 = cohort2019_by_coo2$place_of_birth %in% riskcat_IV
prod_by_strat_mean_risk4 = apply(prod_by_coo_all_array[v_risk4,,], c(2:3), function(x) weighted.mean(x, pop_size[v_risk4]/sum(pop_size[v_risk4]))) 
eff_by_strat_mean_risk4 = apply(eff_by_coo_all_array[v_risk4,,], c(2:3), function(x) weighted.mean(x, pop_size[v_risk4]/sum(pop_size[v_risk4])))

# =================================================================================================
# Useful functions
# =================================================================================================

get_inmbs4_s2_by_riskcost = function(nmb_rdsFileName, prod_df, eff_df, WTP, HCpersp){
  
  dat = readRDS(paste0(CURRENT_RESDIR,"/",nmb_rdsFileName))
  
  if(HCpersp){
    dat2 = lapply(1:1000, function(x) rowSums(dat[[x]][,2:3]) - prod_df[,x])
    dat3 = lapply(1:1000, function(x) eff_df[,x]*WTP - dat2[[x]])
    out = sapply(1:1000, function(x) dat3[[x]][3] - dat3[[x]][1]) # Strategy IV - Strategy II
  } else {
    dat2 = lapply(1:1000, function(x) rowSums(dat[[x]][,2:5]) - - prod_df[,x])
    dat3 = lapply(1:1000, function(x) eff_df[,x]*WTP - dat2[[x]])
    out = sapply(1:1000, function(x) dat3[[x]][3] - dat3[[x]][1]) # Strategy IV - Strategy II
  }
  
  return(out)
} 


get_plot_diff_df = function(nmb_rdsFileName, prod_df, eff_df, WTP, HCpersp){
  
  diff_s4_s2 = get_inmbs4_s2_by_riskcost(nmb_rdsFileName, prod_df, eff_df, WTP, HCpersp)
  
  df = data.frame(matrix(NA, ncol = 4, nrow = 1000))
  colnames(df) = c("diff", "YEARS", "SENS1", "SPEC")
  
  df$diff = diff_s4_s2
  
  PSAtable = read.csv("Data/PSAtable1000_gam.csv") # draft5
  df$YEARS = PSAtable$RNA_TIME_upper
  df$SENS1 = PSAtable$RNA_SENS1
  df$SPEC  = PSAtable$RNA_SPEC
  
  return(df)
}


plot_panel = function(slice, dat) {
  
  rsm = rsm::rsm(diff ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  test = contour(rsm, ~ TWI(SENS1,SPEC), at = data.frame(YEARS=slice))
  
  x = test$`SENS1 ~ SPEC`$x
  y = test$`SENS1 ~ SPEC`$y
  z = test$`SENS1 ~ SPEC`$z
  
  
  # Convert data to a data frame
  df <- data.frame(sens = rep(x, length(y)), spec = rep(y, each = length(x)), inmb = as.vector(z))
  
  # Create the plot
  p = ggplot(df, aes(x = sens, y = spec, z = inmb)) +
    stat_contour_filled(breaks= seq(-4000, 1200, 50), show.legend = F) + # c(min(df$inmb) + (0:10)*(900))) + 
    geom_contour(colour = "black") + 
    metR::geom_text_contour(aes(z = inmb), stroke = 0.2)  + 
    scale_fill_viridis_d(drop = FALSE, option = "plasma") + 
    xlab(paste0("Timeframe = ",slice," yrs")) + 
    ylab(" ") + 
    theme_minimal() + 
    theme(axis.title.x = element_text(margin = margin(t = 20, b = 0)),
          text = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(vjust = 0.1),
          axis.text.y = element_text(hjust = -2))
}

plot_fig = function(dat, title){
    
  p3 = plot_panel(slice = 3, dat)
  p4 = plot_panel(slice = 4, dat)
  p5 = plot_panel(slice = 5, dat)
  p6 = plot_panel(slice = 6, dat)
  p7 = plot_panel(slice = 7, dat)
  p8 = plot_panel(slice = 8, dat)
  p9 = plot_panel(slice = 9, dat)
  p10 = plot_panel(slice = 10, dat)
  
  # different color scheme: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
  
  option1 = ggpubr::ggarrange(p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 4) #plasma
  # option2 = ggpubr::ggarrange(p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 4) #inferno
  # option3 = ggpubr::ggarrange(p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 4) #direction-1 <- no good
  # option4 = ggpubr::ggarrange(p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 4) #plasma + direction-1 <-no good
  # option5 = ggpubr::ggarrange(p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 4) #turbo
  
  option1 
  
  
  library(ggpubr)
  option1_annotated = annotate_figure(option1, 
                                      left = text_grob("P(T- | no TB in lifeimte)", face = "bold", size = 16, rot = 90),
                                      bottom = text_grob("P(T+ | incipient TB)", face = "bold", size = 16, vjust = 1),
                                      top = text_grob(title,face = "bold", size = 18, vjust = -0.5))+ 
    theme(plot.margin = margin(15, 5, 15, 5))
  option1_annotated

return(option1_annotated)

}


# =================================================================================================
# Gather CEA results
# =================================================================================================

# # read in cost data for 1000 param sets
# cost_risk1_300 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_dft5.RDS"))
# cost_risk1_15 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_costItems15_dft5.RDS"))
# cost_risk4_300 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_dft5.RDS"))
# cost_risk4_15 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_costItems15_dft5.RDS"))

# =================================================================================================
# Prepare input data for plotting
# =================================================================================================

# calculate iNMB
ct_df_risk1_s4_hc_300 = get_plot_diff_df(nmb_rdsFileName = "/nmb_by_coo_itm_risk1_dft5.RDS",
                                         prod_df = prod_by_strat_mean_risk1,
                                         eff_df = eff_by_strat_mean_risk1,
                                         WTP = 150000,
                                         HCpersp = T)

plot_fig(ct_df_risk1_s4_hc_300, "Risk Category: > 300, Cost of HrTS: $300") #A9_Fig2-1-2_300

ct_df_risk1_s4_soc_300 = get_plot_diff_df(nmb_rdsFileName = "/nmb_by_coo_itm_risk1_dft5.RDS",
                                          prod_df = prod_by_strat_mean_risk1,
                                          eff_df = eff_by_strat_mean_risk1,
                                          WTP = 150000,
                                          HCpersp = F)
plot_fig(ct_df_risk1_s4_soc_300, "Risk Category: > 300, Cost of HrTS: $300") #A9_Fig2-1-1_300 #1200x680

ct_df_risk1_s4_hc_15 = get_plot_diff_df(nmb_rdsFileName = "nmb_by_coo_itm_risk1_costItems15_dft5.RDS", 
                                        prod_df = prod_by_strat_mean_risk1,
                                        eff_df = eff_by_strat_mean_risk1,
                                        WTP = 150000, 
                                        HCpersp = T)
plot_fig(ct_df_risk1_s4_hc_15, "Risk Category: > 300, Cost of HrTS: $15") #A9_Fig2-1-2_15

ct_df_risk1_s4_soc_15 = get_plot_diff_df(nmb_rdsFileName = "nmb_by_coo_itm_risk1_costItems15_dft5.RDS", 
                                         prod_df = prod_by_strat_mean_risk1,
                                         eff_df = eff_by_strat_mean_risk1,
                                         WTP = 150000, 
                                         HCpersp = F)
plot_fig(ct_df_risk1_s4_soc_15, "Risk Category: > 300, Cost of HrTA: $15") #A9_Fig2-1-1_15 

# calculate iNMB 
ct_df_risk4_s4_hc_300 = get_plot_diff_df(nmb_rdsFileName = "/nmb_by_coo_itm_risk4_dft5.RDS",
                                         prod_df = prod_by_strat_mean_risk4,
                                         eff_df = eff_by_strat_mean_risk4,
                                         WTP = 150000,
                                         HCpersp = T)
plot_fig(ct_df_risk4_s4_hc_300, "Risk Category: 0-10, Cost of HrTS: $300")

ct_df_risk4_s4_soc_300 = get_plot_diff_df(nmb_rdsFileName = "/nmb_by_coo_itm_risk4_dft5.RDS",
                                          prod_df = prod_by_strat_mean_risk4,
                                          eff_df = eff_by_strat_mean_risk4,
                                          WTP = 150000,
                                          HCpersp = F)
plot_fig(ct_df_risk4_s4_soc_300, "Risk Category: 0-10, Cost of HrTS: $300")

ct_df_risk4_s4_hc_15 = get_plot_diff_df(nmb_rdsFileName = "/nmb_by_coo_itm_risk4_dft5.RDS",
                                        prod_df = prod_by_strat_mean_risk4,
                                        eff_df = eff_by_strat_mean_risk4,
                                        WTP = 150000,
                                        HCpersp = T)
plot_fig(ct_df_risk4_s4_hc_15, "Risk Category: 0-10, Cost of HrTS: $15")

ct_df_risk4_s4_soc_15 = get_plot_diff_df(nmb_rdsFileName = "/nmb_by_coo_itm_risk4_dft5.RDS",
                                         prod_df = prod_by_strat_mean_risk4,
                                         eff_df = eff_by_strat_mean_risk4,
                                         WTP = 150000,
                                         HCpersp = F)
plot_fig(ct_df_risk4_s4_soc_15, "Risk Category: 0-10, Cost of RNA: $15")


