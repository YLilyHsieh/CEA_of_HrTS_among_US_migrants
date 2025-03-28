
# ==============================================================================
# APPENDIX 1 TABLE 1 Characteristics of study cohort by country-of-origin
# ==============================================================================

cohort2019_by_coo3 = cohort2019_by_coo2[order(-cohort2019_by_coo2$ltbi_calib),]

S1Table1 = cohort2019_by_coo3 %>% 
  select(c("place_of_birth", "total_pop", "mean_age", "ltbi_calib", "ltbi_calib_lower", "ltbi_calib_upper")) %>% 
  mutate(across(3:6, round, 2)) 

ltbi = lapply(1:nrow(S1Table1), function(x) pretty(c(S1Table1$ltbi_calib[x], S1Table1$ltbi_calib_lower[x], S1Table1$ltbi_calib_upper[x]), 2))
ltbi = do.call("rbind", ltbi)

S1Table1 = cbind(S1Table1[,c("place_of_birth", "total_pop", "mean_age")], ltbi)
colnames(S1Table1) = c("Countries/regions", "Modelled population size", "Mean entry age", "LTBI prevalence, % (95% CI)")

write.table(S1Table1, file = "S1Table.txt", sep = ",", quote = FALSE, row.names = F)

# ==============================================================================
# APPENDIX 6 TABLE 1 incremental gain in QALY by risk category
# ==============================================================================


s1 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(eff_by_riskcat$diff.IGRA_TB[x], 
                                                         eff_by_riskcat$diff.IGRA_TB_lower[x],
                                                         eff_by_riskcat$diff.IGRA_TB_upper[x]), 4))
s1 = do.call("rbind", s1)

s2 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(eff_by_riskcat$diff.IGRA_RNA_TB[x], 
                                                         eff_by_riskcat$diff.IGRA_RNA_TB_lower[x],
                                                         eff_by_riskcat$diff.IGRA_RNA_TB_upper[x]), 4))
s2 = do.call("rbind", s2)

s3 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(eff_by_riskcat$diff.RNA_TB[x], 
                                                         eff_by_riskcat$diff.RNA_TB_lower[x],
                                                         eff_by_riskcat$diff.RNA_TB_upper[x]), 4))
s3 = do.call("rbind", s3)

S6Table1 = cbind(s1, s2, s3)
colnames(S6Table1) = c("IGRA_TB", "IGRA_RNA_TB", "RNA_TB")

write.table(S6Table1, file = "S6Table1.txt", sep = ",", quote = FALSE, row.names = F)


# ==============================================================================
# APPENDIX 6 TABLE 2 incremental cost by risk category
# ==============================================================================

s1 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(cost_by_riskcat$diff.IGRA_TB[x], 
                                                         cost_by_riskcat$diff.IGRA_TB_lower[x],
                                                         cost_by_riskcat$diff.IGRA_TB_upper[x]), 2))
s1 = do.call("rbind", s1)

s2 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(cost_by_riskcat$diff.IGRA_RNA_TB[x], 
                                                         cost_by_riskcat$diff.IGRA_RNA_TB_lower[x],
                                                         cost_by_riskcat$diff.IGRA_RNA_TB_upper[x]), 2))
s2 = do.call("rbind", s2)

s3 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(cost_by_riskcat$diff.RNA_TB[x], 
                                                         cost_by_riskcat$diff.RNA_TB_lower[x],
                                                         cost_by_riskcat$diff.RNA_TB_upper[x]), 2))
s3 = do.call("rbind", s3)

S6Table2 = cbind(s1, s2, s3)
colnames(S6Table2) = c("IGRA_TB", "IGRA_RNA_TB", "RNA_TB")

write.table(S6Table2, file = "S6Table2.txt", sep = ",", quote = FALSE, row.names = F)

# ==============================================================================
# APPENDIX 6 TABLE 3 iNMB by risk category
# ==============================================================================

s1 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(nmb_by_riskcat$iNMB.IGRA_TB[x], 
                                                         nmb_by_riskcat$iNMB.IGRA_TB_lower[x],
                                                         nmb_by_riskcat$iNMB.IGRA_TB_upper[x]), 2))
s1 = do.call("rbind", s1)

s2 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(nmb_by_riskcat$iNMB.IGRA_RNA_TB[x], 
                                                         nmb_by_riskcat$iNMB.IGRA_RNA_TB_lower[x],
                                                         nmb_by_riskcat$iNMB.IGRA_RNA_TB_upper[x]), 2))
s2 = do.call("rbind", s2)

s3 = lapply(1:nrow(nmb_by_riskcat), function(x) pretty(c(nmb_by_riskcat$iNMB.RNA_TB[x], 
                                                         nmb_by_riskcat$iNMB.RNA_TB_lower[x],
                                                         nmb_by_riskcat$iNMB.RNA_TB_upper[x]), 2))
s3 = do.call("rbind", s3)

S6Table3 = cbind(s1, s2, s3)
colnames(S6Table3) = c("IGRA_TB", "IGRA_RNA_TB", "RNA_TB")

write.table(S6Table3, file = "S6Table3.txt", sep = ",", quote = FALSE, row.names = F)

# ==============================================================================
# SUPP RESULTS cost figure for different RNA costs
# ==============================================================================

library('ggpubr')

# Useful functions
get_plot_df.cost = function(newlist){
  
  nmb_by_coo_array2 = array(unlist(newlist), dim = c(97, 53, 1000))
  
  
  nmb_by_coo_mean = apply(nmb_by_coo_array2, c(1:2), mean) # dataframe 97 x 36
  colnames(nmb_by_coo_mean) = colnames(newlist[[1]])
  cost_by_coo_mean = nmb_by_coo_mean %>% 
    as.data.frame() %>%
    select(diff.IGRA_TB, diff.IGRA_RNA_TB, diff.RNA_TB)
  
  
  nmb_by_coo_upper = apply(nmb_by_coo_array2, c(1:2), quantile, p = 0.975) # dataframe 97 x 36
  colnames(nmb_by_coo_upper) = colnames(newlist[[1]])
  cost_by_coo_upper = nmb_by_coo_upper %>% 
    as.data.frame() %>%
    select(diff.IGRA_TB, diff.IGRA_RNA_TB, diff.RNA_TB) %>% 
    rename(., diff.IGRA_TB_upper = diff.IGRA_TB) %>%
    rename(., diff.IGRA_RNA_TB_upper = diff.IGRA_RNA_TB) %>% 
    rename(., diff.RNA_TB_upper = diff.RNA_TB) %>%
    select(diff.IGRA_TB_upper, diff.IGRA_RNA_TB_upper, diff.RNA_TB_upper)
  
  nmb_by_coo_lower = apply(nmb_by_coo_array2, c(1:2), quantile, p = 0.025) # dataframe 97 x 36
  colnames(nmb_by_coo_lower) = colnames(newlist[[1]])
  cost_by_coo_lower = nmb_by_coo_lower %>% 
    as.data.frame() %>%
    select(diff.IGRA_TB, diff.IGRA_RNA_TB, diff.RNA_TB) %>% 
    rename(., diff.IGRA_TB_lower = diff.IGRA_TB) %>%
    rename(., diff.IGRA_RNA_TB_lower = diff.IGRA_RNA_TB) %>% 
    rename(., diff.RNA_TB_lower = diff.RNA_TB) %>%
    select(diff.IGRA_TB_lower, diff.IGRA_RNA_TB_lower, diff.RNA_TB_lower)
  
  
  cost_by_coo_mean$country = cohort2019_by_coo2$place_of_birth
  cost_by_coo_mean$tb_inc = cohort2019_by_coo2$tb_inc
  
  
  cost_by_coo = cbind(cost_by_coo_mean, cost_by_coo_lower, cost_by_coo_upper)
  
  
  coo_by_TBinc_temp = coo_by_TBinc %>% rename(., country = place_of_birth)
  
  temp3 = merge( coo_by_TBinc_temp, cost_by_coo ,by = "country")
  
  
  cost_by_riskcat = temp3 %>% group_by(risk_category) %>%
    summarise(diff.IGRA_TB = weighted.mean(diff.IGRA_TB, total_pop),
              diff.IGRA_RNA_TB = weighted.mean(diff.IGRA_RNA_TB, total_pop),
              diff.RNA_TB = weighted.mean(diff.RNA_TB, total_pop),
              diff.IGRA_TB_lower = weighted.mean(diff.IGRA_TB_lower, total_pop),
              diff.IGRA_RNA_TB_lower = weighted.mean(diff.IGRA_RNA_TB_lower, total_pop),
              diff.RNA_TB_lower = weighted.mean(diff.RNA_TB_lower, total_pop),
              diff.IGRA_TB_upper = weighted.mean(diff.IGRA_TB_upper, total_pop),
              diff.IGRA_RNA_TB_upper = weighted.mean(diff.IGRA_RNA_TB_upper, total_pop),
              diff.RNA_TB_upper = weighted.mean(diff.RNA_TB_upper, total_pop))
  
  cost_by_riskcat.mean = cost_by_riskcat %>% select(risk_category, diff.IGRA_TB, diff.IGRA_RNA_TB, diff.RNA_TB)
  cost_by_riskcat_reshape = reshape2::melt(cost_by_riskcat.mean, id.vars = c("risk_category")) %>% rename(., mean = value)
  
  cost_by_riskcat.upper = cost_by_riskcat %>% select(risk_category, diff.IGRA_TB_upper, diff.IGRA_RNA_TB_upper, diff.RNA_TB_upper)
  cost_by_riskcat.upper_reshape = reshape2::melt(cost_by_riskcat.upper, id.vars = c("risk_category")) %>% rename(., upper = value) %>% select(upper)
  
  
  cost_by_riskcat.lower = cost_by_riskcat %>% select(risk_category, diff.IGRA_TB_lower, diff.IGRA_RNA_TB_lower, diff.RNA_TB_lower)
  cost_by_riskcat.lower_reshape = reshape2::melt(cost_by_riskcat.lower, id.vars = c("risk_category")) %>% rename(., lower = value) %>% select(lower)
  
  cost_by_risk_plot = cbind(cost_by_riskcat_reshape, cost_by_riskcat.upper_reshape, cost_by_riskcat.lower_reshape)
  
  return(cost_by_risk_plot)
}

# RNA COST = 15, 30, 60, 150
nmb_by_coo_list2 = NULL

for (i in c(1:1000)){
  nmb_by_coo_list2[[i]] = nmb_by_coo_list[[i]] %>% 
    mutate(diff.IGRA_TB        =  (cost.IGRA_TB -cost.No_testing),
           diff.IGRA_RNA_TB    =  (cost15.IGRA_RNA_TB - cost.No_testing),
           diff.RNA_TB         =  (cost15.RNA_TB - cost.No_testing))
}

plot.df15 = get_plot_df.cost(newlist = nmb_by_coo_list2)

nmb_by_coo_list2 = NULL

for (i in c(1:1000)){
  nmb_by_coo_list2[[i]] = nmb_by_coo_list[[i]] %>% 
    mutate(diff.IGRA_TB        =  (cost.IGRA_TB -cost.No_testing),
           diff.IGRA_RNA_TB    =  (cost30.IGRA_RNA_TB - cost.No_testing),
           diff.RNA_TB         =  (cost30.RNA_TB - cost.No_testing))
}

plot.df30 = get_plot_df.cost(newlist = nmb_by_coo_list2)

nmb_by_coo_list2 = NULL

for (i in c(1:1000)){
  nmb_by_coo_list2[[i]] = nmb_by_coo_list[[i]] %>% 
    mutate(diff.IGRA_TB        =  (cost.IGRA_TB -cost.No_testing),
           diff.IGRA_RNA_TB    =  (cost60.IGRA_RNA_TB - cost.No_testing),
           diff.RNA_TB         =  (cost60.RNA_TB - cost.No_testing))
}

plot.df60 = get_plot_df.cost(newlist = nmb_by_coo_list2)

nmb_by_coo_list2 = NULL

for (i in c(1:1000)){
  nmb_by_coo_list2[[i]] = nmb_by_coo_list[[i]] %>% 
    mutate(diff.IGRA_TB        =  (cost.IGRA_TB -cost.No_testing),
           diff.IGRA_RNA_TB    =  (cost150.IGRA_RNA_TB - cost.No_testing),
           diff.RNA_TB         =  (cost150.RNA_TB - cost.No_testing))
}

plot.df150 = get_plot_df.cost(newlist = nmb_by_coo_list2)

nmb_by_coo_list2 = NULL

for (i in c(1:1000)){
  nmb_by_coo_list2[[i]] = nmb_by_coo_list[[i]] %>% 
    mutate(diff.IGRA_TB        =  (cost.IGRA_TB -cost.No_testing),
           diff.IGRA_RNA_TB    =  (cost.IGRA_RNA_TB - cost.No_testing),
           diff.RNA_TB         =  (cost.RNA_TB - cost.No_testing))
}

plot.df300 = get_plot_df.cost(newlist = nmb_by_coo_list2)

plot.df_full = rbind(plot.df300, plot.df150, plot.df60, plot.df30, plot.df15)
plot.df_full$cost = rep(c("$300", "$150", "$60", "$30", "15"), each = nrow(plot.df300))
plot.df_full$cost_f = factor(plot.df_full$cost, levels =  c("$300", "$150", "$60", "$30", "15"))


c = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df_full) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df_full, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  xlab("Risk category") +
  ylab("Per-person increase in cost (USD)") +
  # xlab(" ") +
  # ylab(" ") +
  scale_y_continuous(limits = c(0, 450)) +
  facet_wrap(~cost_f, scale = "free", nrow = 1)+
  ggtitle(expression(bold(atop("Increase in cost relative to Strategy I", "by cost of the signature (USD)"))))+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  theme_line() + 
  theme( 
    legend.position = "bottom",
    axis.text.x = element_text(angle =0), 
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    plot.margin = margin(10, 5.5, 0, 5.5),
    title =element_text(size=14, face='bold', margin = margin(b = 30))
  ) 



# ==============================================================================
# APPENDIX 7 Figure 1 iNMB figure for different RNA costs and WTP - rev 
# ==============================================================================

eff_by_coo_risk1 = readRDS(paste0(CURRENT_RESDIR,"/eff_by_coo_risk1.RDS"))
eff_by_coo_risk2 = readRDS(paste0(CURRENT_RESDIR,"/eff_by_coo_risk2.RDS"))
eff_by_coo_risk3 = readRDS(paste0(CURRENT_RESDIR,"/eff_by_coo_risk3.RDS"))
eff_by_coo_risk4 = readRDS(paste0(CURRENT_RESDIR,"/eff_by_coo_risk4.RDS"))

prod_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/prod_by_coo_risk1.RDS"))
prod_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/prod_by_coo_risk2.RDS"))
prod_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/prod_by_coo_risk3.RDS"))
prod_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/prod_by_coo_risk4.RDS"))


get_riskcat_nmb_mean = function(A, B, C, wtp){
  
  A_array = array(unlist(A), dim = c(3, 5, 1000))
  
  z = C[2:4,]
  
  x= apply(A_array, 3, function(x) rowSums(x[,2:5]))
  y = B[2:4,]
  meanSOC = round(rowMeans(z*wtp - (x - y)), digits = 1)
  # ciSOC = round(apply((z*150000 - (x - y)), 1, quantile, c(0.025, 0.975)), digits = 1)
  
  x2= apply(A_array, 3, function(x) rowSums(x[,2:3]))
  meanHC = round(rowMeans(z*wtp - x2), digits = 1)
  # ciHC = round(apply((z*150000 - x2), 1, quantile, c(0.025, 0.975)), digits=1)
  
  return(list = c(meanSOC = meanSOC, meanHC = meanHC))
}



get_inmb_byCost = function(cost, wtp){
  
  # if (cost == 30) {
  #   nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_", cost, ".RDS"))
  #   nmb_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_", cost, ".RDS"))
  #   nmb_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_", cost, ".RDS"))
  #   nmb_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_", cost, ".RDS"))
  # } else 
  if (cost == 300) {
    nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_costItems.RDS"))
    nmb_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_costItems.RDS"))
    nmb_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_costItems.RDS"))
    nmb_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_costItems.RDS"))
  } else {
    nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_costItems", cost, ".RDS"))
    nmb_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_costItems", cost, ".RDS"))
    nmb_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_costItems", cost, ".RDS"))
    nmb_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_costItems", cost, ".RDS"))
  } 
  
  risk1 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk1, B = prod_by_coo_itm_risk1, C = eff_by_coo_risk1, wtp = wtp)
  risk2 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk2, B = prod_by_coo_itm_risk2, C = eff_by_coo_risk2, wtp = wtp)
  risk3 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk3, B = prod_by_coo_itm_risk3, C = eff_by_coo_risk3, wtp = wtp)
  risk4 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk4, B = prod_by_coo_itm_risk4, C = eff_by_coo_risk4, wtp = wtp)
  
  wtp_text = format(wtp, big.mark=",", scientific = F)
  df = data.frame(mean = c(risk1, risk2, risk3, risk4),
                  RiskCat = rep(c("I", "II", "III", "IV"), each = 6),
                  Strategy = rep(rep(c("IGRA_TB", "IGRA_HrTS_TB", "HrTS_TB"), 2), 4), 
                  Perspective = rep(rep(c("SOC", "HC"), each = 3), 4),
                  Cost = rep(paste0("$",cost), 24),
                  WTP = rep(wtp_text, 24))
  
  
}


get_inmb_byWTP = function(wtp){
  
  cost300 = get_inmb_byCost(300, wtp)
  cost150 = get_inmb_byCost(150, wtp)
  cost60  = get_inmb_byCost(60, wtp)
  cost30  = get_inmb_byCost(30, wtp)
  cost15  = get_inmb_byCost(15, wtp)
  
  out = rbind(cost300, cost150, cost60, cost30, cost15)
  
  return(out)
  
}

plot_inmb_df = lapply(c(150000, 100000, 50000, 30000), function(x) get_inmb_byWTP(x))
plot_inmb_df2 = do.call("rbind", plot_inmb_df)

plot_inmb_df2$cost_f = factor(plot_inmb_df2$Cost, levels = c("$300", "$150", "$60", "$30", "$15"))
plot_inmb_df2$wtp_f  = factor(plot_inmb_df2$WTP,  levels = c("150,000", "100,000","50,000", "30,000"))
plot_inmb_df2$strategy_f = factor(plot_inmb_df2$Strategy, levels = c("IGRA_TB", "IGRA_HrTS_TB","HrTS_TB"))


max_y = 550
min_y = -400

ggplot(aes(y = mean, x= as.factor(RiskCat), fill = strategy_f), data = plot_inmb_df2[plot_inmb_df2$Perspective == "HC",]) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  # geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(RiskCat)), data = plot.df.inmb_full2, 
  #               width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  xlab("Risk Category") +
  ylab("Per-person iNMB (USD)") +
  # xlab(" ") +
  # ylab(" ") +
  # labs(title = "$300") +
  #ggtitle(expression(bold(atop("Inc. NMB relative to 'I. No-screening' in a Healthcare sector perspective", "by costs of HrTS and WTP thresholds"))))+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA only", "III: IGRA-HrTS", "IV: HrTS only"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "10-100", "0-10"))+
  scale_y_continuous(limits = c(min_y, max_y)) + 
  facet_grid(wtp_f ~ cost_f) + 
  theme_line() + 
  theme(
    legend.position = "bottom",
    #     legend.justification = c("right", "top"),
    # legend.box.just = "bottom",
    #     legend.margin = margin(6, 6, 6, 6),
    axis.text.x = element_text(angle = 300, hjust = 0, vjust = 0.5, size = 13.5),
    axis.text.y = element_text(size = 13.5),
    axis.title.x = element_text(margin = margin(t = 15), size = 14),
    axis.title.y = element_text(margin = margin(r = 15), size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(15, 15, 15, 15),
    title =element_text(size=14, face='bold')
  ) # 1150x 790


ggplot(aes(y = mean, x= as.factor(RiskCat), fill = strategy_f), data = plot_inmb_df2[plot_inmb_df2$Perspective == "SOC",]) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  # geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(RiskCat)), data = plot.df.inmb_full2, 
  #               width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  xlab("Risk Category") +
  ylab("Per-person iNMB (USD)") +
  # xlab(" ") +
  # ylab(" ") +
  # labs(title = "$300") +
  #ggtitle(expression(bold(atop("Inc. NMB relative to 'I. No-screening' in a Healthcare sector perspective", "by costs of HrTS and WTP thresholds"))))+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA only", "III: IGRA-HrTS", "IV: HrTS only"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "10-100", "0-10"))+
  scale_y_continuous(limits = c(min_y, max_y)) + 
  facet_grid(wtp_f ~ cost_f) + 
  theme_line() + 
  theme(
    legend.position = "bottom",
    #     legend.justification = c("right", "top"),
    # legend.box.just = "bottom",
    #     legend.margin = margin(6, 6, 6, 6),
    axis.text.x = element_text(angle = 300, hjust = 0, vjust = 0.5, size = 13.5),
    axis.text.y = element_text(size = 13.5),
    axis.title.x = element_text(margin = margin(t = 15), size = 14),
    axis.title.y = element_text(margin = margin(r = 15), size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(15, 15, 15, 15),
    title =element_text(size=14, face='bold')
  ) # 1075x735

# ==============================================================================
# APPENDIX 7 Figure 1 iNMB figure for different RNA costs and WTP
# ==============================================================================

# INCREMENTAL iNMB

# useful functions
get_plot_df.inmb = function(nmb_by_riskcat){
  
  nmb_by_riskcat2.mean = nmb_by_riskcat %>% select(risk_category, iNMB.IGRA_TB, iNMB.IGRA_RNA_TB, iNMB.RNA_TB)
  nmb_by_riskcat2.mean_reshape = reshape2::melt(nmb_by_riskcat2.mean, id.vars = c("risk_category")) %>% rename(., mean = value)
  
  nmb_by_riskcat2.lower = nmb_by_riskcat %>% select(risk_category, iNMB.IGRA_TB_lower, iNMB.IGRA_RNA_TB_lower, iNMB.RNA_TB_lower)
  nmb_by_riskcat2.lower_reshape = reshape2::melt(nmb_by_riskcat2.lower, id.vars = c("risk_category")) %>% rename(., lower = value)
  
  nmb_by_riskcat2.upper = nmb_by_riskcat %>% select(risk_category, iNMB.IGRA_TB_upper, iNMB.IGRA_RNA_TB_upper, iNMB.RNA_TB_upper)
  nmb_by_riskcat2.upper_reshape = reshape2::melt(nmb_by_riskcat2.upper, id.vars = c("risk_category")) %>% rename(., upper = value)
  
  nmb_byrisk_reshape = cbind(nmb_by_riskcat2.mean_reshape, nmb_by_riskcat2.lower_reshape[,"lower"], nmb_by_riskcat2.upper_reshape[,"upper"])
  colnames(nmb_byrisk_reshape) = c(colnames(nmb_by_riskcat2.mean_reshape), "lower", "upper")
  
  return(nmb_byrisk_reshape)
}


get_iNMB_byWTP = function(WTP){
# nmb_by_coo_mean2 = read.csv(paste0(CURRENT_RESDIR,"/nmb_by_coo_mean2.csv"))
# Calculate the lower and upper bounds 
nmb_by_coo_array = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_all_array.RDS"))


nmb_by_coo_list = future_lapply(1:1000, function(x) as.data.frame(nmb_by_coo_array[,,x]))  # list of 1000 dataframes 97 x 36


for (i in c(1:1000)){
  colnames(nmb_by_coo_list[[i]]) = c("eff.dQALY.No_testing", "eff.dQALY.IGRA_TB", "eff.dQALY.IGRA_RNA_TB", "eff.dQALY.RNA_TB",
                                     "eff.QALY.No_testing", "eff.QALY.IGRA_TB", "eff.QALY.IGRA_RNA_TB", "eff.QALY.RNA_TB",
                                     "eff.dLY.No_testing", "eff.dLY.IGRA_TB", "eff.dLY.IGRA_RNA_TB", "eff.dLY.RNA_TB",
                                     "eff.LY.No_testing", "eff.LY.IGRA_TB", "eff.LY.IGRA_RNA_TB", "eff.LY.RNA_TB",
                                     "cost.No_testing", "cost.IGRA_TB", "cost.IGRA_RNA_TB", "cost.RNA_TB",
                                     "cost150.No_testing", "cost150.IGRA_TB", "cost150.IGRA_RNA_TB", "cost150.RNA_TB",
                                     "cost60.No_testing", "cost60.IGRA_TB", "cost60.IGRA_RNA_TB", "cost60.RNA_TB",
                                     "cost30.No_testing", "cost30.IGRA_TB", "cost30.IGRA_RNA_TB", "cost30.RNA_TB",
                                     "cost15.No_testing", "cost15.IGRA_TB", "cost15.IGRA_RNA_TB", "cost15.RNA_TB")
}


for (i in c(1:1000)){
  nmb_by_coo_list[[i]] = nmb_by_coo_list[[i]] %>%
    mutate(NMB.No_screening = eff.dQALY.No_testing*WTP - cost.No_testing,
           NMB.IGRA_TB      = eff.dQALY.IGRA_TB*WTP - cost.IGRA_TB,
           NMB.IGRA_RNA_TB  = eff.dQALY.IGRA_RNA_TB*WTP - cost.IGRA_RNA_TB,
           iNMB.IGRA_TB        =  (eff.dQALY.IGRA_TB*WTP - cost.IGRA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB    =  (eff.dQALY.IGRA_RNA_TB*WTP - cost.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB150 =  (eff.dQALY.IGRA_RNA_TB*WTP - cost150.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB60  =  (eff.dQALY.IGRA_RNA_TB*WTP - cost60.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB30  =  (eff.dQALY.IGRA_RNA_TB*WTP - cost30.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB15  =  (eff.dQALY.IGRA_RNA_TB*WTP - cost15.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB    =  (eff.dQALY.RNA_TB*WTP - cost.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB150 =  (eff.dQALY.RNA_TB*WTP - cost150.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB60  =  (eff.dQALY.RNA_TB*WTP - cost60.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB30  =  (eff.dQALY.RNA_TB*WTP - cost30.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB15  =  (eff.dQALY.RNA_TB*WTP - cost15.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing))
}

nmb_by_coo_array2 = array(unlist(nmb_by_coo_list), dim = c(97, 50, 1000))
nmb_by_coo_mean = apply(nmb_by_coo_array2, c(1:2), mean) # dataframe 97 x 36

nmb_by_coo_mean2 = data.frame(country = countrynames,
                              tb_inc  = cohort2019_by_coo2$tb_inc[order(cohort2019_by_coo2$place_of_birth)],
                              pop_size = cohort2019_by_coo2$total_pop[order(cohort2019_by_coo2$place_of_birth)],
                              nmb_by_coo_mean)
colnames(nmb_by_coo_mean2) = c("country", "tb_inc", "pop_size", colnames(nmb_by_coo_list[[1]]))
# 
# write.csv(nmb_by_coo_mean2, paste0(CURRENT_RESDIR,"/nmb_by_coo_mean2.csv"), row.names = F)

# nmb_by_coo_mean2 = read.csv(paste0(CURRENT_RESDIR,"/nmb_by_coo_mean2.csv"))
# Calculate the lower and upper bounds 
nmb_by_coo_upper = apply(nmb_by_coo_array2, c(1:2), quantile, p = 0.975) # dataframe 97 x 36
colnames(nmb_by_coo_upper) = colnames(nmb_by_coo_list[[1]])
nmb_by_coo_upper = nmb_by_coo_upper %>% 
  as.data.frame() %>%
  select(c(iNMB.IGRA_TB, iNMB.IGRA_RNA_TB, iNMB.RNA_TB, 
           iNMB.IGRA_RNA_TB150, iNMB.IGRA_RNA_TB60, iNMB.IGRA_RNA_TB30, iNMB.IGRA_RNA_TB15,
           iNMB.RNA_TB150, iNMB.RNA_TB60, iNMB.RNA_TB30, iNMB.RNA_TB15)) %>% 
  rename(., iNMB.IGRA_TB_upper = iNMB.IGRA_TB) %>%
  rename(., iNMB.IGRA_RNA_TB_upper = iNMB.IGRA_RNA_TB) %>%
  rename(., iNMB.RNA_TB_upper = iNMB.RNA_TB) %>% 
  rename(., iNMB.IGRA_RNA_TB150_upper = iNMB.IGRA_RNA_TB150) %>%
  rename(., iNMB.RNA_TB150_upper = iNMB.RNA_TB150) %>%
  rename(., iNMB.IGRA_RNA_TB60_upper = iNMB.IGRA_RNA_TB60) %>%
  rename(., iNMB.RNA_TB60_upper = iNMB.RNA_TB60) %>%
  rename(., iNMB.IGRA_RNA_TB30_upper = iNMB.IGRA_RNA_TB30) %>%
  rename(., iNMB.RNA_TB30_upper = iNMB.RNA_TB30) %>%
  rename(., iNMB.IGRA_RNA_TB15_upper = iNMB.IGRA_RNA_TB15) %>%
  rename(., iNMB.RNA_TB15_upper = iNMB.RNA_TB15)

nmb_by_coo_lower = apply(nmb_by_coo_array2, c(1:2), quantile, p = 0.025) # dataframe 97 x 36
colnames(nmb_by_coo_lower) = colnames(nmb_by_coo_list[[1]])
nmb_by_coo_lower = nmb_by_coo_lower %>% 
  as.data.frame() %>%
  select(c(iNMB.IGRA_TB, iNMB.IGRA_RNA_TB, iNMB.RNA_TB, 
           iNMB.IGRA_RNA_TB150, iNMB.IGRA_RNA_TB60, iNMB.IGRA_RNA_TB30, iNMB.IGRA_RNA_TB15,
           iNMB.RNA_TB150, iNMB.RNA_TB60, iNMB.RNA_TB30, iNMB.RNA_TB15)) %>% 
  rename(., iNMB.IGRA_TB_lower = iNMB.IGRA_TB) %>%
  rename(., iNMB.IGRA_RNA_TB_lower = iNMB.IGRA_RNA_TB) %>%
  rename(., iNMB.RNA_TB_lower = iNMB.RNA_TB) %>%
  rename(., iNMB.IGRA_RNA_TB150_lower = iNMB.IGRA_RNA_TB150) %>%
  rename(., iNMB.RNA_TB150_lower = iNMB.RNA_TB150) %>%
  rename(., iNMB.IGRA_RNA_TB60_lower = iNMB.IGRA_RNA_TB60) %>%
  rename(., iNMB.RNA_TB60_lower = iNMB.RNA_TB60) %>%
  rename(., iNMB.IGRA_RNA_TB30_lower = iNMB.IGRA_RNA_TB30) %>%
  rename(., iNMB.RNA_TB30_lower = iNMB.RNA_TB30) %>%
  rename(., iNMB.IGRA_RNA_TB15_lower = iNMB.IGRA_RNA_TB15) %>%
  rename(., iNMB.RNA_TB15_lower = iNMB.RNA_TB15)


nmb_by_coo = cbind(nmb_by_coo_mean2, nmb_by_coo_lower, nmb_by_coo_upper)


riskcat_I = coo_by_TBinc[coo_by_TBinc$risk_category == "I_300plus",]$country
riskcat_II = coo_by_TBinc[coo_by_TBinc$risk_category == "II_100_300",]$country
riskcat_III = coo_by_TBinc[coo_by_TBinc$risk_category == "III_10_100",]$country
riskcat_IV = coo_by_TBinc[coo_by_TBinc$risk_category == "VI_0_10",]$country

coo_by_TBinc_temp = coo_by_TBinc %>% rename(., country = place_of_birth)


temp = merge( coo_by_TBinc_temp, nmb_by_coo ,by = "country")

nmb_by_riskcat = temp %>% group_by(risk_category) %>%
  summarise(NMB.No_screening = weighted.mean(NMB.No_screening, total_pop),
            NMB.IGRA_TB      = weighted.mean(NMB.IGRA_TB, total_pop),
            NMB.IGRA_RNA_TB  = weighted.mean(NMB.IGRA_RNA_TB, total_pop),
            iNMB.IGRA_TB     = weighted.mean(iNMB.IGRA_TB, total_pop),
            iNMB.IGRA_RNA_TB = weighted.mean(iNMB.IGRA_RNA_TB, total_pop),
            iNMB.RNA_TB      = weighted.mean(iNMB.RNA_TB, total_pop),
            iNMB.IGRA_TB_lower     = weighted.mean(iNMB.IGRA_TB_lower, total_pop),
            iNMB.IGRA_RNA_TB_lower = weighted.mean(iNMB.IGRA_RNA_TB_lower, total_pop),
            iNMB.RNA_TB_lower      = weighted.mean(iNMB.RNA_TB_lower, total_pop),
            iNMB.IGRA_TB_upper     = weighted.mean(iNMB.IGRA_TB_upper, total_pop),
            iNMB.IGRA_RNA_TB_upper = weighted.mean(iNMB.IGRA_RNA_TB_upper, total_pop),
            iNMB.RNA_TB_upper      = weighted.mean(iNMB.RNA_TB_upper, total_pop))

nmb_by_riskcat15 = temp %>% group_by(risk_category) %>%
  summarise(NMB.No_screening = weighted.mean(NMB.No_screening, total_pop),
            NMB.IGRA_TB      = weighted.mean(NMB.IGRA_TB, total_pop),
            NMB.IGRA_RNA_TB  = weighted.mean(NMB.IGRA_RNA_TB, total_pop),
            iNMB.IGRA_TB     = weighted.mean(iNMB.IGRA_TB, total_pop),
            iNMB.IGRA_RNA_TB = weighted.mean(iNMB.IGRA_RNA_TB15, total_pop),
            iNMB.RNA_TB      = weighted.mean(iNMB.RNA_TB15, total_pop),
            iNMB.IGRA_TB_lower     = weighted.mean(iNMB.IGRA_TB_lower, total_pop),
            iNMB.IGRA_RNA_TB_lower = weighted.mean(iNMB.IGRA_RNA_TB15_lower, total_pop),
            iNMB.RNA_TB_lower      = weighted.mean(iNMB.RNA_TB15_lower, total_pop),
            iNMB.IGRA_TB_upper     = weighted.mean(iNMB.IGRA_TB_upper, total_pop),
            iNMB.IGRA_RNA_TB_upper = weighted.mean(iNMB.IGRA_RNA_TB15_upper, total_pop),
            iNMB.RNA_TB_upper      = weighted.mean(iNMB.RNA_TB15_upper, total_pop))

nmb_by_riskcat30 = temp %>% group_by(risk_category) %>%
  summarise(NMB.No_screening = weighted.mean(NMB.No_screening, total_pop),
            NMB.IGRA_TB      = weighted.mean(NMB.IGRA_TB, total_pop),
            NMB.IGRA_RNA_TB  = weighted.mean(NMB.IGRA_RNA_TB, total_pop),
            iNMB.IGRA_TB     = weighted.mean(iNMB.IGRA_TB, total_pop),
            iNMB.IGRA_RNA_TB = weighted.mean(iNMB.IGRA_RNA_TB30, total_pop),
            iNMB.RNA_TB      = weighted.mean(iNMB.RNA_TB30, total_pop),
            iNMB.IGRA_TB_lower     = weighted.mean(iNMB.IGRA_TB_lower, total_pop),
            iNMB.IGRA_RNA_TB_lower = weighted.mean(iNMB.IGRA_RNA_TB30_lower, total_pop),
            iNMB.RNA_TB_lower      = weighted.mean(iNMB.RNA_TB30_lower, total_pop),
            iNMB.IGRA_TB_upper     = weighted.mean(iNMB.IGRA_TB_upper, total_pop),
            iNMB.IGRA_RNA_TB_upper = weighted.mean(iNMB.IGRA_RNA_TB30_upper, total_pop),
            iNMB.RNA_TB_upper      = weighted.mean(iNMB.RNA_TB30_upper, total_pop))

nmb_by_riskcat60 = temp %>% group_by(risk_category) %>%
  summarise(NMB.No_screening = weighted.mean(NMB.No_screening, total_pop),
            NMB.IGRA_TB      = weighted.mean(NMB.IGRA_TB, total_pop),
            NMB.IGRA_RNA_TB  = weighted.mean(NMB.IGRA_RNA_TB, total_pop),
            iNMB.IGRA_TB     = weighted.mean(iNMB.IGRA_TB, total_pop),
            iNMB.IGRA_RNA_TB = weighted.mean(iNMB.IGRA_RNA_TB60, total_pop),
            iNMB.RNA_TB      = weighted.mean(iNMB.RNA_TB60, total_pop),
            iNMB.IGRA_TB_lower     = weighted.mean(iNMB.IGRA_TB_lower, total_pop),
            iNMB.IGRA_RNA_TB_lower = weighted.mean(iNMB.IGRA_RNA_TB60_lower, total_pop),
            iNMB.RNA_TB_lower      = weighted.mean(iNMB.RNA_TB60_lower, total_pop),
            iNMB.IGRA_TB_upper     = weighted.mean(iNMB.IGRA_TB_upper, total_pop),
            iNMB.IGRA_RNA_TB_upper = weighted.mean(iNMB.IGRA_RNA_TB60_upper, total_pop),
            iNMB.RNA_TB_upper      = weighted.mean(iNMB.RNA_TB60_upper, total_pop))

nmb_by_riskcat150 = temp %>% group_by(risk_category) %>%
  summarise(NMB.No_screening = weighted.mean(NMB.No_screening, total_pop),
            NMB.IGRA_TB      = weighted.mean(NMB.IGRA_TB, total_pop),
            NMB.IGRA_RNA_TB  = weighted.mean(NMB.IGRA_RNA_TB, total_pop),
            iNMB.IGRA_TB     = weighted.mean(iNMB.IGRA_TB, total_pop),
            iNMB.IGRA_RNA_TB = weighted.mean(iNMB.IGRA_RNA_TB150, total_pop),
            iNMB.RNA_TB      = weighted.mean(iNMB.RNA_TB150, total_pop),
            iNMB.IGRA_TB_lower     = weighted.mean(iNMB.IGRA_TB_lower, total_pop),
            iNMB.IGRA_RNA_TB_lower = weighted.mean(iNMB.IGRA_RNA_TB150_lower, total_pop),
            iNMB.RNA_TB_lower      = weighted.mean(iNMB.RNA_TB150_lower, total_pop),
            iNMB.IGRA_TB_upper     = weighted.mean(iNMB.IGRA_TB_upper, total_pop),
            iNMB.IGRA_RNA_TB_upper = weighted.mean(iNMB.IGRA_RNA_TB150_upper, total_pop),
            iNMB.RNA_TB_upper      = weighted.mean(iNMB.RNA_TB150_upper, total_pop))

plot.df.inmb = get_plot_df.inmb(nmb_by_riskcat)
plot.df.inmb15 = get_plot_df.inmb(nmb_by_riskcat15)
plot.df.inmb30 = get_plot_df.inmb(nmb_by_riskcat30)
plot.df.inmb60 = get_plot_df.inmb(nmb_by_riskcat60)
plot.df.inmb150 = get_plot_df.inmb(nmb_by_riskcat150)
plot.df.inmb_full = rbind(plot.df.inmb, plot.df.inmb150, plot.df.inmb60, plot.df.inmb30, plot.df.inmb15)
plot.df.inmb_full$cost = rep(c("$300", "$150", "$60", "$30", "$15"), each = nrow(plot.df.inmb))
plot.df.inmb_full$cost_f = factor(plot.df.inmb_full$cost, levels = c("$300", "$150", "$60", "$30", "$15"))


return(plot.df.inmb_full)
}


inmb150000 = get_iNMB_byWTP(150000)
inmb100000 = get_iNMB_byWTP(100000)
inmb50000  = get_iNMB_byWTP(50000)
inmb30000  = get_iNMB_byWTP(30000)

plot.df.inmb_full2 = rbind(inmb150000, inmb100000, inmb50000, inmb30000)
plot.df.inmb_full2$wtp = rep(c("150,000", "100,000", "50,000", "30,000"), each = nrow(inmb150000))
plot.df.inmb_full2$wtp_f = factor(plot.df.inmb_full2$wtp, levels = c("150,000", "100,000","50,000", "30,000"))

max_y = 7250
min_y = -450

ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df.inmb_full2) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df.inmb_full2, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  xlab("Risk category") +
  ylab("Per-person iNMB (USD)") +
  # xlab(" ") +
  # ylab(" ") +
  # labs(title = "$300") +
  ggtitle(expression(bold(atop("Incremental net monetary benefit relative to Strategy I", "by cost of the signature (USD) and WTP (USD / QALY gained)"))))+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  scale_y_continuous(limits = c(min_y - 100, max_y + 100)) + 
  facet_grid(wtp_f~cost_f) + 
  theme_line() + 
  theme(
    legend.position = "bottom",
    #     legend.justification = c("right", "top"),
        # legend.box.just = "bottom",
    #     legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(angle = 0), 
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        plot.margin = margin(15, 15, 15, 15),
        title =element_text(size=12, face='bold')
  )

# ==============================================================================
# APPENDIX 7 Figure 2 four-way sensitivity analysis
# ==============================================================================

# set results directory to PSA_calib_rnaGAM

CURRENT_RESDIR = "/n/netscratch/menzies_lab/Everyone/ylh202/IncipienTB/PSA_RNAgam"

nmb_by_coo_list = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_all_30.RDS"))


# nmb_by_coo_list = future_lapply(1:1000, function(x) as.data.frame(nmb_by_coo_array[,,x]))  # list of 1000 dataframes 97 x 36


for (i in c(1:1000)){
  colnames(nmb_by_coo_list[[i]]) = c("eff.dQALY.No_testing", "eff.dQALY.IGRA_TB", "eff.dQALY.IGRA_RNA_TB", "eff.dQALY.RNA_TB",
                                     "eff.QALY.No_testing", "eff.QALY.IGRA_TB", "eff.QALY.IGRA_RNA_TB", "eff.QALY.RNA_TB",
                                     "eff.dLY.No_testing", "eff.dLY.IGRA_TB", "eff.dLY.IGRA_RNA_TB", "eff.dLY.RNA_TB",
                                     "eff.LY.No_testing", "eff.LY.IGRA_TB", "eff.LY.IGRA_RNA_TB", "eff.LY.RNA_TB",
                                     "cost.No_testing", "cost.IGRA_TB", "cost.IGRA_RNA_TB", "cost.RNA_TB",
                                     "cost150.No_testing", "cost150.IGRA_TB", "cost150.IGRA_RNA_TB", "cost150.RNA_TB",
                                     "cost60.No_testing", "cost60.IGRA_TB", "cost60.IGRA_RNA_TB", "cost60.RNA_TB",
                                     "cost30.No_testing", "cost30.IGRA_TB", "cost30.IGRA_RNA_TB", "cost30.RNA_TB",
                                     "cost15.No_testing", "cost15.IGRA_TB", "cost15.IGRA_RNA_TB", "cost15.RNA_TB")
}


for (i in c(1:1000)){
  nmb_by_coo_list[[i]] = nmb_by_coo_list[[i]] %>%
    mutate(NMB.No_screening = eff.dQALY.No_testing*WTP - cost.No_testing,
           NMB.IGRA_TB      = eff.dQALY.IGRA_TB*WTP - cost.IGRA_TB,
           NMB.IGRA_RNA_TB    = eff.dQALY.IGRA_RNA_TB*WTP - cost.IGRA_RNA_TB,
           NMB.IGRA_RNA_TB150 = eff.dQALY.IGRA_RNA_TB*WTP - cost150.IGRA_RNA_TB,
           NMB.IGRA_RNA_TB60  = eff.dQALY.IGRA_RNA_TB*WTP -  cost60.IGRA_RNA_TB,
           NMB.IGRA_RNA_TB30  = eff.dQALY.IGRA_RNA_TB*WTP -  cost30.IGRA_RNA_TB,
           NMB.IGRA_RNA_TB15  = eff.dQALY.IGRA_RNA_TB*WTP -  cost15.IGRA_RNA_TB,
           NMB.RNA_TB     = eff.dQALY.RNA_TB*WTP - cost.RNA_TB,
           NMB.RNA_TB150  = eff.dQALY.RNA_TB*WTP - cost150.RNA_TB,
           NMB.RNA_TB60   = eff.dQALY.RNA_TB*WTP - cost60.RNA_TB,
           NMB.RNA_TB30   = eff.dQALY.RNA_TB*WTP - cost30.RNA_TB,
           NMB.RNA_TB15   = eff.dQALY.RNA_TB*WTP - cost15.RNA_TB,
           iNMB.IGRA_TB        =  (eff.dQALY.IGRA_TB*WTP - cost.IGRA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB     =  (eff.dQALY.IGRA_RNA_TB*WTP - cost.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB150  =  (eff.dQALY.IGRA_RNA_TB*WTP - cost150.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB60   =  (eff.dQALY.IGRA_RNA_TB*WTP - cost60.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB30   =  (eff.dQALY.IGRA_RNA_TB*WTP - cost30.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.IGRA_RNA_TB15   =  (eff.dQALY.IGRA_RNA_TB*WTP - cost15.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB    =  (eff.dQALY.RNA_TB*WTP - cost.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB150 =  (eff.dQALY.RNA_TB*WTP - cost150.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB60  =  (eff.dQALY.RNA_TB*WTP - cost60.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB30  =  (eff.dQALY.RNA_TB*WTP - cost30.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing),
           iNMB.RNA_TB15  =  (eff.dQALY.RNA_TB*WTP - cost15.RNA_TB) - (eff.dQALY.No_testing*WTP - cost.No_testing))
}

nmb_by_coo_df = do.call("rbind", nmb_by_coo_list)


# re-import PSAtable
nmb_by_coo_df$paramset = rep(1:1000, each = 97)
nmb_by_coo_df$country  = rep(countrynames, 1000)
nmb_by_coo_df$risk_category = rep(coo_by_TBinc$risk_category[order(coo_by_TBinc$place_of_birth)], 1000)
nmb_by_coo_df$total_pop = rep(cohort2019_by_coo2$total_pop[order(cohort2019_by_coo2$place_of_birth)], 1000)

nmb_by_riskcat = nmb_by_coo_df %>% group_by(paramset, risk_category) %>%
  summarise(NMB.No_screening = weighted.mean(NMB.No_screening, total_pop),
            NMB.IGRA_TB = weighted.mean(NMB.IGRA_TB, total_pop),
            NMB.IGRA_RNA_TB    = weighted.mean(NMB.IGRA_RNA_TB, total_pop),
            NMB.IGRA_RNA_TB150 = weighted.mean(NMB.IGRA_RNA_TB150, total_pop),
            NMB.IGRA_RNA_TB60  = weighted.mean(NMB.IGRA_RNA_TB60, total_pop),
            NMB.IGRA_RNA_TB30  = weighted.mean(NMB.IGRA_RNA_TB30, total_pop),
            NMB.IGRA_RNA_TB15  = weighted.mean(NMB.IGRA_RNA_TB15, total_pop),
            NMB.RNA_TB    = weighted.mean(NMB.RNA_TB, total_pop),
            NMB.RNA_TB150 = weighted.mean(NMB.RNA_TB150, total_pop),
            NMB.RNA_TB60  = weighted.mean(NMB.RNA_TB60, total_pop),
            NMB.RNA_TB30  = weighted.mean(NMB.RNA_TB30, total_pop),
            NMB.RNA_TB15  = weighted.mean(NMB.RNA_TB15, total_pop),
            iNMB.IGRA_TB = weighted.mean(iNMB.IGRA_TB, total_pop),
            iNMB.IGRA_RNA_TB    = weighted.mean(iNMB.IGRA_RNA_TB, total_pop),
            iNMB.IGRA_RNA_TB150 = weighted.mean(iNMB.IGRA_RNA_TB150, total_pop),
            iNMB.IGRA_RNA_TB60  = weighted.mean(iNMB.IGRA_RNA_TB60, total_pop),
            iNMB.IGRA_RNA_TB30  = weighted.mean(iNMB.IGRA_RNA_TB30, total_pop),
            iNMB.IGRA_RNA_TB15  = weighted.mean(iNMB.IGRA_RNA_TB15, total_pop),
            iNMB.RNA_TB    = weighted.mean(iNMB.RNA_TB, total_pop),
            iNMB.RNA_TB150 = weighted.mean(iNMB.RNA_TB150, total_pop),
            iNMB.RNA_TB60  = weighted.mean(iNMB.RNA_TB60, total_pop),
            iNMB.RNA_TB30  = weighted.mean(iNMB.RNA_TB30, total_pop),
            iNMB.RNA_TB15  = weighted.mean(iNMB.RNA_TB15, total_pop))

PSAtable = read.csv("Data/PSAtable1000_gam.csv")
nmb_by_riskcat$YEARS = rep(PSAtable$RNA_TIME_upper, each = 4)
nmb_by_riskcat$SENS1 = rep(PSAtable$RNA_SENS1, each = 4)
nmb_by_riskcat$SPEC  = rep(PSAtable$RNA_SPEC, each = 4)

nmb_by_riskcat = nmb_by_riskcat %>%
  mutate(optimal.300 = NA,
         optimal.150 = NA, 
         optimal.60 = NA, 
         optimal.30 = NA,
         optimal.15 = NA, 
         s2.300 = NA, 
         s2.150 = NA, 
         s2.60 = NA, 
         s2.30 = NA, 
         s2.15 = NA,
         s3.300 = NA, 
         s3.150 = NA, 
         s3.60 = NA, 
         s3.30 = NA, 
         s3.15 = NA)
for(i in 1:nrow(nmb_by_riskcat)){
  nmb_by_riskcat$optimal.300[i] = which.max(c(nmb_by_riskcat$NMB.No_screening[i], nmb_by_riskcat$NMB.IGRA_TB[i], nmb_by_riskcat$NMB.IGRA_RNA_TB[i])) - 1
  nmb_by_riskcat$optimal.150[i] = which.max(c(nmb_by_riskcat$NMB.No_screening[i], nmb_by_riskcat$NMB.IGRA_TB[i], nmb_by_riskcat$NMB.IGRA_RNA_TB150[i])) - 1
  nmb_by_riskcat$optimal.60[i] = which.max(c(nmb_by_riskcat$NMB.No_screening[i], nmb_by_riskcat$NMB.IGRA_TB[i], nmb_by_riskcat$NMB.IGRA_RNA_TB60[i])) - 1
  nmb_by_riskcat$optimal.30[i] = which.max(c(nmb_by_riskcat$NMB.No_screening[i], nmb_by_riskcat$NMB.IGRA_TB[i], nmb_by_riskcat$NMB.IGRA_RNA_TB30[i])) - 1
  nmb_by_riskcat$optimal.15[i] = which.max(c(nmb_by_riskcat$NMB.No_screening[i], nmb_by_riskcat$NMB.IGRA_TB[i], nmb_by_riskcat$NMB.IGRA_RNA_TB15[i])) - 1
  
  nmb_by_riskcat$s2.300[i] = ifelse( nmb_by_riskcat$optimal.300[i] == 2, 1, 0)
  nmb_by_riskcat$s2.150[i] = ifelse( nmb_by_riskcat$optimal.150[i] == 2, 1, 0)
  nmb_by_riskcat$s2.60[i] = ifelse( nmb_by_riskcat$optimal.60[i] == 2, 1, 0)
  nmb_by_riskcat$s2.30[i] = ifelse( nmb_by_riskcat$optimal.30[i] == 2, 1, 0)
  nmb_by_riskcat$s2.15[i] = ifelse( nmb_by_riskcat$optimal.15[i] == 2, 1, 0)
  
  nmb_by_riskcat$s3.300[i] = ifelse( nmb_by_riskcat$optimal.300[i] == 3, 1, 0)
  nmb_by_riskcat$s3.150[i] = ifelse( nmb_by_riskcat$optimal.150[i] == 3, 1, 0)
  nmb_by_riskcat$s3.60[i] = ifelse( nmb_by_riskcat$optimal.60[i] == 3, 1, 0)
  nmb_by_riskcat$s3.30[i] = ifelse( nmb_by_riskcat$optimal.30[i] == 3, 1, 0)
  nmb_by_riskcat$s3.15[i] = ifelse( nmb_by_riskcat$optimal.15[i] == 3, 1, 0)
} 


nmb_by_riskcat1 = nmb_by_riskcat %>% filter(risk_category == "I_300plus")
nmb_by_riskcat2 = nmb_by_riskcat %>% filter(risk_category == "II_100_300")
nmb_by_riskcat3 = nmb_by_riskcat %>% filter(risk_category == "III_10_100")
nmb_by_riskcat4 = nmb_by_riskcat %>% filter(risk_category == "VI_0_10")

# countour plot

get_diff = function(dat){
  
  dat$diff.300 = dat$iNMB.IGRA_RNA_TB - dat$iNMB.IGRA_TB
  dat$diff.150 = dat$iNMB.IGRA_RNA_TB150 - dat$iNMB.IGRA_TB
  dat$diff.60  = dat$iNMB.IGRA_RNA_TB60 - dat$iNMB.IGRA_TB
  dat$diff.30  = dat$iNMB.IGRA_RNA_TB30 - dat$iNMB.IGRA_TB
  dat$diff.15  = dat$iNMB.IGRA_RNA_TB15 - dat$iNMB.IGRA_TB
  
  return(dat)
}


get_diff2 = function(dat){
  
  dat$diff.300 = dat$iNMB.RNA_TB - dat$iNMB.IGRA_TB
  dat$diff.150 = dat$iNMB.RNA_TB150 - dat$iNMB.IGRA_TB
  dat$diff.60  = dat$iNMB.RNA_TB60 - dat$iNMB.IGRA_TB
  dat$diff.30  = dat$iNMB.RNA_TB30 - dat$iNMB.IGRA_TB
  dat$diff.15  = dat$iNMB.RNA_TB15 - dat$iNMB.IGRA_TB
  
  return(dat)
}


get_contour_plot_df = function(dat){
  
  dat = get_diff(dat)
  dat$diff.300 = dat$iNMB.IGRA_RNA_TB - dat$iNMB.IGRA_TB
  dat$diff.150 = dat$iNMB.IGRA_RNA_TB150 - dat$iNMB.IGRA_TB
  dat$diff.60  = dat$iNMB.IGRA_RNA_TB60 - dat$iNMB.IGRA_TB
  dat$diff.30  = dat$iNMB.IGRA_RNA_TB30 - dat$iNMB.IGRA_TB
  dat$diff.15  = dat$iNMB.IGRA_RNA_TB15 - dat$iNMB.IGRA_TB
  
  out.rsm300 = rsm::rsm(diff.300 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm150 = rsm::rsm(diff.150 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm60 = rsm::rsm(diff.60 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm30 = rsm::rsm(diff.30 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm15 = rsm::rsm(diff.15 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  
  return(out.list = list(out.rsm300 = out.rsm300, out.rsm150 = out.rsm150, out.rsm60 = out.rsm60, out.rsm30 = out.rsm30, out.rsm15 = out.rsm15))
}

get_contour_plot_df2 = function(dat){
  
  dat = get_diff2(dat)
  dat$diff.300 = dat$iNMB.RNA_TB - dat$iNMB.IGRA_TB
  dat$diff.150 = dat$iNMB.RNA_TB150 - dat$iNMB.IGRA_TB
  dat$diff.60  = dat$iNMB.RNA_TB60 - dat$iNMB.IGRA_TB
  dat$diff.30  = dat$iNMB.RNA_TB30 - dat$iNMB.IGRA_TB
  dat$diff.15  = dat$iNMB.RNA_TB15 - dat$iNMB.IGRA_TB
  
  out.rsm300 = rsm::rsm(diff.300 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm150 = rsm::rsm(diff.150 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm60 = rsm::rsm(diff.60 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm30 = rsm::rsm(diff.30 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm15 = rsm::rsm(diff.15 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  
  return(out.list = list(out.rsm300 = out.rsm300, out.rsm150 = out.rsm150, out.rsm60 = out.rsm60, out.rsm30 = out.rsm30, out.rsm15 = out.rsm15))
}


plot_riskcat = function(dat){
  op = par(mfrow=c(2,4),mar=c(4,4,2, 2), oma = c(5, 5, 3, 1), cex.axis= 1.2, cex.lab= 1.3)
  
  
  # b = lapply(3:10, function(x) contour(ct_df, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
  #                                      xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
  #                                      xlim = c(0.75, 1), ylim = c(0.75, 1)))
  b = for (x in 3:10) { contour(dat, ~ TWI(SENS1,SPEC), at = data.frame(YEARS = x), image = T,
                                # xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
                                xlabs = c("  ", " "),
                                xlim = c(0.75, 1), ylim = c(0.75, 1), 
                                # labcex = 0.75,
                                lwd = 1)} 
  mtext("P(T+ | incipient TB)",side=1,line=2,outer=TRUE, at = 0.5)
  mtext("P(T- | no TB in lifetime)",side=2,line=0,outer=TRUE, las=0, at = 0.5)
}

ct_df_risk1 = get_contour_plot_df(nmb_by_riskcat1)
ct_df_risk2 = get_contour_plot_df(nmb_by_riskcat2)
ct_df_risk3 = get_contour_plot_df(nmb_by_riskcat3)
ct_df_risk4 = get_contour_plot_df(nmb_by_riskcat4)

ct_df_risk1_s4 = get_contour_plot_df2(nmb_by_riskcat1)
ct_df_risk2_s4 = get_contour_plot_df2(nmb_by_riskcat2)
ct_df_risk3_s4 = get_contour_plot_df2(nmb_by_riskcat3)
ct_df_risk4_s4 = get_contour_plot_df2(nmb_by_riskcat4)

# WTP 150k/QALY, Risk Category 1
png(file=paste0(getwd(),"/ct_risk1_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[1]]) 
title(main = "Risk Category: >300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_150.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[2]]) 
title(main = "Risk Category: >300, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_60.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[3]]) 
title(main = "Risk Category: >300, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_30.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[4]]) 
title(main = "Risk Category: >300, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[5]]) 
title(main = "Risk Category: >300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

# WTP 150k/QALY, Risk Category 2
png(file=paste0(getwd(),"/ct_risk2_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk2[[1]]) 
title(main = "Risk Category: 100-300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk2_150.png"), width=1070, height=600)
plot_riskcat(ct_df_risk2[[2]]) 
title(main = "Risk Category: 100-300, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk2_60.png"), width=1070, height=600)
plot_riskcat(ct_df_risk2[[3]]) 
title(main = "Risk Category: 100-300, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk2_30.png"), width=1070, height=600)
plot_riskcat(ct_df_risk2[[4]]) 
title(main = "Risk Category: 100-300, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk2_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk2[[5]]) 
title(main = "Risk Category: 100-300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()




# WTP 150k/QALY, Risk Category 3
png(file=paste0(getwd(),"/ct_risk3_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk3[[1]]) 
title(main = "Risk Category: 10-100, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk3_150.png"), width=1070, height=600)
plot_riskcat(ct_df_risk3[[2]]) 
title(main = "Risk Category: 10-100, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk3_60.png"), width=1070, height=600)
plot_riskcat(ct_df_risk3[[3]]) 
title(main = "Risk Category: 10-100, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk3_30.png"), width=1070, height=600)
plot_riskcat(ct_df_risk3[[4]]) 
title(main = "Risk Category: 10-100, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk3_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk3[[5]]) 
title(main = "Risk Category: 10-100, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()


# WTP 150k/QALY, Risk Category 4
png(file=paste0(getwd(),"/ct_risk4_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4[[1]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_150.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4[[2]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_60.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4[[3]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_30.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4[[4]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4[[5]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

# --------------

# WTP 150k/QALY, Risk Category 1
png(file=paste0(getwd(),"/ct_risk1_300_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4[[1]]) 
title(main = "Risk Category: >300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_150_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4[[2]]) 
title(main = "Risk Category: >300, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_60_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4[[3]]) 
title(main = "Risk Category: >300, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_30_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4[[4]]) 
title(main = "Risk Category: >300, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk1_15_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4[[5]]) 
title(main = "Risk Category: >300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

# WTP 150k/QALY, Risk Category 4
png(file=paste0(getwd(),"/ct_risk4_300_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4[[1]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_150_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4[[2]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_60_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4[[3]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_30_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4[[4]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/ct_risk4_15_s4.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4[[4]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()


## =============================================================================
## obsolete code 
## =============================================================================

c300 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df300) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df300, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  #xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  #ylab("Per-person increase in cost (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$300")+
  scale_y_continuous(limits = c(0, 450)) +
  #ggtitle("Increase in cost relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  theme_line() + 
  theme( 
    #legend.position = "None",        
    axis.text.x = element_text(angle =0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title =element_text(size=12, face='bold')
  ) #560ox500


c150 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df150) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df150, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  #xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  #ylab("Per-person increase in cost (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$150")+
  scale_y_continuous(limits = c(0, 450)) +
  #ggtitle("Increase in cost relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  theme_line() + 
  theme( 
    #legend.position = "None",        
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title =element_text(size=12, face='bold')
  ) #560ox500

c60 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df60) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df60, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  #xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  #ylab("Per-person increase in cost (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$60")+
  scale_y_continuous(limits = c(0, 450)) +
  #ggtitle("Increase in cost relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  theme_line() + 
  theme( 
    #legend.position = "None",        
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title =element_text(size=12, face='bold')
  ) #560ox500

c30 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df30) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df30, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  #xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  #ylab("Per-person increase in cost (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$30")+
  scale_y_continuous(limits = c(0, 450)) +
  #ggtitle("Increase in cost relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  theme_line() + 
  theme( 
    #legend.position = "None",        
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title =element_text(size=12, face='bold')
  ) #560ox500

c15 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df15) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df15, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  #xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  #ylab("Per-person increase in cost (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$15")+
  scale_y_continuous(limits = c(0, 450)) +
  # ggtitle("Increase in cost relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  theme_line() + 
  theme( 
    #legend.position = "None",        
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title =element_text(size=12, face='bold')
  ) #560ox500

gg_xaxis = cowplot::get_plot_component(ggplot() +
                                         labs(x = expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)")),
                                         ) + 
                                         theme(text = element_text(size = 12)), "xlab-b",)
gg_yaxis = cowplot::get_plot_component(ggplot() + 
                                         labs(y = "Per-person increase in cost (USD)") + 
                                         theme(text = element_text(size = 12)), "ylab-l")
gg_title = cowplot::get_plot_component(c15, "^title")

fig = ggarrange(c300, c150, c60, c30, c15, labels = NULL, ncol=5, nrow=1, common.legend = TRUE, legend="right")
annotate_figure(fig, left = gg_yaxis, bottom= gg_xaxis, 
                top = text_grob("Increase in cost relative to Strategy I, by costs of the signature test (USD)",  size = 18, face = "bold"))





a150 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df.inmb150) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df.inmb150, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  # xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  # ylab("Per-person iNMB (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$150") +
  # ggtitle("Incremental net monetary benefit relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  scale_y_continuous(limits = c(min_y - 100, max_y + 100)) + 
  theme_line() + 
  theme(
    # legend.position = c(.95, .95),
    #     legend.justification = c("right", "top"),
    #     legend.box.just = "right",
    #     legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(angle = 0), 
        plot.margin = margin(10, 5.5, 0, 5.5),
        title=element_text(size = 12, face = 'bold'))

a60 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df.inmb60) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df.inmb60, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  # xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  # ylab("Per-person iNMB (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$60") +
  # ggtitle("Incremental net monetary benefit relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  scale_y_continuous(limits = c(min_y - 100, max_y + 100)) + 
  theme_line() + 
  theme(
    # legend.position = c(.95, .95),
    #     legend.justification = c("right", "top"),
    #     legend.box.just = "right",
    #     legend.margin = margin(6, 6, 6, 6),
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title=element_text(size = 12, face = 'bold'))

a30 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df.inmb30) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df.inmb30, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  # xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  # ylab("Per-person iNMB (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$30") +
  # ggtitle("Incremental net monetary benefit relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  scale_y_continuous(limits = c(min_y - 100, max_y + 100)) + 
  theme_line() + 
  theme(
    # legend.position = c(.95, .95),
    #     legend.justification = c("right", "top"),
    #     legend.box.just = "right",
    #     legend.margin = margin(6, 6, 6, 6),
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title=element_text(size = 12, face = 'bold'))

a15 = ggplot(aes(y = mean, x= as.factor(risk_category), fill = variable), data = plot.df.inmb15) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(risk_category)), data = plot.df.inmb15, 
                width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  # xlab(expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)"))) + 
  # ylab("Per-person iNMB (USD)") +
  xlab(" ") +
  ylab(" ") +
  labs(title = "$15") +
  # ggtitle("Incremental net monetary benefit relative to Strategy I")+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
  scale_y_continuous(limits = c(min_y - 100, max_y + 100)) + 
  theme_line() + 
  theme(
    # legend.position = c(.95, .95),
    #     legend.justification = c("right", "top"),
    #     legend.box.just = "right",
    #     legend.margin = margin(6, 6, 6, 6),
    axis.text.x = element_text(angle = 0), 
    plot.margin = margin(10, 5.5, 0, 5.5),
    title=element_text(size = 12, face = 'bold'))

gg_xaxis = cowplot::get_plot_component(ggplot() +
                                         labs(x = expression(atop("Risk category","(TB incidence per 100k of country-of-origin in 2019)")),
                                         ) + 
                                         theme(text = element_text(size = 12)), "xlab-b",)
gg_yaxis = cowplot::get_plot_component(ggplot() + 
                                         labs(y = "Per-person iNMB (USD)") + 
                                         theme(text = element_text(size = 12)), "ylab-l")
gg_title = cowplot::get_plot_component(a15, "^title")


fig = ggarrange(a300, a150, a60, a30, a15, labels = NULL, ncol=5, nrow=1, common.legend = T)
fig = annotate_figure(fig)


# max_y = max(plot.df.inmb$upper, plot.df.inmb150$upper, plot.df.inmb60$upper, plot.df.inmb30$upper, plot.df.inmb15$upper) # 7201.507
# min_y = min(plot.df.inmb$lower, plot.df.inmb150$lower, plot.df.inmb60$lower, plot.df.inmb30$lower, plot.df.inmb15$lower) # -449.0603


w30000 = get_iNMB_byWTP(30000)
w50000 = get_iNMB_byWTP(50000)
w100000 = get_iNMB_byWTP(100000)
w150000 = get_iNMB_byWTP(150000)

FIG = ggarrange(w150000, w100000, w50000, w30000, labels = NULL, ncol = 1, nrow=4, common.legend = T, legend = 'right')
FIG = annotate_figure(FIG, left = gg_yaxis, bottom= gg_xaxis, 
                      text_grob("Incremental net monetary benefit relative to Strategy I, by costs of the signature test (USD)",  size = 18, face = "bold"))
