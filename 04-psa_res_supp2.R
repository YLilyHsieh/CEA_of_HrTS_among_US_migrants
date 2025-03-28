
CURRENT_RESDIR = "~/IncipienTB/PSA_calib_rnaGAM"


# expected expenditures for one riskcat at one RNA cost
# read in data
get_mean_expend = function(rdsFileName){
  dat = readRDS(paste0(CURRENT_RESDIR,"/",rdsFileName))
  dat_array = array(unlist(dat), dim = c(3, 5, 1000))
  out = apply(dat_array, c(1,2), function(x) mean(x))
  out.df = as.data.frame(out)
  colnames(out.df) =c("strategy", "TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist")
  
  return(out.df)
}

get_mean_cost_by_risk = function(meanExpendDF, meanProd, HCpersp, risknum){
  
  if(HCpersp){
    out = rowSums(meanExpendDF[,2:3]) - t(meanProd[risknum,2:4])
  } else {
    out = rowSums(meanExpendDF[,2:5]) - t(meanProd[risknum,2:4])
  }
  
  return(out)
} 

get_mean_cost = function(rdsFileRisk1, rdsFileRisk2, rdsFileRisk3, rdsFileRisk4, HCpersp){
  
  meanExpendRisk1 = get_mean_expend(rdsFileRisk1)
  meanExpendRisk2 = get_mean_expend(rdsFileRisk2)
  meanExpendRisk3 = get_mean_expend(rdsFileRisk3)
  meanExpendRisk4 = get_mean_expend(rdsFileRisk4)
  
  meanCostRisk1 = get_mean_cost_by_risk(meanExpendRisk1, prod_mean, HCpersp, 1)
  meanCostRisk2 = get_mean_cost_by_risk(meanExpendRisk1, prod_mean, HCpersp, 2)
  meanCostRisk3 = get_mean_cost_by_risk(meanExpendRisk1, prod_mean, HCpersp, 3)
  meanCostRisk4 = get_mean_cost_by_risk(meanExpendRisk1, prod_mean, HCpersp, 4)
  
  meanCost = t(cbind(meanCostRisk1, meanCostRisk2, meanCostRisk3, meanCostRisk4))
  rownames(meanCost) = c("risk1", "risk2", "risk3", "risk4")
  return(meanCost)
}

# effectiveness 
# recalculate to get more decimals 
get_mean_eff = function(){
  dat = read.csv(paste0(CURRENT_RESDIR,"/nmb_by_coo_mean2.csv"))
  
  risk1 = dat %>% filter(country %in% riskcat_I) %>% 
    summarise(IGRA_TB = weighted.mean(eff.dQALY.IGRA_TB - eff.dQALY.No_testing, pop_size),
              IGRA_RNA_TB = weighted.mean(eff.dQALY.IGRA_RNA_TB - eff.dQALY.No_testing, pop_size),
              RNA_TB = weighted.mean(eff.dQALY.RNA_TB - eff.dQALY.No_testing, pop_size))
 
  risk2 = dat %>% filter(country %in% riskcat_II) %>% 
    summarise(IGRA_TB = weighted.mean(eff.dQALY.IGRA_TB - eff.dQALY.No_testing, pop_size),
              IGRA_RNA_TB = weighted.mean(eff.dQALY.IGRA_RNA_TB - eff.dQALY.No_testing, pop_size),
              RNA_TB = weighted.mean(eff.dQALY.RNA_TB - eff.dQALY.No_testing, pop_size))
  
  risk3 = dat %>% filter(country %in% riskcat_III) %>% 
    summarise(IGRA_TB = weighted.mean(eff.dQALY.IGRA_TB - eff.dQALY.No_testing, pop_size),
              IGRA_RNA_TB = weighted.mean(eff.dQALY.IGRA_RNA_TB - eff.dQALY.No_testing, pop_size),
              RNA_TB = weighted.mean(eff.dQALY.RNA_TB - eff.dQALY.No_testing, pop_size))
  
  risk4 = dat %>% filter(country %in% riskcat_IV) %>% 
    summarise(IGRA_TB = weighted.mean(eff.dQALY.IGRA_TB - eff.dQALY.No_testing, pop_size),
              IGRA_RNA_TB = weighted.mean(eff.dQALY.IGRA_RNA_TB - eff.dQALY.No_testing, pop_size),
              RNA_TB = weighted.mean(eff.dQALY.RNA_TB - eff.dQALY.No_testing, pop_size))
  out = rbind(risk1, risk2, risk3, risk4)
  return(out)
}

get_iNMB = function(rdsFileRisk1, rdsFileRisk2, rdsFileRisk3, rdsFileRisk4, HCpersp, WTP){
  
  meanCost = get_mean_cost(rdsFileRisk1, rdsFileRisk2, rdsFileRisk3, rdsFileRisk4, HCpersp)
  meanEff = get_mean_eff()
  
  iNMB_df = meanEff*WTP - meanCost
  iNMB_df$riskcat = c(1:4)
  iNMB_df$wtp = WTP
  
  out = gather(iNMB_df, "strategy", "iNMB", IGRA_TB:RNA_TB)
  return(out)
}

get_iNMB_byWTP = function(HCpersp, WTP){
  
  iNMB_300 = get_iNMB("nmb_by_coo_itm_risk1.RDS", "nmb_by_coo_itm_risk2.RDS", "nmb_by_coo_itm_risk3.RDS", "nmb_by_coo_itm_risk4.RDS", HCpersp = HCpersp, WTP = WTP)
  iNMB_150 = get_iNMB("nmb_by_coo_itm_risk1_costItems150.RDS", "nmb_by_coo_itm_risk2_costItems150.RDS", "nmb_by_coo_itm_risk3_costItems150.RDS", "nmb_by_coo_itm_risk4_costItems150.RDS", HCpersp = HCpersp, WTP = WTP)
  iNMB_60  = get_iNMB("nmb_by_coo_itm_risk1_costItems60.RDS", "nmb_by_coo_itm_risk2_costItems60.RDS", "nmb_by_coo_itm_risk3_costItems60.RDS", "nmb_by_coo_itm_risk4_costItems60.RDS", HCpersp = HCpersp, WTP = WTP)
  iNMB_30  = get_iNMB("nmb_by_coo_itm_risk1_costItems30.RDS", "nmb_by_coo_itm_risk2_costItems30.RDS", "nmb_by_coo_itm_risk3_costItems30.RDS", "nmb_by_coo_itm_risk4_costItems30.RDS", HCpersp = HCpersp, WTP = WTP)
  iNMB_15  = get_iNMB("nmb_by_coo_itm_risk1_costItems15.RDS", "nmb_by_coo_itm_risk2_costItems15.RDS", "nmb_by_coo_itm_risk3_costItems15.RDS", "nmb_by_coo_itm_risk4_costItems15.RDS", HCpersp = HCpersp, WTP = WTP)
  
  out = rbind(iNMB_300, iNMB_150, iNMB_60, iNMB_30, iNMB_15)
  out$RNAcost = c(rep(c(300, 150, 60, 30, 15), each = 12))

  return(out)
}

read_nmb_by_riskcat = function(costLevel){
  assign(paste0("nmb_by_coo_itm_risk1_",costLevel), readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_",costLevel,"_dft5.RDS")))
  assign(paste0("nmb_by_coo_itm_risk2_",costLevel), readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_",costLevel,"_dft5.RDS")))
  assign(paste0("nmb_by_coo_itm_risk3_",costLevel), readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_",costLevel,"_dft5.RDS")))
  assign(paste0("nmb_by_coo_itm_risk4_",costLevel), readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_",costLevel,"_dft5.RDS")))
}


for(i in c("costItems150","costItems60", "costItems30", "costItems15")){
  read_nmb_by_riskcat(costLevel = i)
}


prod_mean = read.csv(paste0(CURRENT_RESDIR,"/prod_df.csv"))


# format data for plotting
iNMB_150k_HC = get_iNMB_byWTP(HCpersp = T, WTP = 150000)
iNMB_100k_HC = get_iNMB_byWTP(HCpersp = T, WTP = 100000)
iNMB_50k_HC = get_iNMB_byWTP(HCpersp = T, WTP = 50000)
iNMB_30k_HC = get_iNMB_byWTP(HCpersp = T, WTP = 30000)

iNMB_150k_SOC = get_iNMB_byWTP(HCpersp = F, WTP = 150000)
iNMB_100k_SOC = get_iNMB_byWTP(HCpersp = F, WTP = 100000)
iNMB_50k_SOC = get_iNMB_byWTP(HCpersp = F, WTP = 50000)
iNMB_30k_SOC = get_iNMB_byWTP(HCpersp = F, WTP = 30000)

plot.iNMB_df_HC = rbind(iNMB_150k_HC, iNMB_100k_HC, iNMB_50k_HC, iNMB_30k_HC)
plot.iNMB_df_SOC = rbind(iNMB_150k_SOC, iNMB_100k_SOC, iNMB_50k_SOC, iNMB_30k_SOC)

# create factor variables for plotting
plot.iNMB_df_HC$wtp_f = as.factor(-plot.iNMB_df_HC$wtp)
levels(plot.iNMB_df_HC$wtp_f) = c("150,000", "100,000", "50,000", "30,000")
plot.iNMB_df_HC$cost_f =  as.factor(-plot.iNMB_df_HC$RNAcost)
levels(plot.iNMB_df_HC$cost_f) = c("$300", "$150","$60", "$30", "$15")
plot.iNMB_df_HC$strategy_f = factor(plot.iNMB_df_HC$strategy, levels = c("IGRA_TB", "IGRA_RNA_TB","RNA_TB"))

# plot! 

max_y = 5000
min_y = -400

ggplot(aes(y = iNMB, x= as.factor(riskcat), fill = strategy_f), data = plot.iNMB_df_HC) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  # geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(riskcat)), data = plot.iNMB_df_HC, 
  #               width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  xlab("Risk Category") +
  ylab("Per-person iNMB (USD)") +
  # xlab(" ") +
  # ylab(" ") +
  # labs(title = "$300") +
  #ggtitle(expression(bold(atop("Inc. NMB relative to 'I. No-screening' in a Healthcare sector perspective", "by costs of HrTS and WTP thresholds"))))+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-only", "III: IGRA-HrTB", "IV: HrTS-only"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
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


# create factor variables for plotting
plot.iNMB_df_SOC$wtp_f = as.factor(-plot.iNMB_df_SOC$wtp)
levels(plot.iNMB_df_SOC$wtp_f) = c("150,000", "100,000", "50,000", "30,000")
plot.iNMB_df_SOC$cost_f =  as.factor(-plot.iNMB_df_SOC$RNAcost)
levels(plot.iNMB_df_SOC$cost_f) = c("$300", "$150","$60", "$30", "$15")
plot.iNMB_df_SOC$strategy_f = factor(plot.iNMB_df_SOC$strategy, levels = c("IGRA_TB", "IGRA_RNA_TB","RNA_TB"))

# plot! 

max_y = 5000
min_y = -400

ggplot(aes(y = iNMB, x= as.factor(riskcat), fill = strategy_f), data = plot.iNMB_df_SOC) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) + 
  # geom_errorbar(aes(ymin = lower, ymax = upper, x= as.factor(riskcat)), data = plot.iNMB_df_HC, 
  #               width=0.1, colour="black",alpha = 0.8, size=0.5, position=position_dodge(0.9)) + 
  xlab("Risk Category") +
  ylab("Per-person iNMB (USD)") +
  # xlab(" ") +
  # ylab(" ") +
  # labs(title = "$300") +
  #ggtitle(expression(bold(atop("Inc. NMB relative to 'I. No-screening' in a Healthcare sector perspective", "by costs of HrTS and WTP thresholds"))))+
  scale_fill_manual(
    name = "Strategy", labels = c( "II: IGRA-only", "III: IGRA-HrTB", "IV: HrTS-only"),values=c("#00BA38", "#619CFF", "grey70")) +
  scale_x_discrete(labels = c("300+", "100-300", "100-10", "0-10"))+
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
# APPENDIX 7 Figure 2 four-way sensitivity analysis
# ==============================================================================

corgi = get_nmb_by_coo(paste0(CURRENT_RESDIR,"/paramset1/cea_SOM.RDS")) # just a random file to get the dimensions
resNames = names(corgi)
corgi = as.data.frame(t(corgi))

# calculate iNMB 

nmb_by_coo_array = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_array.RDS"))
nmb_by_coo_list = future_lapply(1:1000, function(x) as.data.frame(nmb_by_coo_array[,,x]))  # list of 1000 dataframes 97 x 36


for (i in c(1:1000)){
  colnames(nmb_by_coo_list[[i]]) = resNames
}


for (i in c(1:1000)){
  nmb_by_coo_list[[i]] = nmb_by_coo_list[[i]]  %>%
    mutate(NMB.HC.No_screening = eff.dQALY.No_testing*WTP - costHC.No_testing,
           NMB.HC.IGRA_TB      = eff.dQALY.IGRA_TB*WTP - costHC.IGRA_TB,
           NMB.HC.IGRA_RNA_TB  = eff.dQALY.IGRA_RNA_TB*WTP - costHC.IGRA_RNA_TB,
           NMB.HC.RNA_TB       = eff.dQALY.RNA_TB*WTP - costHC.RNA_TB,
           iNMB.HC.IGRA_TB        =  (eff.dQALY.IGRA_TB*WTP - costHC.IGRA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.IGRA_RNA_TB    =  (eff.dQALY.IGRA_RNA_TB*WTP - costHC.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.IGRA_RNA_TB150 =  (eff.dQALY.IGRA_RNA_TB*WTP - costHC150.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.IGRA_RNA_TB60  =  (eff.dQALY.IGRA_RNA_TB*WTP - costHC60.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.IGRA_RNA_TB30  =  (eff.dQALY.IGRA_RNA_TB*WTP - costHC30.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.IGRA_RNA_TB15  =  (eff.dQALY.IGRA_RNA_TB*WTP - costHC15.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.RNA_TB    =  (eff.dQALY.RNA_TB*WTP - costHC.RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.RNA_TB150 =  (eff.dQALY.RNA_TB*WTP - costHC150.RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.RNA_TB60  =  (eff.dQALY.RNA_TB*WTP - costHC60.RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.RNA_TB30  =  (eff.dQALY.RNA_TB*WTP - costHC30.RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing),
           iNMB.HC.RNA_TB15  =  (eff.dQALY.RNA_TB*WTP - costHC15.RNA_TB) - (eff.dQALY.No_testing*WTP - costHC.No_testing)) %>%
    
    mutate(prodGain.IGRA_TB     = prod.IGRA_TB     - prod.No_testing,
           prodGain.IGRA_RNA_TB = prod.IGRA_RNA_TB - prod.No_testing,
           prodGain.RNA_TB      = prod.RNA_TB      - prod.No_testing) %>%
    mutate(costSOC.IGRA_TB        = costSOC.IGRA_TB        - prodGain.IGRA_TB,
           costSOC.IGRA_RNA_TB     = costSOC.IGRA_RNA_TB   - prodGain.IGRA_RNA_TB,
           costSOC150.IGRA_RNA_TB  = costSOC150.IGRA_RNA_TB - prodGain.IGRA_RNA_TB,
           costSOC60.IGRA_RNA_TB   = costSOC60.IGRA_RNA_TB  - prodGain.IGRA_RNA_TB,
           costSOC30.IGRA_RNA_TB   = costSOC30.IGRA_RNA_TB  - prodGain.IGRA_RNA_TB,
           costSOC15.IGRA_RNA_TB   = costSOC15.IGRA_RNA_TB  - prodGain.IGRA_RNA_TB,
           costSOC.RNA_TB    = costSOC.RNA_TB    - prodGain.RNA_TB,
           costSOC150.RNA_TB = costSOC150.RNA_TB - prodGain.RNA_TB,
           costSOC60.RNA_TB  = costSOC60.RNA_TB  - prodGain.RNA_TB,
           costSOC30.RNA_TB  = costSOC30.RNA_TB  - prodGain.RNA_TB,
           costSOC15.RNA_TB  = costSOC15.RNA_TB  - prodGain.RNA_TB) %>%
    mutate(NMB.SOC.No_screening = eff.dQALY.No_testing*WTP - costSOC.No_testing,
           NMB.SOC.IGRA_TB      = eff.dQALY.IGRA_TB*WTP - costSOC.IGRA_TB,
           NMB.SOC.IGRA_RNA_TB  = eff.dQALY.IGRA_RNA_TB*WTP - costSOC.IGRA_RNA_TB,
           NMB.SOC.RNA_TB       = eff.dQALY.RNA_TB*WTP - costSOC.RNA_TB,
           iNMB.SOC.IGRA_TB        =  (eff.dQALY.IGRA_TB*WTP - costSOC.IGRA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.IGRA_RNA_TB    =  (eff.dQALY.IGRA_RNA_TB*WTP - costSOC.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.IGRA_RNA_TB150 =  (eff.dQALY.IGRA_RNA_TB*WTP - costSOC150.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.IGRA_RNA_TB60  =  (eff.dQALY.IGRA_RNA_TB*WTP - costSOC60.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.IGRA_RNA_TB30  =  (eff.dQALY.IGRA_RNA_TB*WTP - costSOC30.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.IGRA_RNA_TB15  =  (eff.dQALY.IGRA_RNA_TB*WTP - costSOC15.IGRA_RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.RNA_TB    =  (eff.dQALY.RNA_TB*WTP - costSOC.RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.RNA_TB150 =  (eff.dQALY.RNA_TB*WTP - costSOC150.RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.RNA_TB60  =  (eff.dQALY.RNA_TB*WTP - costSOC60.RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.RNA_TB30  =  (eff.dQALY.RNA_TB*WTP - costSOC30.RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing),
           iNMB.SOC.RNA_TB15  =  (eff.dQALY.RNA_TB*WTP - costSOC15.RNA_TB) - (eff.dQALY.No_testing*WTP - costSOC.No_testing))
}

nmb_by_coo_df = do.call("rbind", nmb_by_coo_list)

# re-import PSAtable
nmb_by_coo_df$paramset = rep(1:1000, each = 97)
nmb_by_coo_df$country  = rep(countrynames, 1000)
nmb_by_coo_df$risk_category = rep(coo_by_TBinc$risk_category[order(coo_by_TBinc$place_of_birth)], 1000)
nmb_by_coo_df$total_pop = rep(cohort2019_by_coo2$total_pop[order(cohort2019_by_coo2$place_of_birth)], 1000)

NMB.SOC_by_riskcat = nmb_by_coo_df %>% group_by(paramset, risk_category) %>%
  summarise(
    iNMB.SOC.IGRA_TB = weighted.mean(iNMB.SOC.IGRA_TB, total_pop),
    iNMB.SOC.IGRA_RNA_TB    = weighted.mean(iNMB.SOC.IGRA_RNA_TB, total_pop),
    iNMB.SOC.IGRA_RNA_TB150 = weighted.mean(iNMB.SOC.IGRA_RNA_TB150, total_pop),
    iNMB.SOC.IGRA_RNA_TB60  = weighted.mean(iNMB.SOC.IGRA_RNA_TB60, total_pop),
    iNMB.SOC.IGRA_RNA_TB30  = weighted.mean(iNMB.SOC.IGRA_RNA_TB30, total_pop),
    iNMB.SOC.IGRA_RNA_TB15  = weighted.mean(iNMB.SOC.IGRA_RNA_TB15, total_pop),
    iNMB.SOC.RNA_TB    = weighted.mean(iNMB.SOC.RNA_TB, total_pop),
    iNMB.SOC.RNA_TB150 = weighted.mean(iNMB.SOC.RNA_TB150, total_pop),
    iNMB.SOC.RNA_TB60  = weighted.mean(iNMB.SOC.RNA_TB60, total_pop),
    iNMB.SOC.RNA_TB30  = weighted.mean(iNMB.SOC.RNA_TB30, total_pop),
    iNMB.SOC.RNA_TB15  = weighted.mean(iNMB.SOC.RNA_TB15, total_pop))

PSAtable = read.csv("Data/PSAtable1000_gam.csv")
NMB.SOC_by_riskcat$YEARS = rep(PSAtable$RNA_TIME_upper, each = 4)
NMB.SOC_by_riskcat$SENS1 = rep(PSAtable$RNA_SENS1, each = 4)
NMB.SOC_by_riskcat$SPEC  = rep(PSAtable$RNA_SPEC, each = 4)

NMB.HC_by_riskcat = nmb_by_coo_df %>% group_by(paramset, risk_category) %>%
  summarise(
    iNMB.HC.IGRA_TB = weighted.mean(iNMB.HC.IGRA_TB, total_pop),
    iNMB.HC.IGRA_RNA_TB    = weighted.mean(iNMB.HC.IGRA_RNA_TB, total_pop),
    iNMB.HC.IGRA_RNA_TB150 = weighted.mean(iNMB.HC.IGRA_RNA_TB150, total_pop),
    iNMB.HC.IGRA_RNA_TB60  = weighted.mean(iNMB.HC.IGRA_RNA_TB60, total_pop),
    iNMB.HC.IGRA_RNA_TB30  = weighted.mean(iNMB.HC.IGRA_RNA_TB30, total_pop),
    iNMB.HC.IGRA_RNA_TB15  = weighted.mean(iNMB.HC.IGRA_RNA_TB15, total_pop),
    iNMB.HC.RNA_TB    = weighted.mean(iNMB.HC.RNA_TB, total_pop),
    iNMB.HC.RNA_TB150 = weighted.mean(iNMB.HC.RNA_TB150, total_pop),
    iNMB.HC.RNA_TB60  = weighted.mean(iNMB.HC.RNA_TB60, total_pop),
    iNMB.HC.RNA_TB30  = weighted.mean(iNMB.HC.RNA_TB30, total_pop),
    iNMB.HC.RNA_TB15  = weighted.mean(iNMB.HC.RNA_TB15, total_pop))

# PSAtable = read.csv("Data/PSAtable1000_gam.csv") # draft 5
PSAtable = read.csv("Data/PSAtable1000_gam_draft3draft4.csv") 
NMB.HC_by_riskcat$YEARS = rep(PSAtable$RNA_TIME_upper, each = 4)
NMB.HC_by_riskcat$SENS1 = rep(PSAtable$RNA_SENS1, each = 4)
NMB.HC_by_riskcat$SPEC  = rep(PSAtable$RNA_SPEC, each = 4)

# countour plot

get_diff = function(dat, persp){
  
  if (persp == "SOC") {
    dat$diff.300 = dat$iNMB.SOC.IGRA_RNA_TB - dat$iNMB.SOC.IGRA_TB
    dat$diff.150 = dat$iNMB.SOC.IGRA_RNA_TB150 - dat$iNMB.SOC.IGRA_TB
    dat$diff.60  = dat$iNMB.SOC.IGRA_RNA_TB60 - dat$iNMB.SOC.IGRA_TB
    dat$diff.30  = dat$iNMB.SOC.IGRA_RNA_TB30 - dat$iNMB.SOC.IGRA_TB
    dat$diff.15  = dat$iNMB.SOC.IGRA_RNA_TB15 - dat$iNMB.SOC.IGRA_TB
  } else if (persp == "HC"){
    dat$diff.300 = dat$iNMB.HC.IGRA_RNA_TB - dat$iNMB.HC.IGRA_TB
    dat$diff.150 = dat$iNMB.HC.IGRA_RNA_TB150 - dat$iNMB.HC.IGRA_TB
    dat$diff.60  = dat$iNMB.HC.IGRA_RNA_TB60 - dat$iNMB.HC.IGRA_TB
    dat$diff.30  = dat$iNMB.HC.IGRA_RNA_TB30 - dat$iNMB.HC.IGRA_TB
    dat$diff.15  = dat$iNMB.HC.IGRA_RNA_TB15 - dat$iNMB.HC.IGRA_TB
  }

  return(dat)
}


get_diff2 = function(dat, persp){
  
  if (persp == "SOC"){
    dat$diff.300 = dat$iNMB.SOC.RNA_TB - dat$iNMB.SOC.IGRA_TB
    dat$diff.150 = dat$iNMB.SOC.RNA_TB150 - dat$iNMB.SOC.IGRA_TB
    dat$diff.60  = dat$iNMB.SOC.RNA_TB60 - dat$iNMB.SOC.IGRA_TB
    dat$diff.30  = dat$iNMB.SOC.RNA_TB30 - dat$iNMB.SOC.IGRA_TB
    dat$diff.15  = dat$iNMB.SOC.RNA_TB15 - dat$iNMB.SOC.IGRA_TB
  } else if (persp == "HC"){
    dat$diff.300 = dat$iNMB.HC.RNA_TB - dat$iNMB.HC.IGRA_TB
    dat$diff.150 = dat$iNMB.HC.RNA_TB150 - dat$iNMB.HC.IGRA_TB
    dat$diff.60  = dat$iNMB.HC.RNA_TB60 - dat$iNMB.HC.IGRA_TB
    dat$diff.30  = dat$iNMB.HC.RNA_TB30 - dat$iNMB.HC.IGRA_TB
    dat$diff.15  = dat$iNMB.HC.RNA_TB15 - dat$iNMB.HC.IGRA_TB
  }
  
  return(dat)
}

get_contour_plot_df = function(dat, persp){
  
  dat = get_diff(dat, persp)
  
  if (persp == "SOC"){
    dat$diff.300 = dat$iNMB.SOC.IGRA_RNA_TB - dat$iNMB.SOC.IGRA_TB
    dat$diff.150 = dat$iNMB.SOC.IGRA_RNA_TB150 - dat$iNMB.SOC.IGRA_TB
    dat$diff.60  = dat$iNMB.SOC.IGRA_RNA_TB60 - dat$iNMB.SOC.IGRA_TB
    dat$diff.30  = dat$iNMB.SOC.IGRA_RNA_TB30 - dat$iNMB.SOC.IGRA_TB
    dat$diff.15  = dat$iNMB.SOC.IGRA_RNA_TB15 - dat$iNMB.SOC.IGRA_TB
  } else if (persp == "HC"){
    dat$diff.300 = dat$iNMB.HC.IGRA_RNA_TB - dat$iNMB.HC.IGRA_TB
    dat$diff.150 = dat$iNMB.HC.IGRA_RNA_TB150 - dat$iNMB.HC.IGRA_TB
    dat$diff.60  = dat$iNMB.HC.IGRA_RNA_TB60 - dat$iNMB.HC.IGRA_TB
    dat$diff.30  = dat$iNMB.HC.IGRA_RNA_TB30 - dat$iNMB.HC.IGRA_TB
    dat$diff.15  = dat$iNMB.HC.IGRA_RNA_TB15 - dat$iNMB.HC.IGRA_TB
  }
  
  out.rsm300 = rsm::rsm(diff.300 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm150 = rsm::rsm(diff.150 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm60 = rsm::rsm(diff.60 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm30 = rsm::rsm(diff.30 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm15 = rsm::rsm(diff.15 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  
  return(out.list = list(out.rsm300 = out.rsm300, out.rsm150 = out.rsm150, out.rsm60 = out.rsm60, out.rsm30 = out.rsm30, out.rsm15 = out.rsm15))
}

get_contour_plot_df2 = function(dat, persp){
  
  dat = get_diff2(dat, persp)
  
  if (persp == "SOC"){
    dat$diff.300 = dat$iNMB.SOC.RNA_TB - dat$iNMB.SOC.IGRA_TB
    dat$diff.150 = dat$iNMB.SOC.RNA_TB150 - dat$iNMB.SOC.IGRA_TB
    dat$diff.60  = dat$iNMB.SOC.RNA_TB60 - dat$iNMB.SOC.IGRA_TB
    dat$diff.30  = dat$iNMB.SOC.RNA_TB30 - dat$iNMB.SOC.IGRA_TB
    dat$diff.15  = dat$iNMB.SOC.RNA_TB15 - dat$iNMB.SOC.IGRA_TB
  } else if (persp == "HC"){
    dat$diff.300 = dat$iNMB.HC.RNA_TB - dat$iNMB.HC.IGRA_TB
    dat$diff.150 = dat$iNMB.HC.RNA_TB150 - dat$iNMB.HC.IGRA_TB
    dat$diff.60  = dat$iNMB.HC.RNA_TB60 - dat$iNMB.HC.IGRA_TB
    dat$diff.30  = dat$iNMB.HC.RNA_TB30 - dat$iNMB.HC.IGRA_TB
    dat$diff.15  = dat$iNMB.HC.RNA_TB15 - dat$iNMB.HC.IGRA_TB
  }
 
  out.rsm300 = rsm::rsm(diff.300 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm150 = rsm::rsm(diff.150 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm60 = rsm::rsm(diff.60 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm30 = rsm::rsm(diff.30 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm15 = rsm::rsm(diff.15 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  
  return(out.list = list(out.rsm300 = out.rsm300, out.rsm150 = out.rsm150, out.rsm60 = out.rsm60, out.rsm30 = out.rsm30, out.rsm15 = out.rsm15))
}

get_contour_plot_df3 = function(dat, strategy, persp){
  
  if (strategy == 4){
    if (persp == "SOC"){
      dat$diff.300 = dat$iNMB.SOC.RNA_TB 
      dat$diff.150 = dat$iNMB.SOC.RNA_TB150 
      dat$diff.60  = dat$iNMB.SOC.RNA_TB60 
      dat$diff.30  = dat$iNMB.SOC.RNA_TB30 
      dat$diff.15  = dat$iNMB.SOC.RNA_TB15 
    } else if (persp == "HC"){
      dat$diff.300 = dat$iNMB.HC.RNA_TB 
      dat$diff.150 = dat$iNMB.HC.RNA_TB150 
      dat$diff.60  = dat$iNMB.HC.RNA_TB60 
      dat$diff.30  = dat$iNMB.HC.RNA_TB30 
      dat$diff.15  = dat$iNMB.HC.RNA_TB15
      }
    } else if (strategy == 3){
    if (persp == "SOC"){
      dat$diff.300 = dat$iNMB.SOC.IGRA_RNA_TB 
      dat$diff.150 = dat$iNMB.SOC.IGRA_RNA_TB150 
      dat$diff.60  = dat$iNMB.SOC.IGRA_RNA_TB60 
      dat$diff.30  = dat$iNMB.SOC.IGRA_RNA_TB30 
      dat$diff.15  = dat$iNMB.SOC.IGRA_RNA_TB15 
    } else if (persp == "HC"){
      dat$diff.300 = dat$iNMB.HC.IGRA_RNA_TB 
      dat$diff.150 = dat$iNMB.HC.IGRA_RNA_TB150 
      dat$diff.60  = dat$iNMB.HC.IGRA_RNA_TB60 
      dat$diff.30  = dat$iNMB.HC.IGRA_RNA_TB30 
      dat$diff.15  = dat$iNMB.HC.IGRA_RNA_TB15 
    }
  } else if (strategy == 2){
    if (persp == "SOC"){
      dat$diff.300 = dat$iNMB.SOC.IGRA_TB
      dat$diff.150 = dat$iNMB.SOC.IGRA_TB
      dat$diff.60  = dat$iNMB.SOC.IGRA_TB 
      dat$diff.30  = dat$iNMB.SOC.IGRA_TB 
      dat$diff.15  = dat$iNMB.SOC.IGRA_TB 
    } else if (persp == "HC"){
      dat$diff.300 = dat$iNMB.HC.IGRA_TB 
      dat$diff.150 = dat$iNMB.HC.IGRA_TB
      dat$diff.60  = dat$iNMB.HC.IGRA_TB 
      dat$diff.30  = dat$iNMB.HC.IGRA_TB 
      dat$diff.15  = dat$iNMB.HC.IGRA_TB 
    }
  }
  
  out.rsm300 = rsm::rsm(diff.300 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm150 = rsm::rsm(diff.150 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm60 = rsm::rsm(diff.60 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm30 = rsm::rsm(diff.30 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  out.rsm15 = rsm::rsm(diff.15 ~  FO(SENS1,SPEC, YEARS) + TWI(SENS1,SPEC, YEARS), data = dat)
  
  return(out.list = list(out.rsm300 = out.rsm300, out.rsm150 = out.rsm150, out.rsm60 = out.rsm60, out.rsm30 = out.rsm30, out.rsm15 = out.rsm15))
}

plot_riskcat = function(dat){
  op = par(mfrow=c(2,4),mar=c(4,4,2.5, 2.5), oma = c(6, 6, 3, 1), cex.axis= 1.5, cex.lab= 1.8)  #1.2, 1.3
  
  # b = lapply(3:10, function(x) contour(ct_df, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
  #                                      xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
  #                                      xlim = c(0.75, 1), ylim = c(0.75, 1)))
  b = for (x in 3:10) { contour(dat, ~ TWI(SENS1,SPEC), at = data.frame(YEARS = x), image = T,
                                # xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
                                xlabs = c("  ", " "),
                                xlim = c(0.75, 1), ylim = c(0.75, 1), 
                                las = 1,
                                labcex = 1.3,
                                vfont = c("sans serif", "bold"),
                                lwd = 1)} 
  mtext("P(T+ | incipient TB)",side=1,line=3,outer=TRUE, at = 0.5, cex = 2)
  mtext("P(T- | no TB in lifetime)",side=2,line=3,outer=TRUE, las=0, at = 0.5, cex = 2)
}

plot_riskcat_pretty = function(dat){
  op = par(mfrow=c(2,4),mar=c(4,4,2.5, 2.5), oma = c(6, 6, 3, 1), cex.axis= 1.5, cex.lab= 1.8)  #1.2, 1.3
  
  # b = lapply(3:10, function(x) contour(ct_df, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
  #                                      xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
  #                                      xlim = c(0.75, 1), ylim = c(0.75, 1)))
  b = for (x in 3:10) { contour(dat, ~ TWI(SENS1,SPEC), at = data.frame(YEARS = x),  image = T,
                                # xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
                                xlabs = c("  ", " "),
                                xlim = c(0.75, 1), ylim = c(0.75, 1), 
                                las = 1,
                                labcex = 1.3,
                                vfont = c("sans serif", "bold"),
                                lwd = 1,
                                atpos = 0 #,col = hcl.colors(100, "terrain")
                                )
  # image(col = hcl.colors(20, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    #title(main = paste0("Timeframe = ",x," years"))
    mtext(paste0("Timeframe = ",x," yrs"), side = 1, line = 3.5, cex = 1.1)} 
  mtext("P(T+ | incipient TB)",side=1,line=3,outer=TRUE, at = 0.5, cex = 2)
  mtext("P(T- | no TB in lifetime)",side=2,line=3,outer=TRUE, las=0, at = 0.5, cex = 2)
}



## Societal Perspective --------------------------------------------------------
NMB.SOC_by_riskcat1 = NMB.SOC_by_riskcat %>% filter(risk_category == "I_300plus")
NMB.SOC_by_riskcat2 = NMB.SOC_by_riskcat %>% filter(risk_category == "II_100_300")
NMB.SOC_by_riskcat3 = NMB.SOC_by_riskcat %>% filter(risk_category == "III_10_100")
NMB.SOC_by_riskcat4 = NMB.SOC_by_riskcat %>% filter(risk_category == "VI_0_10")

ct_df_risk1_s3_soc = get_contour_plot_df(NMB.SOC_by_riskcat1, "SOC")
# ct_df_risk2_s3 = get_contour_plot_df(NMB.SOC_by_riskcat2, "SOC")
# ct_df_risk3_s3 = get_contour_plot_df(NMB.SOC_by_riskcat3, "SOC")
ct_df_risk4_s3_soc = get_contour_plot_df(NMB.SOC_by_riskcat4, "SOC")

ct_df_risk1_s4_soc = get_contour_plot_df2(NMB.SOC_by_riskcat1, "SOC")
# ct_df_risk2_s4 = get_contour_plot_df2(NMB.SOC_by_riskcat2, "SOC")
# ct_df_risk3_s4 = get_contour_plot_df2(NMB.SOC_by_riskcat3, "SOC")
ct_df_risk4_s4_soc = get_contour_plot_df2(NMB.SOC_by_riskcat4, "SOC")

ct_df_risk4_s4_soc_s2 = get_contour_plot_df3(NMB.SOC_by_riskcat4, 2, "SOC")
ct_df_risk4_s4_soc_s3 = get_contour_plot_df3(NMB.SOC_by_riskcat4, 3, "SOC")
ct_df_risk4_s4_soc_s4 = get_contour_plot_df3(NMB.SOC_by_riskcat4, 4, "SOC")


NMB.HC_by_riskcat1 = NMB.HC_by_riskcat %>% filter(risk_category == "I_300plus")
NMB.HC_by_riskcat2 = NMB.HC_by_riskcat %>% filter(risk_category == "II_100_300")
NMB.HC_by_riskcat3 = NMB.HC_by_riskcat %>% filter(risk_category == "III_10_100")
NMB.HC_by_riskcat4 = NMB.HC_by_riskcat %>% filter(risk_category == "VI_0_10")

ct_df_risk1_s3_hc = get_contour_plot_df(NMB.HC_by_riskcat1, "HC")
ct_df_risk4_s3_hc = get_contour_plot_df(NMB.HC_by_riskcat4, "HC")

ct_df_risk1_s4_hc = get_contour_plot_df2(NMB.HC_by_riskcat1, "HC")
ct_df_risk4_s4_hc = get_contour_plot_df2(NMB.HC_by_riskcat4, "HC")

ct_df_risk4_s4_hc = get_contour_plot_df3(NMB.HC_by_riskcat4, "HC")

## SOCIETAL PERSPECTIVE
# WTP 150k/QALY, Risk Category 1, Strategy 3
png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk1_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s3_soc[[1]]) 
title(main = "Risk Category: >300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk1_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s3_soc[[5]]) 
title(main = "Risk Category: >300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

# png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk1_150.png"), width=1070, height=600)
# plot_riskcat(ct_df_risk1_s3_soc[[2]]) 
# title(main = "Risk Category: >300, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 3)
# dev.off()
# 
# png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk1_60.png"), width=1070, height=600)
# plot_riskcat(ct_df_risk1_s3_soc[[3]]) 
# title(main = "Risk Category: >300, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 3)
# dev.off()
# 
# png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk1_30.png"), width=1070, height=600)
# plot_riskcat(ct_df_risk1_s3_soc[[4]]) 
# title(main = "Risk Category: >300, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 3)
# dev.off()

# WTP 150k/QALY, Risk Category 4, Strategy 3
png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk4_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s3_soc[[1]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s3_risk4_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s3_soc[[5]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

# WTP 150k/QALY, Risk Category 1, Strategy 4
png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s4_risk1_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4_soc[[1]]) 
title(main = "Risk Category: >300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s4_risk1_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4_soc[[5]]) 
title(main = "Risk Category: >300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

# WTP 150k/QALY, Risk Category 4, Strategy 4
png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s4_risk4_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4_soc[[1]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_soc_s4_risk4_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4_soc[[5]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

## HEATHCARE SECTOR PERSPECTIVE
# WTP 150k/QALY, Risk Category 1, Strategy 3
png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk1_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s3_hc[[1]]) 
title(main = "Risk Category: >300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk1_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s3_hc[[5]]) 
title(main = "Risk Category: >300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

# png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk1_150.png"), width=1070, height=600)
# plot_riskcat(ct_df_risk1_s3_hc[[2]]) 
# title(main = "Risk Category: >300, Cost of RNA: 150 USD", outer = TRUE,xpd=NA, cex.main = 3)
# dev.off()
# 
# png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk1_60.png"), width=1070, height=600)
# plot_riskcat(ct_df_risk1_s3_hc[[3]]) 
# title(main = "Risk Category: >300, Cost of RNA: 60 USD", outer = TRUE,xpd=NA, cex.main = 3)
# dev.off()
# 
# png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk1_30.png"), width=1070, height=600)
# plot_riskcat(ct_df_risk1_s3_hc[[4]]) 
# title(main = "Risk Category: >300, Cost of RNA: 30 USD", outer = TRUE,xpd=NA, cex.main = 3)
# dev.off()

# WTP 150k/QALY, Risk Category 4, Strategy 3
png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk4_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s3_hc[[1]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s3_risk4_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s3_hc[[5]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

# WTP 150k/QALY, Risk Category 1, Strategy 4
png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s4_risk1_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4_hc[[1]]) 
title(main = "Risk Category: >300, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s4_risk1_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1_s4_hc[[5]]) 
title(main = "Risk Category: >300, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

# WTP 150k/QALY, Risk Category 4, Strategy 4
png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s4_risk4_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4_hc[[1]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 300 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()

png(file=paste0(getwd(),"/figures.draft3/append7/ct_hc_s4_risk4_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk4_s4_hc[[5]]) 
title(main = "Risk Category: 0-10, Cost of RNA: 15 USD", outer = TRUE,xpd=NA, cex.main = 3)
dev.off()