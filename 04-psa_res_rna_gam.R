# The data compilation code is originated from psa_res.R

# This file extract results from psa results

source("Model_PLoSMed_rev/libraries.R")
source("Model_PLoSMed_rev/gen_helpers.R")
source("Model_PLoSMed_rev/res_helpers.R")
source("Model_PLoSMed_rev/plot_helpers.R")

CURRENT_RESDIR = "/n/netscratch/menzies_lab/Everyone/ylh202/IncipienTB/PSA_RNAgam"

# =================================================================================================
# CHECK IF WE HAVE ALL THE SIMULATION DATA
# =================================================================================================


# Check whether all parameter sets were simulated
dir_list = list.dirs(path = paste0(CURRENT_RESDIR),recursive = FALSE)
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,25, 28)))
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,79, 82)))
ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, 76)))

if (length(ptemp) > 0) stop("Some parameter sets were not simulated.")

# Check whether all countries were simulated for each data set 
# get_file_list = function(dir, type){
#   file_list = list.files(path = paste0(dir), pattern = paste0("\\.",type,"$"))
#   return(length(file_list))
# }

get_file_list = function(dir){
  file_list = list.files(path = paste0(dir), pattern = ".RDS")
  return(length(file_list))
}


library(future.apply)
# ctemp = future_sapply(dir_list, function(x) get_file_list(x, "RDS"), future.seed = NULL)
ctemp = future_sapply(dir_list, function(x) get_file_list(x), future.seed = NULL)

etemp = which(ctemp != 99)
# which paramsets are not simulated completely. 
if(any(length(etemp) != 0))stop("Check simulation.")

resim = str_sub(names(etemp), 79, 82)

# temp = as.numeric(str_sub(dir_list,25, 28))

# exclude = temp[etemp]
# v.paramset = temp[-which(temp %in% exclude)]
# n.paramset = length(v.paramset)


# =================================================================================================
# LOAD NECESSARY DAT
# =================================================================================================


# load countrynames and tbinc data.
PSA = F
source('Data/data.R')
cohort2019_by_coo = cohort2019_by_coo_age %>%
  group_by(place_of_birth) %>%
  summarise(total_pop = sum(pop_size_est),
            mean_age = weighted.mean(entry_age,pop_size_est),
            total_tb_cases = sum(tb_cases_per_person*pop_size_est),
            tb_inc = weighted.mean(tb_cases_per_person, pop_size_est),
            ltbi_calib = weighted.mean(ltbi_calib, pop_size_est),
            ltbi_calib_lower = weighted.mean(p025_calib, pop_size_est),
            ltbi_calib_upper = weighted.mean(p975_calib, pop_size_est))
# load WHO data on TB incidence (per 100k) in origin countries
who = read.csv("Data/TB_burden_countries_2022-09-22.csv")
who_subset = who %>%
  filter(year == 2019, iso3 %in% cohort2019_by_coo$place_of_birth) %>%
  select(country, iso3, year, e_pop_num,
         e_inc_100k, e_inc_100k_lo, e_inc_100k_hi,
         e_inc_num,	e_inc_num_lo,	e_inc_num_hi)
who_subset2 = who_subset %>% select(iso3, e_inc_100k) %>% rename(., place_of_birth = iso3)
rm(who)


# WHO dataset doesn't have TWN, YUG, SUN
# TWN: 37/100k https://www.cdc.gov.tw/En/Category/ListContent/bg0g_VU_Ysrgkes_KRUDgQ?uaid=0WialNbsh7SEGERJLa29FA
# YUG: now MKD: 12/100k
# SUN:


# include the WHO estimate in cohort2019_by_coo dataframe
twn = cohort2019_by_coo %>% filter(place_of_birth == "TWN") %>% mutate(e_inc_100k = 37.0) # pull this row out for now
cohort2019_by_coo2 = cohort2019_by_coo %>% filter(!(place_of_birth %in% c("SUN", "YUG"))) # drop obsolete countries
cohort2019_by_coo2 = merge(cohort2019_by_coo2,who_subset2,by= "place_of_birth")
cohort2019_by_coo2 = rbind(cohort2019_by_coo2, twn)

countrynames = cohort2019_by_coo2$place_of_birth[order(cohort2019_by_coo2$place_of_birth)]

# =================================================================================================
# METHODS: BASELINE CHARACTERISTICS
# =================================================================================================

# Rank countries based on 2019 WHO estimates of country TB incidence
coo_by_TBinc = cohort2019_by_coo2[order(-cohort2019_by_coo2$e_inc_100k),]

coo_by_TBinc = coo_by_TBinc %>%
  mutate(risk_category = case_when(e_inc_100k >= 300 ~ "I_300plus", 
                                   e_inc_100k >= 100 & e_inc_100k < 300 ~ "II_100_300",
                                   e_inc_100k >= 10  & e_inc_100k < 100 ~ "III_10_100",
                                   e_inc_100k < 10                      ~ "VI_0_10")) %>%
  mutate(`300+` = if_else(e_inc_100k >= 300, 1, 0),                       # highly endemic        
         `100_299` = if_else(e_inc_100k >= 100 & e_inc_100k < 300, 1, 0), # endemic
         `10-99` = if_else(e_inc_100k >= 10 & e_inc_100k < 100, 1, 0),    # moderate
         `0-9`   = if_else(e_inc_100k < 10, 1, 0))                        # low

coo_by_riskcat = coo_by_TBinc %>% group_by(risk_category) %>%
  summarise(pop_size = sum(total_pop),
            entry_age.mean = weighted.mean(mean_age, total_pop),
            ltbi_prev = weighted.mean(ltbi_calib, total_pop),
            ltbi_prev_lower = weighted.mean(ltbi_calib_lower, total_pop),
            ltbi_prev_upper = weighted.mean(ltbi_calib_upper, total_pop),
            tb_inc = weighted.mean(tb_inc, total_pop)) %>%
  mutate(across(-risk_category, num, digits = 2))

riskcat_I = coo_by_TBinc[coo_by_TBinc$risk_category == "I_300plus",]$place_of_birth
riskcat_II = coo_by_TBinc[coo_by_TBinc$risk_category == "II_100_300",]$place_of_birth
riskcat_III = coo_by_TBinc[coo_by_TBinc$risk_category == "III_10_100",]$place_of_birth
riskcat_IV = coo_by_TBinc[coo_by_TBinc$risk_category == "VI_0_10",]$place_of_birth

coo_by_riskcat_tot = coo_by_TBinc %>% 
  summarise(pop_size = sum(total_pop),
            entry_age.mean = weighted.mean(mean_age, total_pop),
            ltbi_prev = weighted.mean(ltbi_calib, total_pop),
            ltbi_prev_lower = weighted.mean(ltbi_calib_lower, total_pop),
            ltbi_prev_upper = weighted.mean(ltbi_calib_upper, total_pop),
            tb_inc = weighted.mean(tb_inc, total_pop))

# =================================================================================================
# RESULTS 5 Net Monetary Benefit - rev
# =================================================================================================


eff_by_coo_risk1 = readRDS(paste0(CURRENT_RESDIR,"/eff_by_coo_risk1.RDS"))
prod_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/prod_by_coo_risk1.RDS"))
nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_costItems15.RDS"))


get_riskcat_nmb_paramset = function(A, B, C, wtp){
  
  A_array = array(unlist(A), dim = c(3, 5, 1000))
  
  z = C[2:4,]
  
  x= apply(A_array, 3, function(x) rowSums(x[,2:5]))
  y = B[2:4,]
  nmbSOC = z*wtp - (x - y)
  # ciSOC = round(apply((z*150000 - (x - y)), 1, quantile, c(0.025, 0.975)), digits = 1)
  
  x2= apply(A_array, 3, function(x) rowSums(x[,2:3]))
  nmbHC = z*wtp - x2
  # ciHC = round(apply((z*150000 - x2), 1, quantile, c(0.025, 0.975)), digits=1)
  
  return(list(nmbSOC = nmbSOC, nmbHC = nmbHC))
}



get_inmb_byCost = function(cost, wtp){
  
  # if (cost == 30) {
  #   nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_", cost, ".RDS"))
  #   # nmb_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_", cost, "_rev.RDS"))
  #   # nmb_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_", cost, "_rev.RDS"))
  #   # nmb_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_", cost, "_rev.RDS"))
  # } else 
  if (cost == 300) {
    nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_costItems.RDS"))
    # nmb_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_costItems_rev.RDS"))
    # nmb_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_costItems_rev.RDS"))
    # nmb_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_costItems_rev.RDS"))
  } else {
    nmb_by_coo_itm_risk1 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_costItems", cost, ".RDS"))
    # nmb_by_coo_itm_risk2 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_costItems", cost, "_rev.RDS"))
    # nmb_by_coo_itm_risk3 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_costItems", cost, "_rev.RDS"))
    # nmb_by_coo_itm_risk4 = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_costItems", cost, "_rev.RDS"))
  } 
  
  risk1 = get_riskcat_nmb_paramset(A = nmb_by_coo_itm_risk1, B = prod_by_coo_itm_risk1, C = eff_by_coo_risk1, wtp = wtp)
  # risk2 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk2, B = prod_by_coo_itm_risk2, C = eff_by_coo_risk2, wtp = wtp)
  # risk3 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk3, B = prod_by_coo_itm_risk3, C = eff_by_coo_risk3, wtp = wtp)
  # risk4 = get_riskcat_nmb_mean(A = nmb_by_coo_itm_risk4, B = prod_by_coo_itm_risk4, C = eff_by_coo_risk4, wtp = wtp)
  
  # wtp_text = format(wtp, big.mark=",", scientific = F)
  # df = data.frame(mean = c(risk1, risk2, risk3, risk4),
  #                 RiskCat = rep(c("I", "II", "III", "IV"), each = 6),
  #                 Strategy = rep(rep(c("IGRA_TB", "IGRA_HrTS_TB", "HrTS_TB"), 2), 4), 
  #                 Perspective = rep(rep(c("SOC", "HC"), each = 3), 4),
  #                 Cost = rep(paste0("$",cost), 24),
  #                 WTP = rep(wtp_text, 24))
  
  return(risk1)
}

out = get_inmb_byCost(15, 150000)

PSAtable = read.csv("Data/PSAtable1000_gam.csv")
nmb_by_riskcat1_HC = data.frame(YEARS = PSAtable$RNA_TIME_upper,
                             SENS1 = PSAtable$RNA_SENS1,
                             SPEC  = PSAtable$RNA_SPEC,
                             iNMB.IGRA_TB = as.numeric(out$nmbHC[1,]),
                             iNMB.IGRA_RNA_TB = as.numeric(out$nmbHC[2,]),
                             iNMB.RNA_TB = as.numeric(out$nmbHC[3,]))

nmb_by_riskcat1_SOC = data.frame(YEARS = PSAtable$RNA_TIME_upper,
                                SENS1 = PSAtable$RNA_SENS1,
                                SPEC  = PSAtable$RNA_SPEC,
                                iNMB.IGRA_TB = as.numeric(out$nmbSOC[1,]),
                                iNMB.IGRA_RNA_TB = as.numeric(out$nmbSOC[2,]),
                                iNMB.RNA_TB = as.numeric(out$nmbSOC[3,]))

nmb_by_riskcat1_HC$optimal = NA
for(i in 1:nrow(nmb_by_riskcat1_HC)){
  nmb_by_riskcat1_HC$optimal[i] = which.max(nmb_by_riskcat1_HC[i, 4:6])
    
}
get_plot_diff_df_rev = function(dat){
  
  dat$diff = dat$iNMB.RNA_TB - dat$iNMB.IGRA_TB

  return(dat)
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
    stat_contour_filled(breaks= seq(-200, 50, 5), show.legend = F) + # c(min(df$inmb) + (0:10)*(900))) + 
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
                                      top = text_grob(title,face = "bold", size = 18, vjust = -0.5)) + 
    theme(plot.margin = margin(15, 5, 15, 5))
  option1_annotated
  
  return(option1_annotated)
  
}

ct_df_risk1_s4_hc_15 = get_plot_diff_df_rev(nmb_by_riskcat1_HC)
plot_fig(dat = ct_df_risk1_s4_hc_15, "Risk Category: > 300, Cost of HrTS: $15") #A9_Fig2-1-2_15  1070x 600 1250 x 680

ct_df_risk1_s4_soc_15 = get_plot_diff_df_rev(nmb_by_riskcat1_SOC)
plot_fig(ct_df_risk1_s4_soc_15, "Risk Category: > 300, Cost of HrTS: $15") #A9_Fig2-1-1_15 

# =================================================================================================
# RESULTS 5 Net Monetary Benefit
# =================================================================================================

get_nmb_by_coo = function(rdsfile){
  country_nmb_res = readRDS(rdsfile)
  country_nmb_res = unlist(country_nmb_res)
  return(country_nmb_res)
}

get_nmb_by_coo_df = function(paramsetnum){
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rdsfiles = list.files(paramset_dir, pattern = ".RDS")
  rdsfiles = rdsfiles[-which(rdsfiles %in% c("cea_YUG.RDS", "cea_SUN.RDS"))]
  temp = future_lapply(rdsfiles, function(x) get_nmb_by_coo(paste0(paramset_dir,x)))
  temp = do.call("rbind", temp)
  nmb_by_coo_array[,,paramsetnum] = temp
  
  return(nmb_by_coo_array)
}

# nmb_by_coo_array = array(NA, dim = c(97, 36, 1000))
# for (i in 1:1000){
#   nmb_by_coo_array = get_nmb_by_coo_df(i)
# }
# saveRDS(nmb_by_coo_array, file = paste0(CURRENT_RESDIR,"/nmb_by_coo_array.RDS"))

nmb_by_coo_array = readRDS(paste0(CURRENT_RESDIR,"/nmb_by_coo_array.RDS"))


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

# write.csv(nmb_by_riskcat, paste0(CURRENT_RESDIR,"/nmb_by_riskcat.csv"), row.names = F)

nmb_by_riskcat = read.csv(paste0(CURRENT_RESDIR,"/nmb_by_riskcat.csv"))

# ********************************* 
nmb_by_riskcat1 = nmb_by_riskcat %>% filter(risk_category == "I_300plus")
nmb_by_riskcat2 = nmb_by_riskcat %>% filter(risk_category == "II_100_300")
nmb_by_riskcat3 = nmb_by_riskcat %>% filter(risk_category == "III_10_100")
nmb_by_riskcat4 = nmb_by_riskcat %>% filter(risk_category == "VI_0_10")

# fit gam to the raw nmb data 
# simulate nmb data for full range of sens1, spec, time

# a function to get sens x spec for each time x cost combinations
get_optimalstrat_plot300 = function(time, dat){
  
  temp = dat %>% filter(TIME < time[2] & TIME >= time[1])%>%  as.data.frame %>% dplyr::select(NMB.No_screening, NMB.IGRA_TB, NMB.IGRA_RNA_TB, NMB.RNA_TB, SENS1, SPEC)
  out3 = mgcv::gam(NMB.RNA_TB ~ te(SENS1, SPEC), data=temp)
  out2 = mgcv::gam(NMB.IGRA_RNA_TB ~ te( SENS1, SPEC), data=temp)
  out1 = mgcv::gam(NMB.IGRA_TB ~ te( SENS1, SPEC), data=temp)
  out0 = mgcv::gam(NMB.No_screening ~ te(SENS1, SPEC), data=temp)
  
  dat.new = data.frame(SENS1 = rep(seq(0.75, 1.0, 0.01), each = 26),
                       SPEC = rep(seq(0.75, 1.0, 0.01), 26))
  pred0 = mgcv::predict.gam(out0, newdata = dat.new)
  pred1 = mgcv::predict.gam(out1, newdata = dat.new)
  pred2 = mgcv::predict.gam(out2, newdata = dat.new)
  pred3 = mgcv::predict.gam(out3, newdata = dat.new)
  
  dat.new$pred0 = pred0
  dat.new$pred1 = pred1
  dat.new$pred2 = pred2
  dat.new$pred3 = pred3
  
  dat.new$optimal = NA
  for(i in 1:nrow(dat.new)){
    dat.new$optimal[i]=which.max(c(dat.new$pred0[i], dat.new$pred1[i], dat.new$pred2[i], dat.new$pred3[i])) - 1
  }
  dat.new$optimal = as.factor(dat.new$optimal)
  
  # plot = ggplot(data = dat.new, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal))+  geom_tile(alpha = 0.9) + 
  #   xlab("Specificity") + 
  #   ylab("Sensitivity")
  
  dat.new$time = time[2]
  dat.new$RNAcost = 300
  
  # return(plot)
  return(dat.new)
}

get_optimalstrat_plot150 = function(time, dat){
  
  temp = dat %>% filter(TIME < time[2] & TIME >= time[1])%>%  as.data.frame %>% dplyr::select(NMB.No_screening, NMB.IGRA_TB, NMB.IGRA_RNA_TB150, NMB.RNA_TB150, SENS1, SPEC)
  out3 = mgcv::gam(NMB.RNA_TB150 ~ te(SENS1, SPEC), data=temp)
  out2 = mgcv::gam(NMB.IGRA_RNA_TB150 ~ te( SENS1, SPEC), data=temp)
  out1 = mgcv::gam(NMB.IGRA_TB ~ te( SENS1, SPEC), data=temp)
  out0 = mgcv::gam(NMB.No_screening ~ te(SENS1, SPEC), data=temp)
  
  dat.new = data.frame(SENS1 = rep(seq(0.75, 1.0, 0.01), each = 26),
                       SPEC = rep(seq(0.75, 1.0, 0.01), 26))
  pred0 = mgcv::predict.gam(out0, newdata = dat.new)
  pred1 = mgcv::predict.gam(out1, newdata = dat.new)
  pred2 = mgcv::predict.gam(out2, newdata = dat.new)
  pred3 = mgcv::predict.gam(out3, newdata = dat.new)
  
  dat.new$pred0 = pred0
  dat.new$pred1 = pred1
  dat.new$pred2 = pred2
  dat.new$pred3 = pred3
  
  dat.new$optimal = NA
  for(i in 1:nrow(dat.new)){
    dat.new$optimal[i]=which.max(c(dat.new$pred0[i], dat.new$pred1[i], dat.new$pred2[i], dat.new$pred3[i])) - 1
  }
  dat.new$optimal = as.factor(dat.new$optimal)
  
  # plot = ggplot(data = dat.new, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal))+  geom_tile(alpha = 0.9) + 
  #   xlab("Specificity") + 
  #   ylab("Sensitivity")
  
  dat.new$time = time[2]
  dat.new$RNAcost = 150
  
  # return(plot)
  return(dat.new)
}

get_optimalstrat_plot60 = function(time, dat){
  
  temp = dat %>% filter(TIME < time[2] & TIME >= time[1])%>%  as.data.frame %>% dplyr::select(NMB.No_screening, NMB.IGRA_TB, NMB.IGRA_RNA_TB60, NMB.RNA_TB60, SENS1, SPEC)
  out3 = mgcv::gam(NMB.RNA_TB60 ~ te(SENS1, SPEC), data=temp)
  out2 = mgcv::gam(NMB.IGRA_RNA_TB60 ~ te( SENS1, SPEC), data=temp)
  out1 = mgcv::gam(NMB.IGRA_TB ~ te( SENS1, SPEC), data=temp)
  out0 = mgcv::gam(NMB.No_screening ~ te(SENS1, SPEC), data=temp)
  
  dat.new = data.frame(SENS1 = rep(seq(0.75, 1.0, 0.01), each = 26),
                       SPEC = rep(seq(0.75, 1.0, 0.01), 26))
  pred0 = mgcv::predict.gam(out0, newdata = dat.new)
  pred1 = mgcv::predict.gam(out1, newdata = dat.new)
  pred2 = mgcv::predict.gam(out2, newdata = dat.new)
  pred3 = mgcv::predict.gam(out3, newdata = dat.new)
  
  dat.new$pred0 = pred0
  dat.new$pred1 = pred1
  dat.new$pred2 = pred2
  dat.new$pred3 = pred3
  
  dat.new$optimal = NA
  for(i in 1:nrow(dat.new)){
    dat.new$optimal[i]=which.max(c(dat.new$pred0[i], dat.new$pred1[i], dat.new$pred2[i], dat.new$pred3[i])) - 1
  }
  dat.new$optimal = as.factor(dat.new$optimal)
  
  # plot = ggplot(data = dat.new, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal))+  geom_tile(alpha = 0.9) + 
  #   xlab("Specificity") + 
  #   ylab("Sensitivity")
  
  dat.new$time = time[2]
  dat.new$RNAcost = 60
  
  # return(plot)
  return(dat.new)
}

get_optimalstrat_plot30 = function(time, dat){
  
  temp = dat %>% filter(TIME < time[2] & TIME >= time[1])%>%  as.data.frame %>% dplyr::select(NMB.No_screening, NMB.IGRA_TB, NMB.IGRA_RNA_TB30, NMB.RNA_TB30, SENS1, SPEC)
  out3 = mgcv::gam(NMB.RNA_TB30 ~ te(SENS1, SPEC), data=temp)
  out2 = mgcv::gam(NMB.IGRA_RNA_TB30 ~ te( SENS1, SPEC), data=temp)
  out1 = mgcv::gam(NMB.IGRA_TB ~ te( SENS1, SPEC), data=temp)
  out0 = mgcv::gam(NMB.No_screening ~ te(SENS1, SPEC), data=temp)
  
  dat.new = data.frame(SENS1 = rep(seq(0.75, 1.0, 0.01), each = 26),
                       SPEC = rep(seq(0.75, 1.0, 0.01), 26))
  pred0 = mgcv::predict.gam(out0, newdata = dat.new)
  pred1 = mgcv::predict.gam(out1, newdata = dat.new)
  pred2 = mgcv::predict.gam(out2, newdata = dat.new)
  pred3 = mgcv::predict.gam(out3, newdata = dat.new)
  
  dat.new$pred0 = pred0
  dat.new$pred1 = pred1
  dat.new$pred2 = pred2
  dat.new$pred3 = pred3
  
  dat.new$optimal = NA
  for(i in 1:nrow(dat.new)){
    dat.new$optimal[i]=which.max(c(dat.new$pred0[i], dat.new$pred1[i], dat.new$pred2[i], dat.new$pred3[i])) - 1
  }
  dat.new$optimal = as.factor(dat.new$optimal)
  
  # plot = ggplot(data = dat.new, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal))+  geom_tile(alpha = 0.9) + 
  #   xlab("Specificity") + 
  #   ylab("Sensitivity")
  
  dat.new$time = time[2]
  dat.new$RNAcost = 30
  
  # return(plot)
  return(dat.new)
}

get_optimalstrat_plot15 = function(time, dat){
  
  temp = dat %>% filter(TIME < time[2] & TIME >= time[1])%>%  as.data.frame %>% dplyr::select(NMB.No_screening, NMB.IGRA_TB, NMB.IGRA_RNA_TB15, NMB.RNA_TB15, SENS1, SPEC)
  out3 = mgcv::gam(NMB.RNA_TB15 ~ te(SENS1, SPEC), data=temp)
  out2 = mgcv::gam(NMB.IGRA_RNA_TB15 ~ te( SENS1, SPEC), data=temp)
  out1 = mgcv::gam(NMB.IGRA_TB ~ te( SENS1, SPEC), data=temp)
  out0 = mgcv::gam(NMB.No_screening ~ te(SENS1, SPEC), data=temp)
  
  dat.new = data.frame(SENS1 = rep(seq(0.75, 1.0, 0.01), each = 26),
                       SPEC = rep(seq(0.75, 1.0, 0.01), 26))
  pred0 = mgcv::predict.gam(out0, newdata = dat.new)
  pred1 = mgcv::predict.gam(out1, newdata = dat.new)
  pred2 = mgcv::predict.gam(out2, newdata = dat.new)
  pred3 = mgcv::predict.gam(out3, newdata = dat.new)
  
  dat.new$pred0 = pred0
  dat.new$pred1 = pred1
  dat.new$pred2 = pred2
  dat.new$pred3 = pred3
  
  dat.new$optimal = NA
  for(i in 1:nrow(dat.new)){
    dat.new$optimal[i]=which.max(c(dat.new$pred0[i], dat.new$pred1[i], dat.new$pred2[i], dat.new$pred3[i])) - 1
  }
  dat.new$optimal = as.factor(dat.new$optimal)
  
  # plot = ggplot(data = dat.new, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal))+  geom_tile(alpha = 0.9) + 
  #   xlab("Specificity") + 
  #   ylab("Sensitivity")
  
  dat.new$time = time[2]
  dat.new$RNAcost = 15
  
  # return(plot)
  return(dat.new)
}

timeVec = list(c(2,3), c(3,4), c(4,5), c(5, 6), c(6, 7), c(7, 8), c(8,9), c(9,10))

# get dataframes for plotting 
get300_riskcat1 = lapply(timeVec, function(x) get_optimalstrat_plot300(x, nmb_by_riskcat1)) %>% do.call("rbind", .)
get150_riskcat1 = lapply(timeVec, function(x) get_optimalstrat_plot150(x, nmb_by_riskcat1)) %>% do.call("rbind", .)
get60_riskcat1  = lapply(timeVec, function(x) get_optimalstrat_plot60(x, nmb_by_riskcat1))  %>% do.call("rbind", .)
get30_riskcat1  = lapply(timeVec, function(x) get_optimalstrat_plot30(x, nmb_by_riskcat1))  %>% do.call("rbind", .)
get15_riskcat1  = lapply(timeVec, function(x) get_optimalstrat_plot15(x, nmb_by_riskcat1))  %>% do.call("rbind", .)


df_riskcat1 = rbind(get300_riskcat1, get150_riskcat1, get60_riskcat1, get30_riskcat1, get15_riskcat1)

get300_riskcat2 = lapply(timeVec, function(x) get_optimalstrat_plot300(x, nmb_by_riskcat2)) %>% do.call("rbind", .)
get150_riskcat2 = lapply(timeVec, function(x) get_optimalstrat_plot150(x, nmb_by_riskcat2)) %>% do.call("rbind", .)
get60_riskcat2  = lapply(timeVec, function(x) get_optimalstrat_plot60(x, nmb_by_riskcat2))  %>% do.call("rbind", .)
get30_riskcat2  = lapply(timeVec, function(x) get_optimalstrat_plot30(x, nmb_by_riskcat2))  %>% do.call("rbind", .)
get15_riskcat2  = lapply(timeVec, function(x) get_optimalstrat_plot15(x, nmb_by_riskcat2))  %>% do.call("rbind", .)

df_riskcat2 = rbind(get300_riskcat2, get150_riskcat2, get60_riskcat2, get30_riskcat2, get15_riskcat2)

get300_riskcat3 = lapply(timeVec, function(x) get_optimalstrat_plot300(x, nmb_by_riskcat3)) %>% do.call("rbind", .)
get150_riskcat3 = lapply(timeVec, function(x) get_optimalstrat_plot150(x, nmb_by_riskcat3)) %>% do.call("rbind", .)
get60_riskcat3  = lapply(timeVec, function(x) get_optimalstrat_plot60(x, nmb_by_riskcat3))  %>% do.call("rbind", .)
get30_riskcat3  = lapply(timeVec, function(x) get_optimalstrat_plot30(x, nmb_by_riskcat3))  %>% do.call("rbind", .)
get15_riskcat3  = lapply(timeVec, function(x) get_optimalstrat_plot15(x, nmb_by_riskcat3))  %>% do.call("rbind", .)

df_riskcat3 = rbind(get300_riskcat3, get150_riskcat3, get60_riskcat3, get30_riskcat3, get15_riskcat3)

get300_riskcat4 = lapply(timeVec, function(x) get_optimalstrat_plot300(x, nmb_by_riskcat4)) %>% do.call("rbind", .)
get150_riskcat4 = lapply(timeVec, function(x) get_optimalstrat_plot150(x, nmb_by_riskcat4)) %>% do.call("rbind", .)
get60_riskcat4  = lapply(timeVec, function(x) get_optimalstrat_plot60(x, nmb_by_riskcat4))  %>% do.call("rbind", .)
get30_riskcat4  = lapply(timeVec, function(x) get_optimalstrat_plot30(x, nmb_by_riskcat4))  %>% do.call("rbind", .)
get15_riskcat4  = lapply(timeVec, function(x) get_optimalstrat_plot15(x, nmb_by_riskcat4))  %>% do.call("rbind", .)

df_riskcat4 = rbind(get300_riskcat4, get150_riskcat4, get60_riskcat4, get30_riskcat4, get15_riskcat4)

# new facet label names for predictive time frame 
time.labs <- c("Yr 0-3", "Yr 0-4", "Yr 0-5", "Yr 0-6", "Yr 0-7", "Yr 0-8", "Yr 0-9", "Yr 0-10")
names(time.labs) <- c("3", "4", "5", "6", "7", "8", "9", "10")

# new facet label names for supp variable
cost.labs <- c("$30", "$60", "$150", "$300")
names(cost.labs) <- c("30", "60", "150", "300")

ggplot(data = df_riskcat1, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal)) +  
  geom_tile(alpha = 1) + 
  scale_x_discrete(breaks = seq(0.75,1, 0.05))+
  scale_y_discrete(breaks = seq(0.75,1, 0.05))+
  labs(title = c("Risk Category I (\u2265 300)"),
       x = "Specificity for TB", y = "Sensitivity for incipient TB") + 
  facet_grid(RNAcost ~ time,
             labeller = labeller(RNAcost = cost.labs, time = time.labs)) + 
 scale_fill_manual(name = "Optimal Strategy",
                   labels = c(  "II: IGRA-TB", "III: IGRA-RNA-TB"),
                   values = c("grey", "#619CFF")) + theme(strip.background = element_rect(fill="#efdc75"))

  
ggplot(data = df_riskcat2, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal)) +  
  geom_tile(alpha = 1) + 
  scale_x_discrete(breaks = seq(0.75,1, 0.05))+
  scale_y_discrete(breaks = seq(0.75,1, 0.05))+
  labs(title = c("Risk Category II (100-300)"),
       x = "Specificity for TB", y = "Sensitivity for incipient TB") + 
  facet_grid(RNAcost ~ time,
             labeller = labeller(RNAcost = cost.labs, time = time.labs)) + 
  scale_fill_manual(name = "Optimal Strategy",
                    labels = c(  "II: IGRA-TB", "III: IGRA-RNA-TB"),
                    values = c("grey", "#619CFF")) + theme(strip.background = element_rect(fill="#efdc75"))

ggplot(data = df_riskcat3, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal)) +  
  geom_tile(alpha = 1) + 
  scale_x_discrete(breaks = seq(0.75,1, 0.05))+
  scale_y_discrete(breaks = seq(0.75,1, 0.05))+
  labs(title = c("Risk Category III (10-100)"),
       x = "Specificity for TB", y = "Sensitivity for incipient TB") + 
  facet_grid(RNAcost ~ time,
             labeller = labeller(RNAcost = cost.labs, time = time.labs)) + 
  scale_fill_manual(name = "Optimal Strategy",
                    labels = c(  "II: IGRA-TB", "III: IGRA-RNA-TB"),
                    values = c("grey", "#619CFF")) + theme(strip.background = element_rect(fill="#efdc75"))

ggplot(data = df_riskcat4, aes(y = as.factor(SENS1), x = as.factor(SPEC), fill = optimal)) +  
  geom_tile(alpha = 1) + 
  scale_x_discrete(breaks = seq(0.75,1, 0.05))+
  scale_y_discrete(breaks = seq(0.75,1, 0.05))+
  labs(title = c("Risk Category IV (0-10)"),
       x = "Specificity for TB", y = "Sensitivity for incipient TB") + 
  facet_grid(RNAcost ~ time,
             labeller = labeller(RNAcost = cost.labs, time = time.labs)) + 
  scale_fill_manual(name = "Optimal Strategy",
                    labels = c( "I: No Screening", "II: IGRA-TB", "III: IGRA-RNA-TB"),
                    values = c("white", "grey", "#619CFF")) + theme(strip.background = element_rect(fill="#efdc75")) #1500 x 650





# fit full model
out = mgcv::gam(x ~ s(SENS1, SPEC, TIME), family = binomial(link = "logit"), data = temp)
out = mgcv::gam(x ~ SENS1 + SPEC + TIME + s(SENS1, SPEC, TIME), data = temp)
# generate x and y data for plotting
xord = seq(0.75, 1, 0.001)
yord = seq(0.75, 1, 0.001)

dat.new = data.frame(SPEC = rep(seq(0.75, 1.0, 0.001), 251),
                     SENS1 = rep(seq(0.75, 1.0, 0.001), each = 251))
pred0 = predict(out, newdata = dat.new, type = "response")

dat.new$pred0 = pred0
pmatu  = matrix(rev(pred0), 251, 251)

temp = nmb_by_riskcat1
temp$x = temp$iNMB.IGRA_RNA_TB - temp$iNMB.IGRA_TB

out.rsm = rsm::rsm(s3.60 ~  FO(SENS1,SPEC)+TWI(SENS1, TIME) + +TWI(SENS1, TIME), data = temp)
out.rsm = rsm::rsm(x ~  FO(SENS1,SPEC, TIME) + TWI(SENS1,SPEC, TIME), data = temp)


contour(yord, xord, pmatu) #sens vs spec
contour(out.rsm,  ~ TWI(SENS1,SPEC), at = data.frame(TIME = 10), image = T)
contour(out.rsm,  ~ TWI(TIME,SPEC), at = data.frame(SENS1 = 1), image = T)
persp(out.rsm,  ~ TWI(TIME,SPEC), at = data.frame(SENS1 = 0.3), image = T)


persp(out.rsm, SPEC ~ TIME, at = data.frame(SPEC= 0.2))
persp(out.rsm, SPEC ~ SENS1, at = data.frame(TIME = 9))

persp(out.rsm, SENS1 ~ SPEC, at = data.frame(TIME = 10))


# countour plot

get_diff = function(dat){

    dat$diff.300 = dat$iNMB.IGRA_RNA_TB - dat$iNMB.IGRA_TB
    dat$diff.150 = dat$iNMB.IGRA_RNA_TB150 - dat$iNMB.IGRA_TB
    dat$diff.60  = dat$iNMB.IGRA_RNA_TB60 - dat$iNMB.IGRA_TB
    dat$diff.30  = dat$iNMB.IGRA_RNA_TB30 - dat$iNMB.IGRA_TB
    dat$diff.15  = dat$iNMB.IGRA_RNA_TB15 - dat$iNMB.IGRA_TB
  
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

get_persp_plot = function(ct_df, dat){
  a = lapply(3:10, function(x) persp(ct_df, SENS1 ~ SPEC, at = data.frame(YEARS = x),
                                     ~ TWI(SENS1,SPEC), at = data.frame(YEARS = x), image = T,
                                     xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
                                     # ylabs = c("P(T+ | incipient TB)"),
                                     zlab = "iNMB(III vs II)",
                                     # main = paste0("Predictive timeframe = ", x, " years"),
                                     xlim = c(0.75, 1), ylim = c(0.75, 1), zlim = c(min(ct_df$fitted.values)-500, max(ct_df$fitted.values)+500)))
  
}

get_contour = function(ct_df){
  par(mfrow=c(2,4),c(10, 10, 4.1, 2.1))
  
  # op = par(mfrow=c(2,4))
  # b = lapply(3:10, function(x) contour(ct_df, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
  #                                      xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
  #                                      xlim = c(0.75, 1), ylim = c(0.75, 1)))
  for (x in 3:10) { contour(ct_df, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
                                       xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
                                       xlim = c(0.75, 1), ylim = c(0.75, 1))}
  mtext("P(T+ | incipient TB)",side=1,line=0,outer=TRUE,cex=1.3)
  mtext("P(T- | no TB in lifetime)",side=1,line=0,outer=TRUE,cex=1.3, las=0)
  # title(xlab = "P(T+ | incipient TB)",
  #       ylab = "P(T- | no TB in lifetime)",
  #       outer = TRUE, line = 3)
  par(op)
}

ct_df_risk1 = get_contour_plot_df(nmb_by_riskcat1)
ct_df_risk2 = get_contour_plot_df(nmb_by_riskcat2)
ct_df_risk3 = get_contour_plot_df(nmb_by_riskcat3)
ct_df_risk4 = get_contour_plot_df(nmb_by_riskcat4)





plot_riskcat = function(dat){
op = par(mfrow=c(2,4),mar=c(4,4,2, 2), oma = c(5, 5, 3, 1), cex.axis= 1.2, cex.lab= 1.3)


# b = lapply(3:10, function(x) contour(ct_df, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
#                                      xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
#                                      xlim = c(0.75, 1), ylim = c(0.75, 1)))
b = for (x in 3:10) { contour(dat, ~ TWI(SENS1,SPEC), at = data.frame(TIME = x), image = T,
                         # xlabs = c("P(T+ | incipient TB)", "P(T- | no TB in lifetime)"),
                         xlabs = c("  ", " "),
                          xlim = c(0.75, 1), ylim = c(0.75, 1), 
                        # labcex = 0.75,
                         lwd = 1)} 
mtext("P(T+ | incipient TB)",side=1,line=2,outer=TRUE, at = 0.5)
mtext("P(T- | no TB in lifetime)",side=2,line=0,outer=TRUE, las=0, at = 0.5)
}


# WTP 150k/QALY, Risk Category 1
png(file=paste0(getwd(),"/contour_risk1_300.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[1]]) 
title(main = "Risk Category I, Cost = 300 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/contour_risk1_150.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[2]]) 
title(main = "Risk Category I, Cost = 150 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/contour_risk1_60.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[3]]) 
title(main = "Risk Category I, Cost = 60 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/contour_risk1_30.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[4]]) 
title(main = "Risk Category I, Cost = 30 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

png(file=paste0(getwd(),"/contour_risk1_15.png"), width=1070, height=600)
plot_riskcat(ct_df_risk1[[5]]) 
title(main = "Risk Category I, Cost = 15 USD", outer = TRUE,xpd=NA, cex.main = 2)
dev.off()

