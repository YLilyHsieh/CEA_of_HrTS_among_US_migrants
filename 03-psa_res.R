#!/usr/bin/env Rscript

SA = commandArgs(trailingOnly=TRUE)[1]
SA = as.character(SA)

# This file extract results from psa results

source("Model/libraries.R")
source("Model/gen_helpers.R")
source("Model/res_helpers.R")

library(future.apply)
plan(multisession)

CURRENT_RESDIR = paste0("/n/netscratch/menzies_lab/Everyone/ylh202/IncipienTB/",SA)



# =================================================================================================
# CHECK IF WE HAVE ALL THE SIMULATION DATA
# =================================================================================================


# # Check whether all parameter sets were simulated
# dir_list = list.dirs(path = paste0(CURRENT_RESDIR),recursive = FALSE)
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,71, nchar(dir_list)))) #PSA_main
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,74, nchar(dir_list)))) #PSA_noreinf
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, nchar(dir_list)))) #PSA_reinf2
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, nchar(dir_list)))) #PSA_RNAgam
# 
# if (length(ptemp) > 0) stop("Some parameter sets were not simulated.")
# 
# # Check whether all countries were simulated for each data set
# get_file_list = function(dir, type){
#   file_list = list.files(path = paste0(dir), pattern = paste0("\\.",type,"$"))
#   return(length(file_list))
# }
# 
# ctemp = future_sapply(dir_list, function(x) get_file_list(x, "RDS"), future.seed = NULL)
# etemp = which(ctemp != 99)
# # which paramsets are not simulated completely.
# if(any(length(etemp) != 0))stop("Check simulation.")
# resim = str_sub(names(etemp), 71, nchar(names(etemp)))



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
# YUG: country no longer exists
# SUN: country no longer exists

# include the WHO estimate in cohort2019_by_coo dataframe
twn = cohort2019_by_coo %>% filter(place_of_birth == "TWN") %>% mutate(e_inc_100k = 37.0) # pull this row out for now
cohort2019_by_coo2 = cohort2019_by_coo %>% filter(!(place_of_birth %in% c("SUN", "YUG"))) # drop obsolete countries
cohort2019_by_coo2 = merge(cohort2019_by_coo2,who_subset2,by= "place_of_birth")
cohort2019_by_coo2 = rbind(cohort2019_by_coo2, twn)


# =================================================================================================
# Useful functions
# =================================================================================================


extract_n_per_country = function(csvfile, strategynum, event){
  df = read.csv(csvfile)
  n = sum(df[df$Strategy == strategynum & df$Event %in%  event,]$n)
  n = ifelse( length(n) ==0, 0, n)
  return(n)
}

extract_n_per_country2 = function(csvfile, strategynum, event = "TREAT.TB.GUARANTEE"){
  df = read.csv(csvfile)
  n = sum(df[df$Strategy == strategynum & df$Event %in%  event & Time < POSTARRIVAL_SCREENING,]$n)
  n = ifelse( length(n) ==0, 0, n)
  return(n)
}

extract_n_per_paramset = function(paramsetnum, event){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*event_tab_"))
  csvfiles = csvfiles[-c(which(str_sub(csvfiles, 11,13) %in% c("YUG", "SUN")))]
  
  n1 = unlist(sapply(csvfiles, function(x) extract_n_per_country(paste0(paramset_dir,x), 0, event)))
  n2 = unlist(sapply(csvfiles, function(x) extract_n_per_country(paste0(paramset_dir,x), 1, event)))
  n3 = unlist(sapply(csvfiles, function(x) extract_n_per_country(paste0(paramset_dir,x), 2, event)))
  n4 = unlist(sapply(csvfiles, function(x) extract_n_per_country(paste0(paramset_dir,x), 3, event)))
  
  return(c(n1, n2, n3, n4))
}

extract_n_per_paramset2 = function(paramsetnum, event = "TREAT.TB.GUARANTEE"){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*event_tab_"))
  csvfiles = csvfiles[-c(which(str_sub(csvfiles, 11,13) %in% c("YUG", "SUN")))]
  
  n1 = unlist(sapply(csvfiles, function(x) extract_n_per_country2(paste0(paramset_dir,x), 0, event)))
  n2 = unlist(sapply(csvfiles, function(x) extract_n_per_country2(paste0(paramset_dir,x), 1, event)))
  n3 = unlist(sapply(csvfiles, function(x) extract_n_per_country2(paste0(paramset_dir,x), 2, event)))
  n4 = unlist(sapply(csvfiles, function(x) extract_n_per_country(paste0(paramset_dir,x), 3, event)))
  
  return(c(n1, n2, n3, n4))
}

get_abs_by_riskcat = function(dat){
  r1 = colSums(dat[dat$country %in% riskcat_I,1:4])
  r2 = colSums(dat[dat$country %in% riskcat_II, 1:4])
  r3 = colSums(dat[dat$country %in% riskcat_III, 1:4])
  r4 = colSums(dat[dat$country %in% riskcat_IV, 1:4])
  
  cat_abs = rbind(r1, r2, r3, r4)
  
  return(cat_abs)
}

pretty = function(x0,r, sep) {
  
  x = gsub(" ","",format(round(x0,r),nsmall=r,big.mark= ","))
  
  return(paste0(x[1]," (",x[2],sep,x[3],")"))
}

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



#& ---Table 1; Text --- &#
coo_by_riskcat_tot
coo_by_riskcat

riskcat_I
riskcat_II
riskcat_III
riskcat_IV
#& -------------------- &#



# =================================================================================================
# RESULTS: TB OUTCOMES
# =================================================================================================

# -------------------------------------------------------------------------------------------------
# TIME-TO-TB
# -------------------------------------------------------------------------------------------------

# Useful functions
read_time2TB_files = function(paramsetnum){
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*time2TB_"))
  csvfiles = csvfiles[-c(which(str_sub(csvfiles, 9,11) %in% c("YUG", "SUN")))]
  
  
  # read in data
  dat = lapply(csvfiles, function(x) read.csv(paste0(paramset_dir, x)))
  
  # rbind data
  dat = do.call("rbind", dat)
  
  # calculate country-level n
  w_df = data.frame(Year = rep(seq(0, 90),4),
                    Strategy = rep(c(0:3), each = 91),
                    cases = rep(NA, 91*4))
  
  
  for (j in 0:3){
    
    w = dat %>%
      filter(Time>0) %>%
      filter(Strategy == j) %>%
      mutate(Time = floor(Time)) %>%
      count(Time)
    
    for(i in 1:91){
      w_df$cases[j*91+i] = ifelse(w_df$Year[j*91+i] %in% w$Time, w$n[w$Time == w_df$Year[i]], 0)
    }
  }
  
  return(w_df$cases)
}

read_time2TB_files_all = function(paramsetnum){
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*time2TB_"))
  csvfiles = csvfiles[-c(which(str_sub(csvfiles, 9,11) %in% c("YUG", "SUN")))]
  
  
  # read in data
  dat = lapply(csvfiles, function(x) read.csv(paste0(paramset_dir, x)))
  
  # rbind data
  dat = do.call("rbind", dat)
  
  # calculate country-level n
  w_df = data.frame(Year = rep(seq(0, 90),4),
                    Strategy = rep(c(0:3), each = 91),
                    cases = rep(NA, 91*4))
  
  
  for (j in 0:3){
    
    w = dat %>%
      #filter(Time>0) %>%
      filter(Strategy == j) %>%
      mutate(Time = floor(Time)) %>%
      count(Time)
    
    for(i in 1:91){
      w_df$cases[j*91+i] = ifelse(w_df$Year[j*91+i] %in% w$Time, w$n[w$Time == w_df$Year[i]], 0)
    }
  }
  
  return(w_df$cases)
}


read_time2TB_files_raw = function(paramsetnum){
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*time2TB_"))
  csvfiles = csvfiles[-c(which(str_sub(csvfiles, 9,11) %in% c("YUG", "SUN")))]
  
  
  # read in data
  dat = lapply(csvfiles, function(x) read.csv(paste0(paramset_dir, x)))
  
  # rbind data
  dat = do.call("rbind", dat)
  
  # calculate summary stats
  dat2 = dat %>% filter(Strategy == 0) %>% select(Time)
  
  return(dat2)
}

read_time2TB_files_all_coo = function(paramsetnum, countries){
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*time2TB_"))
  csvfiles = csvfiles[c(which(str_sub(csvfiles, 9,11) %in% countries))]
  
  
  # read in data
  dat = lapply(csvfiles, function(x) read.csv(paste0(paramset_dir, x)))
  
  # rbind data
  dat = do.call("rbind", dat)
  
  # calculate country-level n
  w_df = data.frame(Year = rep(seq(0, 90),4),
                    Strategy = rep(c(0:3), each = 91),
                    cases = rep(NA, 91*4))
  
  
  for (j in 0:3){
    
    w = dat %>%
      #filter(Time>0) %>%
      filter(Strategy == j) %>%
      mutate(Time = floor(Time)) %>%
      count(Time)
    
    for(i in 1:91){
      w_df$cases[j*91+i] = ifelse(w_df$Year[j*91+i] %in% w$Time, w$n[w$Time == w_df$Year[i]], 0)
    }
  }
  
  return(w_df$cases)
}

time2TB_sum = function(x){
  # dat = x %>% filter(Strategy == 0)
  dat = data.frame(mean = mean(x$Time),
                   median = median(x$Time),
                   lower = quantile(x$Time, p = 0.025),
                   upper = quantile(x$Time, p = 0.975))
  return(dat)
}

num_cases_entry = function(x){
  n = sum(x$Time <= 1/12)
  return(n)
}

read_time2death_files_raw = function(paramsetnum){
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  csvfiles = list.files(paramset_dir, pattern = paste0("*time2death_"))
  csvfiles = csvfiles[-c(which(str_sub(csvfiles, 12,14) %in% c("YUG", "SUN")))]
  
  
  # read in data
  dat = lapply(csvfiles, function(x) read.csv(paste0(paramset_dir, x)))
  
  # rbind data
  dat = do.call("rbind", dat)
  
  # calculate summary stats
  dat2 = dat %>% filter(Strategy == 0) %>% select(Time)
  
  return(dat2)
}

# ------------------------------------------------------------------------------
#  Summary stats on time to TB
# ------------------------------------------------------------------------------

library(future.apply)
plan(multisession)

# time2TB_raw = future_lapply(1:1000, function(x) read_time2TB_files_raw(x))
# saveRDS(time2TB_raw, paste0(CURRENT_RESDIR,"/time2TB_raw.RDS"))
# rm(time2TB_raw)
# 
# time2TB_all = future_lapply(1:1000, function(x) read_time2TB_files_all(x))
# saveRDS(time2TB_all, paste0(CURRENT_RESDIR,"/time2TB_all.RDS"))
# rm(time2TB_all)
# 
# time2TB_inc = future_lapply(1:1000, function(x) read_time2TB_files(x))
# saveRDS(time2TB_inc, paste0(CURRENT_RESDIR,"/time2TB_inc.RDS"))
# rm(time2TB_inc)


time2TB_raw = readRDS(paste0(CURRENT_RESDIR,"/time2TB_raw.RDS"))
median = mean(sapply(1:1000, function(x) median(time2TB_raw[[x]]$Time))) # take the mean of the medians
lower  = mean(sapply(1:1000, function(x) quantile(time2TB_raw[[x]]$Time, p = 0.25)))
upper  = mean(sapply(1:1000, function(x) quantile(time2TB_raw[[x]]$Time, p = 0.75)))

#& --- text --- &#
median
lower
upper
# -------------- #


# ------------------------------------------------------------------------------
#  Summary stats on number of TB by year
# ------------------------------------------------------------------------------

# number of lifetime TB cases
cases = sapply(1:1000, function(x) nrow(time2TB_raw[[x]]))
mean = mean(cases)
ci = quantile(cases, p = c(0.025, 0.975))

#Summary stats on the number of cases that occurred prior to screening

nTB_entry = mclapply(c(1:1000), function(x) num_cases_entry(time2TB_raw[[x]]))
nTB_entry = do.call("rbind",nTB_entry)

time2TB_national = read.csv(paste0(CURRENT_RESDIR, "/time2TB_national.csv"))
tt = time2TB_national %>% filter(Strategy == 0)
total = colSums(tt)

mean(total[3:1002])
quantile(total[3:1002], p = c(0.025, 0.975))

# Summary stats on the number of cases in the highest riskcat
time2TB_all_risk1 = readRDS(paste0(CURRENT_RESDIR,"/time2TB_riskcat1.RDS"))
time2TB_riskcat1 = do.call("cbind", time2TB_all_risk1)
tbinc_first1yr_riskcat1 = mean(sapply(1:1000, function(x) sum(time2TB_riskcat1[1,x])))/c(coo_by_riskcat$pop_size[1])

# Time-to-TB plot

time2TB = readRDS(paste0(CURRENT_RESDIR,"/time2TB_all.RDS"))
# time2TB = do.call("cbind", time2TB)
# colnames(time2TB) = c(1:1000)
# time2TB_national = data.frame(Year = rep(seq(0, 90),4),
#                               Strategy = rep(c(0:3), each = 91),
#                               time2TB)
#
# write.csv(time2TB_national,paste0(CURRENT_RESDIR, "/time2TB_national.csv"), row.names = F)

time2TB_national = read.csv(paste0(CURRENT_RESDIR, "/time2TB_national.csv"))

time2TB_reshape = reshape2::melt(time2TB_national, id.vars = c("Year", "Strategy"))

time2TB_national.summary = data.frame(Year = time2TB_national$Year,
                                      Strategy = time2TB_national$Strategy,
                                      mean = rowMeans(time2TB_national[,3:1002]),
                                      median = apply(time2TB_national[,3:1002], 1, quantile, probs = 0.5),
                                      lower = apply(time2TB_national[,3:1002],1, quantile, probs = 0.025),
                                      upper = apply(time2TB_national[,3:1002],1, quantile, probs = 0.975))

prop_first2yrs = mean(sapply(3:1002, function(x) sum(time2TB_national[1:2,x])/sum(time2TB_national[1:91,x])))
prop_first2yrs_lower = quantile(sapply(3:1002, function(x) sum(time2TB_national[1:2,x])/sum(time2TB_national[1:91,x])), 0.025)
prop_first2yrs_upper = quantile(sapply(3:1002, function(x) sum(time2TB_national[1:2,x])/sum(time2TB_national[1:91,x])), 0.975)

prop_first5yrs = mean(sapply(3:1002, function(x) sum(time2TB_national[1:5,x])/sum(time2TB_national[,x])))
prop_first10yrs = mean(sapply(3:1002, function(x) sum(time2TB_national[1:10,x])/sum(time2TB_national[,x])))

prop_first1yrs = mean(sapply(3:1002, function(x) sum(time2TB_national[1,x])))
prop_first1yrs_lower = quantile(sapply(3:1002, function(x) sum(time2TB_national[1,x])/sum(time2TB_national[1:91,x])), 0.025)
prop_first1yrs_upper = quantile(sapply(3:1002, function(x) sum(time2TB_national[1,x])/sum(time2TB_national[1:91,x])), 0.975)


#& ---text discussion ---&#
time2TB_national$prop.preventable.tb = 0.5^(time2TB_national$Year/20)
time2TB_national.summary = data.frame(Year = time2TB_national$Year,
                                      Strategy = time2TB_national$Strategy,
                                      mean = rowMeans(time2TB_national[,3:1002])*time2TB_national$prop.preventable.tb,
                                      median = apply(time2TB_national[,3:1002], 1, quantile, probs = 0.5)*time2TB_national$prop.preventable.tb,
                                      lower = apply(time2TB_national[,3:1002], 1, quantile, probs = 0.025)*time2TB_national$prop.preventable.tb,
                                      upper = apply(time2TB_national[,3:1002], 1, quantile, probs = 0.975)*time2TB_national$prop.preventable.tb)

prop_first2yrs = mean(sapply(3:1002, function(x) sum(time2TB_national[1:2,x]*time2TB_national$prop.preventable.tb[1:2])/sum(time2TB_national[1:91,x])*time2TB_national$prop.preventable.tb[1:91]))
prop_first2yrs_lower = quantile(sapply(3:1002, function(x) sum(time2TB_national[1:2,x])/sum(time2TB_national[1:91,x])), 0.025)
prop_first2yrs_upper = quantile(sapply(3:1002, function(x) sum(time2TB_national[1:2,x])/sum(time2TB_national[1:91,x])), 0.975)

colSums(time2TB_national.summary[time2TB_national.summary$Strategy==0,])


#& ---text--- &#
colSums(time2TB_national.summary[time2TB_national.summary$Strategy==0,])
# ------------ #

# gam_time2TB_s0 = mgcv::gam(mean ~ s(Year) + I(Year == 0), family = poisson(link = log),
#                            data = time2TB_national.summary[time2TB_national.summary$Strategy == 0,])
#
# gam_time2TB_s0_sq = mgcv::gam(mean ~ s(sqrt(Year)) + I(sqrt(Year) == 0), family = poisson(link = log),
#                            data = time2TB_national.summary[time2TB_national.summary$Strategy == 0,])
#
# gam_time2TB = mgcv::gam(value ~ as.factor(Strategy) + s(Year) + I(Year == 0) + s(Year, by = as.factor(Strategy)) , family = poisson(link = log),
#                   data = time2TB_reshape)
#
# gam_time2TB_s0 = mgcv::gam(value ~ I(Year == 0) + s(Year) , family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 0,])
# gam_time2TB_s1 = mgcv::gam(value ~ I(Year == 0) + s(Year) , family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 1,])
# gam_time2TB_s2 = mgcv::gam(value ~ I(Year == 0) + s(Year) , family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 2,])
# gam_time2TB_s3 = mgcv::gam(value ~ I(Year == 0) + s(Year) , family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 3,])
#
# preds_s0 = tidymv::get_gam_predictions(gam_time2TB_s0, Year, transform = exp, series_length = 91)
# preds_s1 = tidymv::get_gam_predictions(gam_time2TB_s1, Year, transform = exp, series_length = 91)
# preds_s2 = tidymv::get_gam_predictions(gam_time2TB_s2, Year, transform = exp, series_length = 91)
# preds_s3 = tidymv::get_gam_predictions(gam_time2TB_s3, Year, transform = exp, series_length = 91)
#
#
# preds = rbind(preds_s0, preds_s1, preds_s2, preds_s3)
# preds$strategy = rep(c(0:3), each = nrow(preds_s0))
#
# ggplot() +
#   geom_line(aes(x = Year, y = value, col = as.factor(strategy)), data = preds, size = 0.75) +
#   geom_ribbon(aes(x = Year, ymin = CI_lower, ymax = CI_upper, fill = as.factor(strategy)), data = preds) +
#   # geom_point(aes(x = Year, y = mean, col = as.factor(Strategy)), data = time2TB_national.summary, size = 2, shape = 1) +
#   scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
#   scale_x_continuous(limits = c(0, 90), breaks = seq(0 , 90, 10)) +
#   xlab("Years since entry") +
#   ylab("Count") +
#   ggtitle("Number of incident TB by year since entry") +
#   scale_color_discrete(name = "Strategy", labels = c("I: No screening", "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"))  +
#   theme_line() +
#   guides(fill = "none") + theme(
#     legend.position = c(.95, .95),
#     legend.justification = c("right", "top"),
#     legend.box.just = "right",
#     legend.margin = margin(6, 6, 6, 6)
#   ) # 580x500
#
# ggplot() +
#   geom_line(aes(x = Year, y = value, col = as.factor(strategy)), data = preds, size = 0.75) +
#   geom_ribbon(aes(x = Year, ymin = CI_lower, ymax = CI_upper, fill = as.factor(strategy)), data = preds) +
#   # geom_point(aes(x = Year, y = mean, col = as.factor(Strategy)), data = time2TB_national.summary, size = 2, shape = 1) +
#   scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
#   scale_x_continuous(limits = c(0, 10), breaks = seq(0 , 10, 10)) +
#   xlab("Years since entry") +
#   ylab("Count") +
#   ggtitle("Number of incident TB by year since entry") +
#   scale_color_discrete(name = "Strategy", labels = c("I: No screening", "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"))  +
#   theme_line() +
#   guides(fill = "none") + theme(
#     legend.position = c(.95, .95),
#     legend.justification = c("right", "top"),
#     legend.box.just = "right",
#     legend.margin = margin(6, 6, 6, 6)
#   )
#
# gam_time2TB_s0_sq = mgcv::gam(value ~ I(Year == 0) + s(sqrt(Year)) , family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 0,])
# gam_time2TB_s1_sq = mgcv::gam(value ~ I(Year == 0) + s(sqrt(Year)) , family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 1,])
# gam_time2TB_s2_sq = mgcv::gam(value ~ I(Year == 0) + s(sqrt(Year)), family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 2,])
# gam_time2TB_s3_sq = mgcv::gam(value ~ I(Year == 0) + s(sqrt(Year)), family = poisson(link = log),
#                            data = time2TB_reshape[time2TB_reshape$Strategy == 3,])
# 
# saveRDS(gam_time2TB_s0_sq, paste0(CURRENT_RESDIR, "/gam_time2TB_s0_sq_tot.RDS"))
# saveRDS(gam_time2TB_s1_sq, paste0(CURRENT_RESDIR, "/gam_time2TB_s1_sq_tot.RDS"))
# saveRDS(gam_time2TB_s2_sq, paste0(CURRENT_RESDIR, "/gam_time2TB_s2_sq_tot.RDS"))
# saveRDS(gam_time2TB_s3_sq, paste0(CURRENT_RESDIR, "/gam_time2TB_s3_sq_tot.RDS"))

gam_time2TB_s0_sq = readRDS(paste0(CURRENT_RESDIR, "/gam_time2TB_s0_sq_tot.RDS"))
gam_time2TB_s1_sq = readRDS(paste0(CURRENT_RESDIR, "/gam_time2TB_s1_sq_tot.RDS"))
gam_time2TB_s2_sq = readRDS(paste0(CURRENT_RESDIR, "/gam_time2TB_s2_sq_tot.RDS"))
gam_time2TB_s3_sq = readRDS(paste0(CURRENT_RESDIR, "/gam_time2TB_s3_sq_tot.RDS"))

preds_s0_sq = tidymv::get_gam_predictions(gam_time2TB_s0_sq, sqrt(Year), transform = exp, series_length = 91)
preds_s1_sq = tidymv::get_gam_predictions(gam_time2TB_s1_sq, sqrt(Year), transform = exp, series_length = 91)
preds_s2_sq = tidymv::get_gam_predictions(gam_time2TB_s2_sq, sqrt(Year), transform = exp, series_length = 91)
preds_s3_sq = tidymv::get_gam_predictions(gam_time2TB_s3_sq, sqrt(Year), transform = exp, series_length = 91)

preds_sq = rbind(preds_s0_sq, preds_s1_sq, preds_s2_sq, preds_s3_sq)
preds_sq$strategy = rep(c(0:3), each = nrow(preds_s0_sq))

main = ggplot() +
  geom_line(aes(x = Year + 0.5, y = value, col = as.factor(strategy)), data = preds_sq, size = 0.7) +
  geom_ribbon(aes(x = Year + 0.5, ymin = CI_lower, ymax = CI_upper, fill = as.factor(strategy)), alpha = 0.3, data = preds_sq) +
  geom_point(aes(x = Year + 0.5, y = mean, col = as.factor(Strategy)), data = time2TB_national.summary, size = 2, shape = 1) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 100)) +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0 , 90, 10)) +
  xlab("Years since entry") +
  ylab("Number of cases") +
  # ggtitle("Number of incident TB by year since entry") +
  scale_color_manual(name = "Strategy", labels = c("I: No screening", "II: IGRA only", "III: IGRA-HrTS", "IV: HrTS only"),
                     values = c("darkorange", "#00BA38", "#619CFF", "grey40"))  +
  theme_line()+
  guides(fill = "none") +
  theme(axis.text.y = element_text( margin = margin(r = 12)),
        axis.text.x = element_text( margin = margin(t = 12)),
        axis.title.y = element_text(vjust = + 3),
        axis.title.x = element_text(vjust = -0.75),
        plot.margin = margin(16, 16, 16, 16),
        legend.position="bottom")

# theme(
#   legend.position = c(.95, .95),
#   legend.justification = c("right", "top"),
#   legend.box.just = "right",
#   legend.margin = margin(6, 6, 6, 6)
# ) # 1240 x 500 # 580x500

inset = ggplot() +
  geom_line(aes(x = Year + 0.5, y = value, col = as.factor(strategy)), data = preds_sq, size = 0.75) +
  geom_ribbon(aes(x = Year + 0.5, ymin = CI_lower, ymax = CI_upper, fill = as.factor(strategy)), alpha = 0.3, data = preds_sq) +
  geom_point(aes(x = Year + 0.5, y = mean, col = as.factor(Strategy)), data = time2TB_national.summary, size = 2, shape = 1) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 100)) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0 , 10, 1)) +
  xlab("Years since entry") +
  ylab("Number of cases") +
  # ggtitle("Number of incident TB by year since entry") +
  scale_color_manual(name = "Strategy", labels = c("I: No screening", "II: IGRA-TB", "III: IGRA-RNA-TB", "IV: RNA-TB"),
                     values = c("darkorange", "#00BA38", "#619CFF", "grey40"))  +
  theme_line() +
  guides(fill = "none") +
  guides(col = "none") #+ theme(
#   legend.position = c(.95, .95),
#   legend.justification = c("right", "top"),
#   legend.box.just = "right",
#   legend.margin = margin(6, 6, 6, 6)
# ) # 580x500

#& --- Figure2 --- &#
library(cowplot)
ggdraw() + draw_plot(main)+ draw_plot(inset, x = 0.45, y = 0.5, width = 0.5, height = 0.45) #935x540 main-res-fig1
# ----------------- #

# Other smoothing methods - not used in final write-up

# tidymv::plot_smooths(
#   model = gam_time2TB_s0,
#   series = Year,
#   series_length = 90,
#   transform = exp) +
#   theme(legend.position = "top")
#
# plot(time2TB_nationa.summary[time2TB_nationa.summary$Strategy == 0,]$mean, col = "blue")
# lines(gam_time2TB_s0$fitted.values)
#
# lines(gam_time2TB_s0_sq$fitted.values, col = "red")
#
# plot(gam_time2TB_s0, residuals = T, pch = 16, shift = time2TB_national.summary[time2TB_national.summary$Strategy == 0,]$mean, ylab = "TB cases")

# -------------------------------------------------------------------------------------------------
# REDUCTION IN TB CASES OVER LIFETIME
# -------------------------------------------------------------------------------------------------

# tb_tab_array = array(NA, dim = c(97,1000, 4))
# tab.tb = lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = "TB"))
# tab.tb = do.call("cbind", tab.tb)
# tb_tab_array[,,1] = tab.tb[1:97,]
# tb_tab_array[,,2] = tab.tb[(97*1 +1): (97*2),]
# tb_tab_array[,,3] = tab.tb[(97*2 +1): (97*3),]
# tb_tab_array[,,4] = tab.tb[(97*3 +1): (97*4),]
# saveRDS(tb_tab_array, file = paste0(CURRENT_RESDIR,"/tb_tab_array.RDS"))

tb_tab_array = readRDS(paste0(CURRENT_RESDIR,"/tb_tab_array.RDS"))
tb_tab = apply(tb_tab_array, c(3), rowMeans)
tb_tab = as.data.frame(tb_tab)
tb_tab$country = countrynames

tb_tab_lower = apply(tb_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
tb_tab_lower = as.data.frame(tb_tab_lower)
tb_tab_lower$country = countrynames

tb_tab_upper = apply(tb_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
tb_tab_upper = as.data.frame(tb_tab_upper)
tb_tab_upper$country = countrynames

#& --- Table3 --- &#
# Use print( , digits = 6) to control the digits

# total
total = colSums(tb_tab[,1:4])
total_lower = colSums(tb_tab_lower[,1:4])
total_upper = colSums(tb_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(tb_tab)
cat_abs_lower = get_abs_by_riskcat(tb_tab_lower)
cat_abs_upper = get_abs_by_riskcat(tb_tab_upper)

# relative reduction
total_rel = 1- total/total[1]
total_rel_lower = 1 - total_upper/total_upper[1]
total_rel_upper = 1 - total_lower/total_lower[1]
cat_rel = 1 - cat_abs/cat_abs[,1]
cat_rel_lower = 1 - cat_abs_upper/cat_abs_upper[,1] # this is correct
cat_rel_upper = 1 - cat_abs_lower/cat_abs_lower[,1]
# ---------------- #

total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
cat(apply(total.df, 2, function(x) pretty(x*100, 1, ", ")), sep = "\n")

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)
cat(apply(cat.df[c(1, 1+4, 1+4*2),], 2, function(x) pretty(x*100, 1, ", ")), sep = "\n")
cat(apply(cat.df[c(2, 2+4, 2+4*2),], 2, function(x) pretty(x*100, 1, ", ")), sep = "\n")
cat(apply(cat.df[c(3, 3+4, 3+4*2),], 2, function(x) pretty(x*100, 1, ", ")), sep = "\n")
cat(apply(cat.df[c(4, 4+4, 4+4*2),], 2, function(x) pretty(x*100, 1, ", ")), sep = "\n")


# =======================================================================================
# TESTING RESOURCE UTILIZATION
# =======================================================================================

# TEST.IGRA
igra_tab_array = array(NA, dim = c(97, 1000, 4))
tab.igra = future_lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = "TEST.IGRA"))
tab.igra = do.call("cbind", tab.igra)
igra_tab_array[,,1] = tab.igra[1:97,]
igra_tab_array[,,2] = tab.igra[(97*1 +1): (97*2),]
igra_tab_array[,,3] = tab.igra[(97*2 +1): (97*3),]
igra_tab_array[,,4] = tab.igra[(97*3 +1): (97*4),]

saveRDS(igra_tab_array, file = paste0(CURRENT_RESDIR,"/igra_tab_array.RDS"))

print("Done writing: igra_tab_array.RDS")

# TEST.RNA
rna_tab_array = array(NA, dim = c(97, 1000, 4))
tab.rna = future_lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = "TEST.RNA"))
tab.rna = do.call("cbind", tab.rna)
rna_tab_array[,,1] = tab.rna[1:97,]
rna_tab_array[,,2] = tab.rna[(97*1 +1): (97*2),]
rna_tab_array[,,3] = tab.rna[(97*2 +1): (97*3),]
rna_tab_array[,,4] = tab.rna[(97*3 +1): (97*4),]

saveRDS(rna_tab_array, file = paste0(CURRENT_RESDIR,"/rna_tab_array.RDS"))

print("Done writing: rna_tab_array.RDS")


# TEST.XRAY
xray_tab_array = array(NA, dim = c(97, 1000, 4))
tab.xray = future_lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = "TEST.XRAY"))
tab.xray = do.call("cbind", tab.xray)
xray_tab_array[,,1] = tab.xray[1:97,]
xray_tab_array[,,2] = tab.xray[(97*1 +1): (97*2),]
xray_tab_array[,,3] = tab.xray[(97*2 +1): (97*3),]
xray_tab_array[,,4] = tab.xray[(97*3 +1): (97*4),]

saveRDS(xray_tab_array, file = paste0(CURRENT_RESDIR,"/xray_tab_array.RDS"))

print("Done writing: xray_tab_array.RDS")

# IGRA -------------------------------------------------------------------------------
igra_tab_array = readRDS(paste0(CURRENT_RESDIR,"/igra_tab_array.RDS"))
igra_tab = apply(igra_tab_array, c(3), rowMeans)
igra_tab = as.data.frame(igra_tab)
igra_tab$country = countrynames

igra_tab_lower = apply(igra_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
igra_tab_lower = as.data.frame(igra_tab_lower)
igra_tab_lower$country = countrynames

igra_tab_upper = apply(igra_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
igra_tab_upper = as.data.frame(igra_tab_upper)
igra_tab_upper$country = countrynames

#& --- Table4 --- &#
# total
total = colSums(igra_tab[,1:4])
total_lower = colSums(igra_tab_lower[,1:4])
total_upper = colSums(igra_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(igra_tab)
cat_abs_lower = get_abs_by_riskcat(igra_tab_lower)
cat_abs_upper = get_abs_by_riskcat(igra_tab_upper)

# relative
total_rel = total/sum(cohort2019_by_coo2$total_pop)
total_rel_lower = total_upper/sum(cohort2019_by_coo2$total_pop)
total_rel_upper = total_lower/sum(cohort2019_by_coo2$total_pop)
cat_rel = cat_abs/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_lower = cat_abs_lower/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_upper = cat_abs_upper/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
# ----------------- #

total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
apply(total.df, 2, function(x) pretty(x*100, 1, ", "))

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)

cat(apply(cat.df[c(1, 1+4, 1+4*2),] *100, 2, function(x) pretty(x, 1, ", ")), sep = "\n")
cat(apply(cat.df[c(2, 2+4, 2+4*2),] *100, 2, function(x) pretty(x, 1, ", ")), sep = "\n")
cat(apply(cat.df[c(3, 3+4, 3+4*2),] *100, 2, function(x) pretty(x, 1, ", ")), sep = "\n")
cat(apply(cat.df[c(4, 4+4, 4+4*2),] *100, 2, function(x) pretty(x, 1, ", ")), sep = "\n")
# RNA -------------------------------------------------------------------------------------------
rna_tab_array = readRDS(paste0(CURRENT_RESDIR,"/rna_tab_array.RDS"))
rna_tab = apply(rna_tab_array, c(3), rowMeans)
rna_tab = as.data.frame(rna_tab)
rna_tab$country = countrynames

rna_tab_lower = apply(rna_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
rna_tab_lower = as.data.frame(rna_tab_lower)
rna_tab_lower$country = countrynames

rna_tab_upper = apply(rna_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
rna_tab_upper = as.data.frame(rna_tab_upper)
rna_tab_upper$country = countrynames

#& --- Table4 --- &#
# total
total = colSums(rna_tab[,1:4])
total_lower = colSums(rna_tab_lower[,1:4])
total_upper = colSums(rna_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(rna_tab)
cat_abs_lower = get_abs_by_riskcat(rna_tab_lower)
cat_abs_upper = get_abs_by_riskcat(rna_tab_upper)

# relative
total_rel = total/sum(cohort2019_by_coo2$total_pop)
total_rel_lower = total_lower/sum(cohort2019_by_coo2$total_pop)
total_rel_upper = total_upper/sum(cohort2019_by_coo2$total_pop)
cat_rel = cat_abs/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_lower = cat_abs_lower/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_upper = cat_abs_upper/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
# ----------------- #


total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
apply(total.df, 2, function(x) pretty(x*100, 1, paste0(" ", expression( "\U2012" )," ")))

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)

cat(apply(cat.df[c(1, 1+4, 1+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(2, 2+4, 2+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(3, 3+4, 3+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(4, 4+4, 4+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")


# XRAY -------------------------------------------------------------------------------------------

xray_tab_array = readRDS(paste0(CURRENT_RESDIR,"/xray_tab_array.RDS"))
xray_tab = apply(xray_tab_array, c(3), rowMeans)
xray_tab = as.data.frame(xray_tab)
xray_tab$country = countrynames

xray_tab_lower = apply(xray_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
xray_tab_lower = as.data.frame(xray_tab_lower)
xray_tab_lower$country = countrynames

xray_tab_upper = apply(xray_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
xray_tab_upper = as.data.frame(xray_tab_upper)
xray_tab_upper$country = countrynames

#& --- Table4 --- &#
# total
total = colSums(xray_tab[,1:4])
total_lower = colSums(xray_tab_lower[,1:4])
total_upper = colSums(xray_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(xray_tab)
cat_abs_lower = get_abs_by_riskcat(xray_tab_lower)
cat_abs_upper = get_abs_by_riskcat(xray_tab_upper)

# relative
total_rel = total/sum(cohort2019_by_coo2$total_pop)
total_rel_lower = total_lower/sum(cohort2019_by_coo2$total_pop)
total_rel_upper = total_upper/sum(cohort2019_by_coo2$total_pop)
cat_rel = cat_abs/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_lower = cat_abs_lower/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_upper = cat_abs_upper/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
# ----------------- #


total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
apply(total.df, 2, function(x) pretty(x*100, 2, paste0(" ", expression( "\U2012" )," ")))

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)

cat(apply(cat.df[c(1, 1+4, 1+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(2, 2+4, 2+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(3, 3+4, 3+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(4, 4+4, 4+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")



# ==============================================================================
# TREATMENT RESOURCE UTILIZATION
# ==============================================================================


# LTBI TREATMENT
ltbiTx_array = array(NA, dim = c(97, 1000, 4))
tab.ltbiTx = future_lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = "TREAT.LTBI_INIT.SUCCESS"))
tab.ltbiTx = do.call("cbind", tab.ltbiTx)
ltbiTx_array[,,1] = tab.ltbiTx[1:97,]
ltbiTx_array[,,2] = tab.ltbiTx[(97*1 +1): (97*2),]
ltbiTx_array[,,3] = tab.ltbiTx[(97*2 +1): (97*3),]
ltbiTx_array[,,4] = tab.ltbiTx[(97*3 +1): (97*4),]

saveRDS(ltbiTx_array, file = paste0(CURRENT_RESDIR,"/ltbiTx_array.RDS"))
rm(ltbiTx_array); rm(tab.ltbiTx)

print("Done writing: ltbiTx_array.RDS")


# INCIPIENT TB TREATMENT
incipTx_array = array(NA, dim = c(97, 1000, 4))
tab.incipTx = future_lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = "TREAT.INCIPIENT_INIT.SUCCESS"))
tab.incipTx = do.call("cbind", tab.incipTx)
incipTx_array[,,1] = tab.incipTx[1:97,]
incipTx_array[,,2] = tab.incipTx[(97*1 +1): (97*2),]
incipTx_array[,,3] = tab.incipTx[(97*2 +1): (97*3),]
incipTx_array[,,4] = tab.incipTx[(97*3 +1): (97*4),]

saveRDS(incipTx_array, file = paste0(CURRENT_RESDIR,"/incipTx_array.RDS"))
rm(incipTx_array); rm(tab.incipTx)

print("Done writing: incipTx_array.RDS")


# TB TREATMENT
tbTx_array = array(NA, dim = c(97, 1000, 4))
tab.tbTx = future_lapply(1:1000, function(x) extract_n_per_paramset(paramsetnum = x, event = c("TREAT.TB.GUARANTEE", "TREAT.TB_INIT.SUCCESS")))
tab.tbTx = do.call("cbind", tab.tbTx)
tbTx_array[,,1] = tab.tbTx[1:97,]
tbTx_array[,,2] = tab.tbTx[(97*1 +1): (97*2),]
tbTx_array[,,3] = tab.tbTx[(97*2 +1): (97*3),]
tbTx_array[,,4] = tab.tbTx[(97*3 +1): (97*4),]

saveRDS(tbTx_array, file = paste0(CURRENT_RESDIR,"/tbTx_array.RDS"))
rm(tbTx_array); rm(tab.tbTx)

print("Done writing: tbTx_array.RDS")


# LTBI TREATMENT ----------------------------------------------------------------

ltbiTx_tab_array = readRDS(paste0(CURRENT_RESDIR,"/ltbiTx_array.RDS"))
ltbiTx_tab = apply(ltbiTx_tab_array, c(3), rowMeans)
ltbiTx_tab = as.data.frame(ltbiTx_tab)
ltbiTx_tab$country = countrynames

ltbiTx_tab_lower = apply(ltbiTx_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
ltbiTx_tab_lower = as.data.frame(ltbiTx_tab_lower)
ltbiTx_tab_lower$country = countrynames

ltbiTx_tab_upper = apply(ltbiTx_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
ltbiTx_tab_upper = as.data.frame(ltbiTx_tab_upper)
ltbiTx_tab_upper$country = countrynames

#& --- Table5 --- &#
# total
total = colSums(ltbiTx_tab[,1:4])
total_lower = colSums(ltbiTx_tab_lower[,1:4])
total_upper = colSums(ltbiTx_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(ltbiTx_tab)
cat_abs_lower = get_abs_by_riskcat(ltbiTx_tab_lower)
cat_abs_upper = get_abs_by_riskcat(ltbiTx_tab_upper)

# relative
total_rel = total/sum(cohort2019_by_coo2$total_pop)
total_rel_lower = total_lower/sum(cohort2019_by_coo2$total_pop)
total_rel_upper = total_upper/sum(cohort2019_by_coo2$total_pop)
cat_rel = cat_abs/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_lower = cat_abs_lower/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_upper = cat_abs_upper/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
# ----------------- #


total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
cat(apply(total.df, 2, function(x) pretty(x*100, 1, paste0(" ", expression( "\U2012" )," "))), sep = "\n")

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)

cat(apply(cat.df[c(1, 1+4, 1+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = "\n")
cat(apply(cat.df[c(2, 2+4, 2+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = "\n")
cat(apply(cat.df[c(3, 3+4, 3+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = "\n")
cat(apply(cat.df[c(4, 4+4, 4+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = "\n")



# INCIPIENT TB TREATMENT --------------------------------------------------------


incipTx_tab_array = readRDS(paste0(CURRENT_RESDIR,"/incipTx_array.RDS"))
incipTx_tab = apply(incipTx_tab_array, c(3), rowMeans)
incipTx_tab = as.data.frame(incipTx_tab)
incipTx_tab$country = countrynames

incipTx_tab_lower = apply(incipTx_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
incipTx_tab_lower = as.data.frame(incipTx_tab_lower)
incipTx_tab_lower$country = countrynames

incipTx_tab_upper = apply(incipTx_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
incipTx_tab_upper = as.data.frame(incipTx_tab_upper)
incipTx_tab_upper$country = countrynames

#& --- Table5 --- &#
# total
total = colSums(incipTx_tab[,1:4])
total_lower = colSums(incipTx_tab_lower[,1:4])
total_upper = colSums(incipTx_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(incipTx_tab)
cat_abs_lower = get_abs_by_riskcat(incipTx_tab_lower)
cat_abs_upper = get_abs_by_riskcat(incipTx_tab_upper)

# relative
total_rel = total/sum(cohort2019_by_coo2$total_pop)
total_rel_lower = total_lower/sum(cohort2019_by_coo2$total_pop)
total_rel_upper = total_upper/sum(cohort2019_by_coo2$total_pop)
cat_rel = cat_abs/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_lower = cat_abs_lower/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_upper = cat_abs_upper/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
# ----------------- #


total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
cat(apply(total.df, 2, function(x) pretty(x*100, 1, paste0(" ", expression( "\U2012" )," "))))

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)

cat(apply(cat.df[c(1, 1+4, 1+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(2, 2+4, 2+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(3, 3+4, 3+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")
cat(apply(cat.df[c(4, 4+4, 4+4*2),] *100, 2, function(x) pretty(x, 1, paste0(" ", expression( "\U2012" )," "))), sep = " ")


# ACTIVE TB TREATMENT -----------------------------------------------------------
tbTx_tab_array = readRDS(paste0(CURRENT_RESDIR,"/tbTx_array.RDS"))
tbTx_tab = apply(tbTx_tab_array, c(3), rowMeans)
tbTx_tab = as.data.frame(tbTx_tab)
tbTx_tab$country = countrynames

tbTx_tab_lower = apply(tbTx_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.025))
tbTx_tab_lower = as.data.frame(tbTx_tab_lower)
tbTx_tab_lower$country = countrynames

tbTx_tab_upper = apply(tbTx_tab_array, c(3), function(x) matrixStats::rowQuantiles(x, p = 0.975))
tbTx_tab_upper = as.data.frame(tbTx_tab_upper)
tbTx_tab_upper$country = countrynames

#& --- Table5 --- &#
# total
total = colSums(tbTx_tab[,1:4])
total_lower = colSums(tbTx_tab_lower[,1:4])
total_upper = colSums(tbTx_tab_upper[,1:4])

# by riskcat
cat_abs = get_abs_by_riskcat(tbTx_tab)
cat_abs_lower = get_abs_by_riskcat(tbTx_tab_lower)
cat_abs_upper = get_abs_by_riskcat(tbTx_tab_upper)

# relative
total_rel = total/sum(cohort2019_by_coo2$total_pop)
total_rel_lower = total_lower/sum(cohort2019_by_coo2$total_pop)
total_rel_upper = total_upper/sum(cohort2019_by_coo2$total_pop)
cat_rel = cat_abs/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_lower = cat_abs_lower/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
cat_rel_upper = cat_abs_upper/as.numeric(coo_by_riskcat[order(coo_by_riskcat$tb_inc, decreasing = T),]$pop_size)
# ----------------- #


total.df = rbind(total_rel, total_rel_lower, total_rel_upper)
cat(apply(total.df, 2, function(x) pretty(x*100, 2, paste0(" ", expression( "\U2012" )," "))), sep = "\n")

cat.df = rbind(cat_rel, cat_rel_lower, cat_rel_upper)

cat(apply(cat.df[c(1, 1+4, 1+4*2),] *100, 2, function(x) pretty(x, 2, paste0(" ", expression( "\U2012" )," "))), sep = "\n")
cat(apply(cat.df[c(2, 2+4, 2+4*2),] *100, 2, function(x) pretty(x, 2, paste0(" ", expression( "\U2012" )," "))), sep = "\n")
cat(apply(cat.df[c(3, 3+4, 3+4*2),] *100, 2, function(x) pretty(x, 2, paste0(" ", expression( "\U2012" )," "))), sep = "\n")
cat(apply(cat.df[c(4, 4+4, 4+4*2),] *100, 2, function(x) pretty(x, 2, paste0(" ", expression( "\U2012" )," "))), sep = "\n")


# =================================================================================================
# RESULTS: PPV, NPV
# =================================================================================================

# Check whether all parameter sets were simulated
# dir_list = list.dirs(path = paste0(CURRENT_RESDIR),recursive = FALSE)
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,71, nchar(dir_list)))) #PSA_main
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,74, nchar(dir_list)))) #PSA_noreinf
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, nchar(dir_list)))) #PSA_reinf2
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, nchar(dir_list)))) #PSA_RNAgam
# 
# if (length(ptemp) > 0) stop("Some parameter sets were not simulated.")
# 
# # Check whether all countries were simulated for each data set
# get_file_list2 = function(dir){
#   file_list = list.files(path = paste0(dir), pattern = paste0("ppvnpv.csv"))
#   return(length(file_list))
# }
# 
# ctemp = future_sapply(dir_list, function(x) get_file_list2(x), future.seed = NULL)
# etemp = str_sub(dir_list[which(ctemp != 1)], 71, nchar(dir_list[which(ctemp != 1)]))
# etemp = str_sub(dir_list[which(ctemp != 1)], 73, nchar(dir_list[which(ctemp != 1)]))
# etemp = str_sub(dir_list[which(ctemp != 1)], 74, nchar(dir_list[which(ctemp != 1)]))

# Useful functions
get_ppvnpv_paramset = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv.csv'))
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  
  dat = dat %>%
    # rename(., testres.s2 = testres.s2_mod) %>%
    # rename(., testres.s3 = testres.s3_mod) %>%
    # rename(., testres.s4 = testres.s4_mod) %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4))
  
  # ppv.s2 = sum(dat$lifetimeTB == 1 & dat$testres.s2 == 1)/sum(dat$testres.s2 == 1)
  # ppv.s3 = sum(dat$lifetimeTB == 1 & dat$testres.s3 == 1)/sum(dat$testres.s3 == 1)
  # ppv.s4 = sum(dat$lifetimeTB == 1 & dat$testres.s4 == 1)/sum(dat$testres.s4 == 1)
  
  ppvLT.s2 = mean(dat[dat$testres.s2 == 1,]$lifetimeTB)
  ppvLT.s3 = mean(dat[dat$testres.s3 == 1,]$lifetimeTB)
  ppvLT.s4 = mean(dat[dat$testres.s4 == 1,]$lifetimeTB)
  
  # npv.s2 = sum(dat$lifetimeTB == 0 & dat$testres.s2 == 0)/sum(dat$testres.s2 == 0)
  # npv.s3 = sum(dat$lifetimeTB == 0 & dat$testres.s3 == 0)/sum(dat$testres.s3 == 0)
  # npv.s4 = sum(dat$lifetimeTB == 0 & dat$testres.s4 == 0)/sum(dat$testres.s4 == 0)
  
  npvLT.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$lifetimeTB)
  npvLT.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$lifetimeTB)
  npvLT.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$lifetimeTB)
  
  
  ppvY2.s2 = mean(dat[dat$testres.s2 == 1,]$TBin2years)
  ppvY2.s3 = mean(dat[dat$testres.s3 == 1,]$TBin2years)
  ppvY2.s4 = mean(dat[dat$testres.s4 == 1,]$TBin2years)
  
  npvY2.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$TBin2years)
  npvY2.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$TBin2years)
  npvY2.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$TBin2years)
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4)
  
  return(out)
}


get_ppvnpv_paramset_newinfec = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR, "/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv_rev.csv'))
  # rawdat2 = read.csv(paste0(paramset_dir,'ppvnpv_tb2.csv'))
  # 
  # rawdat = cbind(rawdat1, rawdat2[,c("lifetimeTB_notReInf", "TBin2years_notReInf")])
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  dat = dat %>%
    rename(., testres.s2 = testres.s2_mod) %>%
    rename(., testres.s3 = testres.s3_mod) %>%
    rename(., testres.s4 = testres.s4_mod) %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4))
  
  dat_s2_1 = dat[dat$testres.s2 == 1,]
  dat_s3_1 = dat[dat$testres.s3 == 1,]
  dat_s4_1 = dat[dat$testres.s4 == 1,]
  
  dat_s2_0 = dat[dat$testres.s2 == 0,]
  dat_s3_0 = dat[dat$testres.s3 == 0,]
  dat_s4_0 = dat[dat$testres.s4 == 0,]
  
  ppvLT.s2 = mean(dat_s2_1$lifetimeTB)
  ppvLT.s3 = mean(dat_s3_1$lifetimeTB)
  ppvLT.s4 = mean(dat_s4_1$lifetimeTB)
  
  npvLT.s2 = 1 - mean(dat_s2_0$lifetimeTB)
  npvLT.s3 = 1 - mean(dat_s3_0$lifetimeTB)
  npvLT.s4 = 1 - mean(dat_s4_0$lifetimeTB)
  
  ppvY2.s2 = mean(dat_s2_1$TBin2years)
  ppvY2.s3 = mean(dat_s3_1$TBin2years)
  ppvY2.s4 = mean(dat_s4_1$TBin2years)
  
  npvY2.s2 = 1 - mean(dat_s2_0$TBin2years)
  npvY2.s3 = 1 - mean(dat_s3_0$TBin2years)
  npvY2.s4 = 1 - mean(dat_s4_0$TBin2years)
  
  ppvLT.s2.nri = mean(dat_s2_1$lifetimeTB_notReInf)
  ppvLT.s3.nri = mean(dat_s3_1$lifetimeTB_notReInf)
  ppvLT.s4.nri = mean(dat_s4_1$lifetimeTB_notReInf)
  
  npvLT.s2.nri = 1 - mean(dat_s2_0$lifetimeTB_notReInf)
  npvLT.s3.nri = 1 - mean(dat_s3_0$lifetimeTB_notReInf)
  npvLT.s4.nri = 1 - mean(dat_s4_0$lifetimeTB_notReInf)
  
  ppvY2.s2.nri = mean(dat_s2_1$TBin2years_notReInf)
  ppvY2.s3.nri = mean(dat_s3_1$TBin2years_notReInf)
  ppvY2.s4.nri = mean(dat_s4_1$TBin2years_notReInf)
  
  npvY2.s2.nri = 1 - mean(dat_s2_0$TBin2years_notReInf)
  npvY2.s3.nri = 1 - mean(dat_s3_0$TBin2years_notReInf)
  npvY2.s4.nri = 1 - mean(dat_s4_0$TBin2years_notReInf)
  
  
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4,
          ppvLT.s2.nri, ppvLT.s3.nri, ppvLT.s4.nri, npvLT.s2.nri, npvLT.s3.nri, npvLT.s4.nri,
          ppvY2.s2.nri, ppvY2.s3.nri, ppvY2.s4.nri, npvY2.s2.nri, npvY2.s3.nri, npvY2.s4.nri)
  
  return(out)
}

get_riskcat_ppvnpv = function(dat, countries_in_riskcat){
  
  dat = dat %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4)) %>%
    filter(CountryName %in% countries_in_riskcat)
  
  # ppv.s2 = sum(dat$lifetimeTB == 1 & dat$testres.s2 == 1)/sum(dat$testres.s2 == 1)
  # ppv.s3 = sum(dat$lifetimeTB == 1 & dat$testres.s3 == 1)/sum(dat$testres.s3 == 1)
  # ppv.s4 = sum(dat$lifetimeTB == 1 & dat$testres.s4 == 1)/sum(dat$testres.s4 == 1)
  
  ppvLT.s2 = mean(dat[dat$testres.s2 == 1,]$lifetimeTB)
  ppvLT.s3 = mean(dat[dat$testres.s3 == 1,]$lifetimeTB)
  ppvLT.s4 = mean(dat[dat$testres.s4 == 1,]$lifetimeTB)
  
  # npv.s2 = sum(dat$lifetimeTB == 0 & dat$testres.s2 == 0)/sum(dat$testres.s2 == 0)
  # npv.s3 = sum(dat$lifetimeTB == 0 & dat$testres.s3 == 0)/sum(dat$testres.s3 == 0)
  # npv.s4 = sum(dat$lifetimeTB == 0 & dat$testres.s4 == 0)/sum(dat$testres.s4 == 0)
  
  npvLT.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$lifetimeTB)
  npvLT.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$lifetimeTB)
  npvLT.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$lifetimeTB)
  
  ppvY2.s2 = mean(dat[dat$testres.s2 == 1,]$TBin2years)
  ppvY2.s3 = mean(dat[dat$testres.s3 == 1,]$TBin2years)
  ppvY2.s4 = mean(dat[dat$testres.s4 == 1,]$TBin2years)
  
  npvY2.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$TBin2years)
  npvY2.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$TBin2years)
  npvY2.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$TBin2years)
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4)
  
  return(out)
}

get_riskcat_ppvnpv_newinfec = function(dat, countries_in_riskcat){
  
  dat = dat %>%
    rename(., testres.s2 = testres.s2_mod) %>%
    rename(., testres.s3 = testres.s3_mod) %>%
    rename(., testres.s4 = testres.s4_mod) %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4)) %>%
    filter(CountryName %in% countries_in_riskcat)
  
  dat_s2_1 = dat[dat$testres.s2 == 1,]
  dat_s3_1 = dat[dat$testres.s3 == 1,]
  dat_s4_1 = dat[dat$testres.s4 == 1,]
  
  dat_s2_0 = dat[dat$testres.s2 == 0,]
  dat_s3_0 = dat[dat$testres.s3 == 0,]
  dat_s4_0 = dat[dat$testres.s4 == 0,]
  
  ppvLT.s2 = mean(dat_s2_1$lifetimeTB)
  ppvLT.s3 = mean(dat_s3_1$lifetimeTB)
  ppvLT.s4 = mean(dat_s4_1$lifetimeTB)
  
  npvLT.s2 = 1 - mean(dat_s2_0$lifetimeTB)
  npvLT.s3 = 1 - mean(dat_s3_0$lifetimeTB)
  npvLT.s4 = 1 - mean(dat_s4_0$lifetimeTB)
  
  ppvY2.s2 = mean(dat_s2_1$TBin2years)
  ppvY2.s3 = mean(dat_s3_1$TBin2years)
  ppvY2.s4 = mean(dat_s4_1$TBin2years)
  
  npvY2.s2 = 1 - mean(dat_s2_0$TBin2years)
  npvY2.s3 = 1 - mean(dat_s3_0$TBin2years)
  npvY2.s4 = 1 - mean(dat_s4_0$TBin2years)
  
  ppvLT.s2.nri = mean(dat_s2_1$lifetimeTB_notReInf)
  ppvLT.s3.nri = mean(dat_s3_1$lifetimeTB_notReInf)
  ppvLT.s4.nri = mean(dat_s4_1$lifetimeTB_notReInf)
  
  npvLT.s2.nri = 1 - mean(dat_s2_0$lifetimeTB_notReInf)
  npvLT.s3.nri = 1 - mean(dat_s3_0$lifetimeTB_notReInf)
  npvLT.s4.nri = 1 - mean(dat_s4_0$lifetimeTB_notReInf)
  
  ppvY2.s2.nri = mean(dat_s2_1$TBin2years_notReInf)
  ppvY2.s3.nri = mean(dat_s3_1$TBin2years_notReInf)
  ppvY2.s4.nri = mean(dat_s4_1$TBin2years_notReInf)
  
  npvY2.s2.nri = 1 - mean(dat_s2_0$TBin2years_notReInf)
  npvY2.s3.nri = 1 - mean(dat_s3_0$TBin2years_notReInf)
  npvY2.s4.nri = 1 - mean(dat_s4_0$TBin2years_notReInf)
  
  
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4,
          ppvLT.s2.nri, ppvLT.s3.nri, ppvLT.s4.nri, npvLT.s2.nri, npvLT.s3.nri, npvLT.s4.nri,
          ppvY2.s2.nri, ppvY2.s3.nri, ppvY2.s4.nri, npvY2.s2.nri, npvY2.s3.nri, npvY2.s4.nri)
  
  return(out)
}


get_riskcat_ppvnpv_paramset = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv.csv'))
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  # out = future_lapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv(dat, x))
  out = mclapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv(dat, x), mc.cores = numCores)
  
  out = do.call("rbind", out)
  out = cbind(rep(paramsetnum, 4), c(1:4), out)
  colnames(out) = c("paramsetnum", "riskcat",
                    "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                    "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4")
  
  return(out)
}

get_riskcat_ppvnpv_paramset_newinfec = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv_rev.csv'))
  # rawdat2 = read.csv(paste0(paramset_dir,'ppvnpv_tb2.csv'))
  # 
  # rawdat = cbind(rawdat1, rawdat2[,c("lifetimeTB_notReInf", "TBin2years_notReInf")])
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  # out = future_lapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv(dat, x))
  out = mclapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv_newinfec(dat, x), mc.cores = numCores)
  
  out = do.call("rbind", out)
  out = cbind(rep(paramsetnum, 4), c(1:4), out)
  colnames(out) = c("paramsetnum", "riskcat",
                    "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                    "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4",
                    "ppvLT.s2.nri", "ppvLT.s3.nri", "ppvLT.s4.nri", "npvLT.s2.nri", "npvLT.s3.nri", "npvLT.s4.nri",
                    "ppvY2.s2.nri", "ppvY2.s3.nri", "ppvY2.s4.nri", "npvY2.s2.nri", "npvY2.s3.nri", "npvY2.s4.nri")
  
  return(out)
}




if (SA == "PSA_noreinf"){
  # # ppvnpv_allparamsets = future_lapply(1:1000, function(x) get_ppvnpv_paramset(x)) # this takes a while
  ppvnpv_allparamsets = mclapply(1:1000, function(x) get_ppvnpv_paramset(x), mc.cores = numCores) # this takes a while
  ppvnpv_allparamsets = do.call('rbind', ppvnpv_allparamsets)
  colnames(ppvnpv_allparamsets)  = c("ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                                     "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4")
  write.csv(ppvnpv_allparamsets, paste0(CURRENT_RESDIR,"/ppvnpv_tot.csv"), row.names = F)
  
  
  # ppvnpv_riskcat = future_lapply(1:1000, function(x)get_riskcat_ppvnpv_paramset(x))
  ppvnpv_riskcat = mclapply(1:1000, function(x)get_riskcat_ppvnpv_paramset(x), mc.cores=numCores)
  ppvnpv_riskcat = do.call('rbind', ppvnpv_riskcat)
  colnames(ppvnpv_riskcat) = c("paramsetnum", "riskcat",
                               "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                               "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4")
  write.csv(ppvnpv_riskcat, paste0(CURRENT_RESDIR, "/ppvnpv_riskcat.csv"), row.names = F)
  
} else {
  ppvnpv_allparamsets = mclapply(1:1000, function(x) get_ppvnpv_paramset_newinfec(x), mc.cores = numCores) # this takes a while
  
  ppvnpv_allparamsets = do.call('rbind', ppvnpv_allparamsets)
  colnames(ppvnpv_allparamsets) = c("ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                                    "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4",
                                    "ppvLT.s2.nri", "ppvLT.s3.nri", "ppvLT.s4.nri", "npvLT.s2.nri", "npvLT.s3.nri", "npvLT.s4.nri",
                                    "ppvY2.s2.nri", "ppvY2.s3.nri", "ppvY2.s4.nri", "npvY2.s2.nri", "npvY2.s3.nri", "npvY2.s4.nri")
  write.csv(ppvnpv_allparamsets, paste0(CURRENT_RESDIR,"/ppvnpv_tot_rev.csv"), row.names = F)
  
  # ppvnpv_riskcat = mclapply(1:1000, function(x) get_riskcat_ppvnpv_paramset_newinfec(x), mc.cores=numCores)
  # 
  # ppvnpv_riskcat = do.call('rbind', ppvnpv_riskcat)
  # colnames(ppvnpv_riskcat) = c("paramsetnum", "riskcat",
  #                              "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
  #                              "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4",
  #                              "ppvLT.s2.nri", "ppvLT.s3.nri", "ppvLT.s4.nri", "npvLT.s2.nri", "npvLT.s3.nri", "npvLT.s4.nri",
  #                              "ppvY2.s2.nri", "ppvY2.s3.nri", "ppvY2.s4.nri", "npvY2.s2.nri", "npvY2.s3.nri", "npvY2.s4.nri")
  # write.csv(ppvnpv_riskcat, paste0(CURRENT_RESDIR, "/ppvnpv_riskcat_rev.csv"), row.names = F)
}

# =================================================================================================
# RESULTS: PPV, NPV
# =================================================================================================

# Check whether all parameter sets were simulated
# dir_list = list.dirs(path = paste0(CURRENT_RESDIR),recursive = FALSE)
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,71, nchar(dir_list)))) #PSA_main
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,74, nchar(dir_list)))) #PSA_noreinf
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, nchar(dir_list)))) #PSA_reinf2
# ptemp = which(!c(1:1000) %in% c(str_sub(dir_list,73, nchar(dir_list)))) #PSA_RNAgam
# 
# if (length(ptemp) > 0) stop("Some parameter sets were not simulated.")
# 
# # Check whether all countries were simulated for each data set
# get_file_list2 = function(dir){
#   file_list = list.files(path = paste0(dir), pattern = paste0("ppvnpv.csv"))
#   return(length(file_list))
# }
# 
# ctemp = future_sapply(dir_list, function(x) get_file_list2(x), future.seed = NULL)
# etemp = str_sub(dir_list[which(ctemp != 1)], 71, nchar(dir_list[which(ctemp != 1)]))
# etemp = str_sub(dir_list[which(ctemp != 1)], 73, nchar(dir_list[which(ctemp != 1)]))
# etemp = str_sub(dir_list[which(ctemp != 1)], 74, nchar(dir_list[which(ctemp != 1)]))

# Useful functions
get_ppvnpv_paramset = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv.csv'))
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  
  dat = dat %>%
    # rename(., testres.s2 = testres.s2_mod) %>%
    # rename(., testres.s3 = testres.s3_mod) %>%
    # rename(., testres.s4 = testres.s4_mod) %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4))
  
  # ppv.s2 = sum(dat$lifetimeTB == 1 & dat$testres.s2 == 1)/sum(dat$testres.s2 == 1)
  # ppv.s3 = sum(dat$lifetimeTB == 1 & dat$testres.s3 == 1)/sum(dat$testres.s3 == 1)
  # ppv.s4 = sum(dat$lifetimeTB == 1 & dat$testres.s4 == 1)/sum(dat$testres.s4 == 1)
  
  ppvLT.s2 = mean(dat[dat$testres.s2 == 1,]$lifetimeTB)
  ppvLT.s3 = mean(dat[dat$testres.s3 == 1,]$lifetimeTB)
  ppvLT.s4 = mean(dat[dat$testres.s4 == 1,]$lifetimeTB)
  
  # npv.s2 = sum(dat$lifetimeTB == 0 & dat$testres.s2 == 0)/sum(dat$testres.s2 == 0)
  # npv.s3 = sum(dat$lifetimeTB == 0 & dat$testres.s3 == 0)/sum(dat$testres.s3 == 0)
  # npv.s4 = sum(dat$lifetimeTB == 0 & dat$testres.s4 == 0)/sum(dat$testres.s4 == 0)
  
  npvLT.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$lifetimeTB)
  npvLT.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$lifetimeTB)
  npvLT.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$lifetimeTB)
  
  
  ppvY2.s2 = mean(dat[dat$testres.s2 == 1,]$TBin2years)
  ppvY2.s3 = mean(dat[dat$testres.s3 == 1,]$TBin2years)
  ppvY2.s4 = mean(dat[dat$testres.s4 == 1,]$TBin2years)
  
  npvY2.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$TBin2years)
  npvY2.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$TBin2years)
  npvY2.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$TBin2years)
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4)
  
  return(out)
}


get_ppvnpv_paramset_newinfec = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR, "/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv_rev.csv'))
  # rawdat2 = read.csv(paste0(paramset_dir,'ppvnpv_tb2.csv'))
  # 
  # rawdat = cbind(rawdat1, rawdat2[,c("lifetimeTB_notReInf", "TBin2years_notReInf")])
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  dat = dat %>%
    rename(., testres.s2 = testres.s2_mod) %>%
    rename(., testres.s3 = testres.s3_mod) %>%
    rename(., testres.s4 = testres.s4_mod) %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4))
  
  dat_s2_1 = dat[dat$testres.s2 == 1,]
  dat_s3_1 = dat[dat$testres.s3 == 1,]
  dat_s4_1 = dat[dat$testres.s4 == 1,]
  
  dat_s2_0 = dat[dat$testres.s2 == 0,]
  dat_s3_0 = dat[dat$testres.s3 == 0,]
  dat_s4_0 = dat[dat$testres.s4 == 0,]
  
  ppvLT.s2 = mean(dat_s2_1$lifetimeTB)
  ppvLT.s3 = mean(dat_s3_1$lifetimeTB)
  ppvLT.s4 = mean(dat_s4_1$lifetimeTB)
  
  npvLT.s2 = 1 - mean(dat_s2_0$lifetimeTB)
  npvLT.s3 = 1 - mean(dat_s3_0$lifetimeTB)
  npvLT.s4 = 1 - mean(dat_s4_0$lifetimeTB)
  
  ppvY2.s2 = mean(dat_s2_1$TBin2years)
  ppvY2.s3 = mean(dat_s3_1$TBin2years)
  ppvY2.s4 = mean(dat_s4_1$TBin2years)
  
  npvY2.s2 = 1 - mean(dat_s2_0$TBin2years)
  npvY2.s3 = 1 - mean(dat_s3_0$TBin2years)
  npvY2.s4 = 1 - mean(dat_s4_0$TBin2years)
  
  ppvLT.s2.nri = mean(dat_s2_1$lifetimeTB_notReInf)
  ppvLT.s3.nri = mean(dat_s3_1$lifetimeTB_notReInf)
  ppvLT.s4.nri = mean(dat_s4_1$lifetimeTB_notReInf)
  
  npvLT.s2.nri = 1 - mean(dat_s2_0$lifetimeTB_notReInf)
  npvLT.s3.nri = 1 - mean(dat_s3_0$lifetimeTB_notReInf)
  npvLT.s4.nri = 1 - mean(dat_s4_0$lifetimeTB_notReInf)
  
  ppvY2.s2.nri = mean(dat_s2_1$TBin2years_notReInf)
  ppvY2.s3.nri = mean(dat_s3_1$TBin2years_notReInf)
  ppvY2.s4.nri = mean(dat_s4_1$TBin2years_notReInf)
  
  npvY2.s2.nri = 1 - mean(dat_s2_0$TBin2years_notReInf)
  npvY2.s3.nri = 1 - mean(dat_s3_0$TBin2years_notReInf)
  npvY2.s4.nri = 1 - mean(dat_s4_0$TBin2years_notReInf)
  
  
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4,
          ppvLT.s2.nri, ppvLT.s3.nri, ppvLT.s4.nri, npvLT.s2.nri, npvLT.s3.nri, npvLT.s4.nri,
          ppvY2.s2.nri, ppvY2.s3.nri, ppvY2.s4.nri, npvY2.s2.nri, npvY2.s3.nri, npvY2.s4.nri)
  
  return(out)
}

get_riskcat_ppvnpv = function(dat, countries_in_riskcat){
  
  dat = dat %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4)) %>%
    filter(CountryName %in% countries_in_riskcat)
  
  # ppv.s2 = sum(dat$lifetimeTB == 1 & dat$testres.s2 == 1)/sum(dat$testres.s2 == 1)
  # ppv.s3 = sum(dat$lifetimeTB == 1 & dat$testres.s3 == 1)/sum(dat$testres.s3 == 1)
  # ppv.s4 = sum(dat$lifetimeTB == 1 & dat$testres.s4 == 1)/sum(dat$testres.s4 == 1)
  
  ppvLT.s2 = mean(dat[dat$testres.s2 == 1,]$lifetimeTB)
  ppvLT.s3 = mean(dat[dat$testres.s3 == 1,]$lifetimeTB)
  ppvLT.s4 = mean(dat[dat$testres.s4 == 1,]$lifetimeTB)
  
  # npv.s2 = sum(dat$lifetimeTB == 0 & dat$testres.s2 == 0)/sum(dat$testres.s2 == 0)
  # npv.s3 = sum(dat$lifetimeTB == 0 & dat$testres.s3 == 0)/sum(dat$testres.s3 == 0)
  # npv.s4 = sum(dat$lifetimeTB == 0 & dat$testres.s4 == 0)/sum(dat$testres.s4 == 0)
  
  npvLT.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$lifetimeTB)
  npvLT.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$lifetimeTB)
  npvLT.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$lifetimeTB)
  
  ppvY2.s2 = mean(dat[dat$testres.s2 == 1,]$TBin2years)
  ppvY2.s3 = mean(dat[dat$testres.s3 == 1,]$TBin2years)
  ppvY2.s4 = mean(dat[dat$testres.s4 == 1,]$TBin2years)
  
  npvY2.s2 = 1 - mean(dat[dat$testres.s2 == 0,]$TBin2years)
  npvY2.s3 = 1 - mean(dat[dat$testres.s3 == 0,]$TBin2years)
  npvY2.s4 = 1 - mean(dat[dat$testres.s4 == 0,]$TBin2years)
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4)
  
  return(out)
}

get_riskcat_ppvnpv_newinfec = function(dat, countries_in_riskcat){
  
  dat = dat %>%
    rename(., testres.s2 = testres.s2_mod) %>%
    rename(., testres.s3 = testres.s3_mod) %>%
    rename(., testres.s4 = testres.s4_mod) %>%
    filter(!is.na(testres.s2)) %>%
    filter(!is.na(testres.s3)) %>%
    filter(!is.na(testres.s4)) %>%
    filter(CountryName %in% countries_in_riskcat)
  
  dat_s2_1 = dat[dat$testres.s2 == 1,]
  dat_s3_1 = dat[dat$testres.s3 == 1,]
  dat_s4_1 = dat[dat$testres.s4 == 1,]
  
  dat_s2_0 = dat[dat$testres.s2 == 0,]
  dat_s3_0 = dat[dat$testres.s3 == 0,]
  dat_s4_0 = dat[dat$testres.s4 == 0,]
  
  ppvLT.s2 = mean(dat_s2_1$lifetimeTB)
  ppvLT.s3 = mean(dat_s3_1$lifetimeTB)
  ppvLT.s4 = mean(dat_s4_1$lifetimeTB)
  
  npvLT.s2 = 1 - mean(dat_s2_0$lifetimeTB)
  npvLT.s3 = 1 - mean(dat_s3_0$lifetimeTB)
  npvLT.s4 = 1 - mean(dat_s4_0$lifetimeTB)
  
  ppvY2.s2 = mean(dat_s2_1$TBin2years)
  ppvY2.s3 = mean(dat_s3_1$TBin2years)
  ppvY2.s4 = mean(dat_s4_1$TBin2years)
  
  npvY2.s2 = 1 - mean(dat_s2_0$TBin2years)
  npvY2.s3 = 1 - mean(dat_s3_0$TBin2years)
  npvY2.s4 = 1 - mean(dat_s4_0$TBin2years)
  
  ppvLT.s2.nri = mean(dat_s2_1$lifetimeTB_notReInf)
  ppvLT.s3.nri = mean(dat_s3_1$lifetimeTB_notReInf)
  ppvLT.s4.nri = mean(dat_s4_1$lifetimeTB_notReInf)
  
  npvLT.s2.nri = 1 - mean(dat_s2_0$lifetimeTB_notReInf)
  npvLT.s3.nri = 1 - mean(dat_s3_0$lifetimeTB_notReInf)
  npvLT.s4.nri = 1 - mean(dat_s4_0$lifetimeTB_notReInf)
  
  ppvY2.s2.nri = mean(dat_s2_1$TBin2years_notReInf)
  ppvY2.s3.nri = mean(dat_s3_1$TBin2years_notReInf)
  ppvY2.s4.nri = mean(dat_s4_1$TBin2years_notReInf)
  
  npvY2.s2.nri = 1 - mean(dat_s2_0$TBin2years_notReInf)
  npvY2.s3.nri = 1 - mean(dat_s3_0$TBin2years_notReInf)
  npvY2.s4.nri = 1 - mean(dat_s4_0$TBin2years_notReInf)
  
  
  
  out = c(ppvLT.s2, ppvLT.s3, ppvLT.s4, npvLT.s2, npvLT.s3, npvLT.s4,
          ppvY2.s2, ppvY2.s3, ppvY2.s4, npvY2.s2, npvY2.s3, npvY2.s4,
          ppvLT.s2.nri, ppvLT.s3.nri, ppvLT.s4.nri, npvLT.s2.nri, npvLT.s3.nri, npvLT.s4.nri,
          ppvY2.s2.nri, ppvY2.s3.nri, ppvY2.s4.nri, npvY2.s2.nri, npvY2.s3.nri, npvY2.s4.nri)
  
  return(out)
}


get_riskcat_ppvnpv_paramset = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv.csv'))
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  # out = future_lapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv(dat, x))
  out = mclapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv(dat, x), mc.cores = numCores)
  
  out = do.call("rbind", out)
  out = cbind(rep(paramsetnum, 4), c(1:4), out)
  colnames(out) = c("paramsetnum", "riskcat",
                    "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                    "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4")
  
  return(out)
}

get_riskcat_ppvnpv_paramset_newinfec = function(paramsetnum){ # output should be a vector of length 99
  
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rawdat = read.csv(paste0(paramset_dir,'ppvnpv_rev.csv'))
  # rawdat2 = read.csv(paste0(paramset_dir,'ppvnpv_tb2.csv'))
  # 
  # rawdat = cbind(rawdat1, rawdat2[,c("lifetimeTB_notReInf", "TBin2years_notReInf")])
  
  temp1 = rawdat %>% dplyr::filter(CountryName %in% c("MEX", "IND", "CHN")) %>%
    slice(rep(1:n(), each = 10))
  
  temp2 = rawdat %>% dplyr::filter(CountryName %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = rawdat %>% dplyr::filter(!(CountryName %in% c("MEX", "IND", "CHN", "PHL")))
  
  dat = rbind(temp1, temp2, temp3)
  
  # out = future_lapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv(dat, x))
  out = mclapply(list(riskcat_I, riskcat_II, riskcat_III, riskcat_IV), function(x) get_riskcat_ppvnpv_newinfec(dat, x), mc.cores = numCores)
  
  out = do.call("rbind", out)
  out = cbind(rep(paramsetnum, 4), c(1:4), out)
  colnames(out) = c("paramsetnum", "riskcat",
                    "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                    "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4",
                    "ppvLT.s2.nri", "ppvLT.s3.nri", "ppvLT.s4.nri", "npvLT.s2.nri", "npvLT.s3.nri", "npvLT.s4.nri",
                    "ppvY2.s2.nri", "ppvY2.s3.nri", "ppvY2.s4.nri", "npvY2.s2.nri", "npvY2.s3.nri", "npvY2.s4.nri")
  
  return(out)
}




if (SA == "PSA_noreinf"){
  # # ppvnpv_allparamsets = future_lapply(1:1000, function(x) get_ppvnpv_paramset(x)) # this takes a while
  ppvnpv_allparamsets = mclapply(1:1000, function(x) get_ppvnpv_paramset(x), mc.cores = numCores) # this takes a while
  ppvnpv_allparamsets = do.call('rbind', ppvnpv_allparamsets)
  colnames(ppvnpv_allparamsets)  = c("ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                                     "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4")
  write.csv(ppvnpv_allparamsets, paste0(CURRENT_RESDIR,"/ppvnpv_tot.csv"), row.names = F)
  
  
  # ppvnpv_riskcat = future_lapply(1:1000, function(x)get_riskcat_ppvnpv_paramset(x))
  ppvnpv_riskcat = mclapply(1:1000, function(x)get_riskcat_ppvnpv_paramset(x), mc.cores=numCores)
  ppvnpv_riskcat = do.call('rbind', ppvnpv_riskcat)
  colnames(ppvnpv_riskcat) = c("paramsetnum", "riskcat",
                               "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                               "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4")
  write.csv(ppvnpv_riskcat, paste0(CURRENT_RESDIR, "/ppvnpv_riskcat.csv"), row.names = F)
  
} else {
  ppvnpv_allparamsets = mclapply(1:1000, function(x) get_ppvnpv_paramset_newinfec(x), mc.cores = numCores) # this takes a while
  
  ppvnpv_allparamsets = do.call('rbind', ppvnpv_allparamsets)
  colnames(ppvnpv_allparamsets) = c("ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
                                    "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4",
                                    "ppvLT.s2.nri", "ppvLT.s3.nri", "ppvLT.s4.nri", "npvLT.s2.nri", "npvLT.s3.nri", "npvLT.s4.nri",
                                    "ppvY2.s2.nri", "ppvY2.s3.nri", "ppvY2.s4.nri", "npvY2.s2.nri", "npvY2.s3.nri", "npvY2.s4.nri")
  write.csv(ppvnpv_allparamsets, paste0(CURRENT_RESDIR,"/ppvnpv_tot_rev.csv"), row.names = F)
  
  # ppvnpv_riskcat = mclapply(1:1000, function(x) get_riskcat_ppvnpv_paramset_newinfec(x), mc.cores=numCores)
  # 
  # ppvnpv_riskcat = do.call('rbind', ppvnpv_riskcat)
  # colnames(ppvnpv_riskcat) = c("paramsetnum", "riskcat",
  #                              "ppvLT.s2", "ppvLT.s3", "ppvLT.s4", "npvLT.s2", "npvLT.s3", "npvLT.s4",
  #                              "ppvY2.s2", "ppvY2.s3", "ppvY2.s4", "npvY2.s2", "npvY2.s3", "npvY2.s4",
  #                              "ppvLT.s2.nri", "ppvLT.s3.nri", "ppvLT.s4.nri", "npvLT.s2.nri", "npvLT.s3.nri", "npvLT.s4.nri",
  #                              "ppvY2.s2.nri", "ppvY2.s3.nri", "ppvY2.s4.nri", "npvY2.s2.nri", "npvY2.s3.nri", "npvY2.s4.nri")
  # write.csv(ppvnpv_riskcat, paste0(CURRENT_RESDIR, "/ppvnpv_riskcat_rev.csv"), row.names = F)
}

# =======================================================================================
# Expenditures
# =======================================================================================


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

get_prod_by_coo = function(rdsfile){
  country_nmb_res = readRDS(rdsfile)
  country_name    = str_sub(rdsfile, -17, -15)
  df = country_nmb_res[["prod"]]
  
  return(df)
}

get_prod_by_coo_df = function(paramsetnum){
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rdsfiles = list.files(paramset_dir, pattern = ".RDS")
  rdsfiles = rdsfiles[-which(rdsfiles %in% c("cea_YUG.RDS", "cea_SUN.RDS"))]
  # temp = future_lapply(rdsfiles, function(x) get_prod_by_coo(paste0(paramset_dir,x)))
  temp = future_lapply(rdsfiles, function(x) get_prod_by_coo(paste0(paramset_dir,x)))
  
  temp = do.call("rbind", temp)
  
  out = temp[,2:4] - temp[,1]
  
  return(out)
}

get_itemizedcosts_by_coo = function(rdsfile, cost.level, varName){
  
  country_nmb_res = readRDS(rdsfile)
  country_name    = str_sub(rdsfile, -17, -15)
  df = country_nmb_res[[cost.level]]
  out = data.frame(strategy = df[, "strategy"], country = country_name, df[, varName])
  colnames(out) = c("strategy", "country", varName)
  
  if (country_name %in% c("MEX", "IND", "CHN")){
    out = out %>% slice(rep(1:n(), each = 10))
  } else if (country_name == "PHL"){
    out = out %>% slice(rep(1:n(), each = 2))
  }
  
  return(out)
}

get_itemizedcosts_df = function(paramsetnum, cost.level, varName){
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rdsfiles = list.files(paramset_dir, pattern = ".RDS")
  rdsfiles = rdsfiles[-which(rdsfiles %in% c("cea_YUG.RDS", "cea_SUN.RDS"))]
  # temp = future_lapply(rdsfiles, function(x) get_itemizedcosts_by_coo(paste0(paramset_dir,x), cost.level, varName))
  temp = future_lapply(rdsfiles, function(x) get_itemizedcosts_by_coo(paste0(paramset_dir,x), cost.level, varName))
  
  temp = do.call("rbind", temp)
  
  out = get_itemized_summary(temp)
  
  s1 = out %>%
    slice(rep(1,3))
  
  out2 = out[2:4,]
  out2[,2:5] = out2[,2:5] - s1[,2:5]
  
  return(out2)
}

get_itemizedcosts_riskcat_df = function(paramsetnum, cost.level, varName, riskcat){
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rdsfiles = list.files(paramset_dir, pattern = ".RDS")
  rdsfiles = rdsfiles[-which(rdsfiles %in% c("cea_YUG.RDS", "cea_SUN.RDS"))]
  # temp = future_lapply(rdsfiles, function(x) get_itemizedcosts_by_coo(paste0(paramset_dir,x), cost.level, varName))
  
  rdsfiles = rdsfiles[str_sub(rdsfiles, 5,7) %in% riskcat]
  
  temp = future_lapply(rdsfiles, function(x) get_itemizedcosts_by_coo(paste0(paramset_dir,x), cost.level, varName))
  
  temp = do.call("rbind", temp)
  
  out = get_itemized_summary(temp)
  
  s1 = out %>%
    slice(rep(1,3))
  
  out2 = out[2:4,]
  out2[,2:5] = out2[,2:5] - s1[,2:5]
  
  return(out2)
}

get_itemized_summary = function(dat){
  
  summary = dat %>% group_by(strategy) %>%
    summarise_at(vars(TBHC_costs.dist:nonHC_costs.dist), mean, na.rm = TRUE)
  
  return(summary)
}



## overall
nmb_by_coo_itm_all = future_lapply(1:1000, function(x) get_itemizedcosts_df(x, "costItems30", c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist")))
saveRDS(nmb_by_coo_itm_all, file = paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_all_30.RDS"))
print("Done writing nmb_by_coo_itm_all_30.RDS")
rm(nmb_by_coo_itm_all)


## subgroup
compile_nmb_by_riskcat = function(costLevel){
  print(paste0("Running nmb_by_coo_itm_risk1_", costLevel,".RDS"))
  nmb_by_coo_itm_risk1 = future_lapply(1:1000, function(x) get_itemizedcosts_riskcat_df(x, costLevel, c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist"), riskcat_I))
  saveRDS(nmb_by_coo_itm_risk1, file = paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk1_",costLevel,".RDS"))
  print(paste0("Done writing nmb_by_coo_itm_risk1_", costLevel,".RDS"))

  print(paste0("Running nmb_by_coo_itm_risk2_", costLevel,".RDS"))
  nmb_by_coo_itm_risk2 = future_lapply(1:1000, function(x) get_itemizedcosts_riskcat_df(x,costLevel, c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist"), riskcat_II))
  saveRDS(nmb_by_coo_itm_risk2, file = paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk2_",costLevel,".RDS"))
  print(paste0("Done writing nmb_by_coo_itm_risk2_", costLevel,".RDS"))

  print(paste0("Running writing nmb_by_coo_itm_risk3_", costLevel,".RDS"))
  nmb_by_coo_itm_risk3 = future_lapply(1:1000, function(x) get_itemizedcosts_riskcat_df(x, costLevel, c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist"), riskcat_III))
  saveRDS(nmb_by_coo_itm_risk3, file = paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk3_",costLevel,".RDS"))
  print(paste0("Done writing nmb_by_coo_itm_risk3_", costLevel,".RDS"))
  
  print(paste0("Running writing nmb_by_coo_itm_risk4_", costLevel,".RDS"))
  nmb_by_coo_itm_risk4 = future_lapply(1:1000, function(x) get_itemizedcosts_riskcat_df(x, costLevel, c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist"), riskcat_IV))
  saveRDS(nmb_by_coo_itm_risk4, file = paste0(CURRENT_RESDIR,"/nmb_by_coo_itm_risk4_",costLevel,".RDS"))
  print(paste0("Done writing nmb_by_coo_itm_risk4_", costLevel,".RDS"))
}


for(i in "costItems15"){ #c("costItems","costItems150", "costItems60", "costItems30", "costItems15")){
  compile_nmb_by_riskcat(i)
}
print("Done with all cost levels")


# =================================================================================================
# Productivity gain
# =================================================================================================

get_prod_by_coo = function(rdsfile){
  country_nmb_res = readRDS(rdsfile)
  country_name    = str_sub(rdsfile, -7, -5)
  df = country_nmb_res[["prod"]]
  
  out = data.frame(strategy = df[, "strategy"], country = country_name, df[, "prod.dist"])
  colnames(out) = c("strategy", "country", "prod.dist")
  
  if (country_name %in% c("MEX", "IND", "CHN")){
    out = out %>% slice(rep(1:n(), each = 10))
  } else if (country_name == "PHL"){
    out = out %>% slice(rep(1:n(), each = 2))
  }
  
  return(out)
}

get_prod_by_coo_df = function(paramsetnum, countries){
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rdsfiles = list.files(paramset_dir, pattern = ".RDS")
  ind = str_sub(rdsfiles, 5, 7)
  rdsfiles = rdsfiles[which(ind %in% countries)]
  # temp = future_lapply(rdsfiles, function(x) get_prod_by_coo(paste0(paramset_dir,x)))
  temp = mclapply(rdsfiles, function(x) get_prod_by_coo(paste0(paramset_dir,x)), mc.cores = numCores)
  
  temp = do.call("rbind", temp)
  temp2 = temp %>% group_by(strategy) %>%
    summarise(mean = mean(prod.dist))
  
  out = temp2 %>% mutate(inc = mean - mean[1]) %>% select(inc)
  
  return(out)
}


get_prod_summary = function(dat){
  
  summary = dat %>% group_by(strategy) %>%
    summarise_at(vars(prod.dist), mean)
  
  return(summary)
}


## overall
prod_by_coo_all = future_lapply(1:1000, function(x) get_prod_by_coo_df(x, countrynames))
prod_by_coo_all_array = array(unlist(prod_by_coo_all), dim = c(97, 3, 1000))
saveRDS(prod_by_coo_all_array, file = paste0(CURRENT_RESDIR,"/prod_by_coo_all_array.RDS"))

pop_size = cohort2019_by_coo2$total_pop[order(cohort2019_by_coo2$place_of_birth)]
prod_by_strat_mean = apply(prod_by_coo_all_array, c(2:3), function(x) weighted.mean(x, pop_size/sum(pop_size)))
prod_by_strat_mean = apply(prod_by_strat_mean, 1, function(x) mean(x))

## subgroup
v = cohort2019_by_coo2$place_of_birth %in% riskcat_I
prod_by_strat_mean_risk1 = apply(prod_by_coo_all_array[v,,], c(2:3), function(x) weighted.mean(x, pop_size[v]/sum(pop_size[v])))
prod_by_strat_mean_risk1 = apply(prod_by_strat_mean_risk1, 1, function(x) mean(x))

v = cohort2019_by_coo2$place_of_birth %in% riskcat_II
prod_by_strat_mean_risk2 = apply(prod_by_coo_all_array[v,,], c(2:3), function(x) weighted.mean(x, pop_size[v]/sum(pop_size[v])))
prod_by_strat_mean_risk2 = apply(prod_by_strat_mean_risk2, 1, function(x) mean(x))

v = cohort2019_by_coo2$place_of_birth %in% riskcat_III
prod_by_strat_mean_risk3 = apply(prod_by_coo_all_array[v,,], c(2:3), function(x) weighted.mean(x, pop_size[v]/sum(pop_size[v])))
prod_by_strat_mean_risk3 = apply(prod_by_strat_mean_risk3, 1, function(x) mean(x))

v = cohort2019_by_coo2$place_of_birth %in% riskcat_IV
prod_by_strat_mean_risk4 = apply(prod_by_coo_all_array[v,,], c(2:3), function(x) weighted.mean(x, pop_size[v]/sum(pop_size[v])))
prod_by_strat_mean_risk4 = apply(prod_by_strat_mean_risk4, 1, function(x) mean(x))

prod_df = rbind(prod_by_strat_mean_risk1, prod_by_strat_mean_risk2, prod_by_strat_mean_risk3, prod_by_strat_mean_risk4, prod_by_strat_mean)
colnames(prod_df) = c("s2", "s3", "s4")
write.csv(prod_df, paste0(CURRENT_RESDIR,"/prod_df.csv"), row.names = rownames(prod_df))


print("Running prod_by_coo_all.RDS ")
prod_by_coo_all = future_lapply(1:1000, function(x) get_prod_by_coo_df(x, countrynames))
prod_by_coo_all = do.call("cbind", prod_by_coo_all)
colnames(prod_by_coo_all) = c(1:1000)
saveRDS(prod_by_coo_all, file = paste0(CURRENT_RESDIR,"/prod_by_coo_all.RDS"))


print("Running prod_by_coo_all_risk1.RDS ")
prod_by_coo_all = future_lapply(1:1000, function(x) get_prod_by_coo_df(x, riskcat_I))
prod_by_coo_all = do.call("cbind", prod_by_coo_all)
colnames(prod_by_coo_all) = c(1:1000)
saveRDS(prod_by_coo_all, file = paste0(CURRENT_RESDIR,"/prod_by_coo_risk1.RDS"))
#
# rm(prod_by_coo_all)
#
print("Running prod_by_coo_all_risk2.RDS ")
prod_by_coo_all = future_lapply(1:1000, function(x) get_prod_by_coo_df(x, riskcat_II))
prod_by_coo_all = do.call("cbind", prod_by_coo_all)
colnames(prod_by_coo_all) = c(1:1000)
saveRDS(prod_by_coo_all, file = paste0(CURRENT_RESDIR,"/prod_by_coo_risk2.RDS"))
#
# rm(prod_by_coo_all)
#
print("Running prod_by_coo_all_risk3.RDS ")
prod_by_coo_all = future_lapply(1:1000, function(x) get_prod_by_coo_df(x, riskcat_III))
prod_by_coo_all = do.call("cbind", prod_by_coo_all)
colnames(prod_by_coo_all) = c(1:1000)
saveRDS(prod_by_coo_all, file = paste0(CURRENT_RESDIR,"/prod_by_coo_risk3.RDS"))
#
# rm(prod_by_coo_all)
#
print("Running prod_by_coo_all_risk4.RDS ")
prod_by_coo_all = future_lapply(1:1000, function(x) get_prod_by_coo_df(x, riskcat_IV))
prod_by_coo_all = do.call("cbind", prod_by_coo_all)
colnames(prod_by_coo_all) = c(1:1000)
saveRDS(prod_by_coo_all, file = paste0(CURRENT_RESDIR,"/prod_by_coo_risk4.RDS"))


# =================================================================================================
# Effectiveness
# =================================================================================================

get_eff_by_coo = function(rdsfile){
  country_nmb_res = readRDS(rdsfile)
  # country_name    = str_sub(rdsfile, -18, -16)
  country_name    = str_sub(rdsfile, nchar(rdsfile)-6, nchar(rdsfile)-4)
  
  df = country_nmb_res[["eff.dQALY"]]
  
  out = data.frame(strategy = df[, "strategy"], country = country_name, disQALY_sum = df[, "disQALY_sum"])
  
  if (country_name %in% c("MEX", "IND", "CHN")){
    out = out %>% slice(rep(1:n(), each = 10))
  } else if (country_name == "PHL"){
    out = out %>% slice(rep(1:n(), each = 2))
  }
  
  return(out)
}

get_eff_by_coo_df = function(paramsetnum, countries){
  paramset_dir = paste0(CURRENT_RESDIR,"/paramset",paramsetnum,"/")
  rdsfiles = list.files(paramset_dir, pattern = ".RDS")
  # ind = str_sub(rdsfiles, -18, -16)
  ind =str_sub(paste0(paramset_dir,rdsfiles), nchar(paste0(paramset_dir,rdsfiles))-6, nchar(paste0(paramset_dir,rdsfiles))-4)
  
  rdsfiles = rdsfiles[which(ind %in% countries)]
  # temp = future_lapply(rdsfiles, function(x) get_eff_by_coo(paste0(paramset_dir,x)))
  temp = mclapply(rdsfiles, function(x) get_eff_by_coo(paste0(paramset_dir,x)), mc.cores = numCores)
  
  temp = do.call("rbind", temp)
  temp2 = temp %>% group_by(strategy) %>%
    summarise(mean = mean(disQALY_sum))
  
  out = temp2 %>% mutate(inc = mean - mean[1]) %>% select(inc)
  
  return(out)
}

countrynames = cohort2019_by_coo2$place_of_birth[order(cohort2019_by_coo2$place_of_birth)]


print("eff_by_coo_all starting...")
eff_by_coo_all = future_lapply(1:1000, function(x) get_eff_by_coo_df(x, countrynames))
saveRDS(eff_by_coo_all, file = paste0(CURRENT_RESDIR,"/eff_by_coo_all.RDS"))
rm(eff_by_coo_all)

print("eff_by_coo_risk1.RDS starting...")
eff_by_coo_all = future_lapply(1:1000, function(x) get_eff_by_coo_df(x, riskcat_I))
eff_by_coo_all = do.call("cbind", eff_by_coo_all)
colnames(eff_by_coo_all) = c(1:1000)
saveRDS(eff_by_coo_all, file = paste0(CURRENT_RESDIR,"/eff_by_coo_risk1.RDS"))
rm(eff_by_coo_all)

print("eff_by_coo_risk2.RDS starting...")
eff_by_coo_all = future_lapply(1:1000, function(x) get_eff_by_coo_df(x, riskcat_II))
eff_by_coo_all = do.call("cbind", eff_by_coo_all)
colnames(eff_by_coo_all) = c(1:1000)
saveRDS(eff_by_coo_all, file = paste0(CURRENT_RESDIR,"/eff_by_coo_risk2.RDS"))
rm(eff_by_coo_all)

print("eff_by_coo_risk3.RDS starting...")
eff_by_coo_all = future_lapply(1:1000, function(x) get_eff_by_coo_df(x, riskcat_III))
eff_by_coo_all = do.call("cbind", eff_by_coo_all)
colnames(eff_by_coo_all) = c(1:1000)
saveRDS(eff_by_coo_all, file = paste0(CURRENT_RESDIR,"/eff_by_coo_risk3.RDS"))
rm(eff_by_coo_all)

print("eff_by_coo_risk4.RDS starting...")
eff_by_coo_all = future_lapply(1:1000, function(x) get_eff_by_coo_df(x, riskcat_IV))
eff_by_coo_all = do.call("cbind", eff_by_coo_all)
colnames(eff_by_coo_all) = c(1:1000)
saveRDS(eff_by_coo_all, file = paste0(CURRENT_RESDIR,"/eff_by_coo_risk4.RDS"))



