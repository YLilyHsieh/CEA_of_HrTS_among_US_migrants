#!/usr/bin/env Rscript

PARAMset_num  = commandArgs(trailingOnly=TRUE)[1]
PARAMset_num = as.numeric(PARAMset_num)
print(paste0("PARAMset_num = ",PARAMset_num, "; class = ", class(PARAMset_num)))

SA = commandArgs(trailingOnly=TRUE)[2]
SA = as.character(SA)
print(paste0("SA = ", SA, "; class = ", class(SA)))

## Load settings and functions ================================================
PSA = TRUE

# Read in the PSAtable  # created in psa_helpers.R
# PSAtable = read.csv("Data/PSAtable80.csv")
PSAtable = read.csv("Data/PSAtable1000.csv")

# Subset the PSAtable to rows included in this slurm job
PARAMset = PSAtable[PSAtable$n == PARAMset_num,]; rm(PSAtable)

# Read in the sim_preds0 table # created in data.R
sim_preds0 = read.csv("Data/sim_preds0_1000_calib.csv")
predTBcases = sim_preds0[,PARAMset_num]; rm(sim_preds0)

# Call helper functions and initiate the cohort
source("Model_PLoSMed_rev/init_cohort.R")
source("Model_PLoSMed_rev/sim_helpers.R")

# library(future.apply)
# plan(multisession)

# Point to the folder that will store results for this parameter set
ResultsDir_PSA_current = paste0("/n/netscratch/menzies_lab/Everyone/ylh202/IncipienTB/",SA,"/paramset",PARAMset_num,"/")


get_ppvnpv_df = function(CountryName){
  
  if (CountryName %in% c("MEX","CHN","IND")){
    fraction = 0.1
  } else if (CountryName == "PHL") {
    fraction = 0.5
  } else {
    fraction = 1
  }
  
  # Subset data to the country-of-origin of interest
  cohort_df_subset = NULL
  cohort_df_subset = cohort_df[cohort_df$place_of_birth == CountryName,]
  
  # Create a dataframe of agents for that country-of-origin
  cohort = create_cohort(CountryName = CountryName, cohort_raw = cohort_df_subset, fraction = fraction, PSA = PSA)
  
  temp = as.data.frame(matrix(nrow=nrow(cohort),ncol=1)); colnames(temp) = "lifetimeTB"
  temp$lifetimeTB = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] > cohort$time_to_TB[x], 1, 0))
  temp$TBin2years = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] > cohort$time_to_TB[x] & cohort$time_to_TB[x] <= 25/12, 1, 0)) # 2 years from testing
  temp$rna_time.s3 = 1/12 + 1/52 
  temp$rna_time.s4 = 1/12
  
  
  # Strategy II
  temp.s2 = temp
  temp.s2$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12, "TB", cohort$TB_status[x]))
  temp.s2$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] < 1/12, "Dead", cohort$TB_status[x]))
  temp.s2$igra_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s2$healthstate_igra[x] == "Dead", NA, 
                                                                get_igra_result(temp.s2$healthstate_igra[x],
                                                                                IGRA_SPEC_HEALTHY, IGRA_SENS_LTBI, IGRA_SENS_TB, cohort$rnum.igra[x])))
  temp.s2$test.res = ifelse(temp.s2$igra_pos, 1, 0)
  # temp.s2$healthstate_xray = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x] < 1/12 + 1/52, "TB", temp.s2$healthstate_igra[x]))
  # temp.s2$xray_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s2$healthstate_xray[x] == "DEAD", NA,
  #                                                               get_xray_result(temp.s2$healthstate_xray[x], 
  #                                                                               XRAY_SPEC_HEALTHY, XRAY_SPEC_LTBI, XRAY_SENS_TB,cohort$rnum.xray[x])))
  
  # Strategy III
  temp.s3 = temp.s2 %>% select(-test.res)
  temp.s3$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12, "TB", cohort$TB_status[x]))
  temp.s3$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] < 1/12, "Dead", cohort$TB_status[x]))
  
  temp.s3$healthstate_rna  = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12 + 1/52, "TB", temp.s3$healthstate_igra[x]))
  temp.s3$rna_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s3$healthstate_igra[x] == "Dead", NA,
                                                               get_rna_result(health_status = temp.s3$healthstate_rna[x],
                                                                              time_to_tb = cohort$time_to_TB[x] - temp$rna_time.s3[x], 
                                                                              agent_rnum = cohort$rnum.rna[x])))
  temp.s3$rna_pos = ifelse(temp.s3$igra_pos == F, F, temp.s3$rna_pos)
  
  temp.s3$test.res = ifelse(temp.s3$rna_pos, 1, 0)
  # temp.s3$healthstate_xray = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x] < 1/12 + 1/52*2, "TB", temp.s3$healthstate_rna[x]))
  # temp.s3$xray_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s3$healthstate_xray == "DEAD", NA,
  #                                                               get_xray_result(temp.s3$healthstate_xray[x],
  #                                                                               XRAY_SPEC_HEALTHY, XRAY_SPEC_LTBI, XRAY_SENS_TB,cohort$rnum.xray[x])))
  
  
  # Strategy IV
  temp.s4 = temp
  temp.s4$healthstate_rna = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12, "TB", cohort$TB_status[x]))
  temp.s4$healthstate_rna = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] < 1/12, "Dead", cohort$TB_status[x]))
  temp.s4$rna_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s4$healthstate_rna[x] == "Dead", NA,
                                                               get_rna_result(health_status = temp.s4$healthstate_rna[x],
                                                                              time_to_tb = cohort$time_to_TB[x] - temp$rna_time.s4[x], 
                                                                              agent_rnum = cohort$rnum.rna[x])))
  temp.s4$test.res = ifelse(temp.s4$rna_pos, 1, 0)
  # temp.s4$healthstate_xray = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x] < 1/12 + 1/52, "TB", temp.s4$healthstate_rna[x]))
  # temp.s4$xray_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s4$healthstate_xray == "DEAD", NA,
  #                                                               get_xray_result(temp.s4$healthstate_xray[x],
  #                                                                               XRAY_SPEC_HEALTHY, XRAY_SPEC_LTBI, XRAY_SENS_TB,cohort$rnum.xray[x])))
  
  df_out = data.frame(matrix(NA, nrow = nrow(cohort), ncol = 6))
  df_out = data.frame(
    PARAMset_num = rep(PARAMset_num, nrow(cohort)),
    CountryName  = rep(CountryName, nrow(cohort)),
    lifetimeTB   = temp.s2$lifetimeTB,
    TBin2years   = temp.s2$TBin2years,
    testres.s2   = temp.s2$test.res,
    testres.s3   = temp.s3$test.res,
    testres.s4   = temp.s4$test.res)
  print(paste0(CountryName, ": Done!"))
  return(df_out)
}



get_ppvnpv_mod_df = function(CountryName){
  
  if (CountryName %in% c("MEX","CHN","IND")){
    fraction = 0.1
  } else if (CountryName == "PHL") {
    fraction = 0.5
  } else {
    fraction = 1
  }
  
  # Subset data to the country-of-origin of interest
  cohort_df_subset = NULL
  cohort_df_subset = cohort_df[cohort_df$place_of_birth == CountryName,]
  
  # Create a dataframe of agents for that country-of-origin
  cohort = create_cohort(CountryName = CountryName, cohort_raw = cohort_df_subset, fraction = fraction, PSA = PSA)
  
  # temp = as.data.frame(matrix(nrow=nrow(cohort),ncol=1)); colnames(temp) = "lifetimeTB"
  temp = cohort
  temp$lifetimeTB = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] > cohort$time_to_TB[x], 1, 0))
  temp$TBin2years = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] > cohort$time_to_TB[x] & cohort$time_to_TB[x] <= 25/12, 1, 0)) # 2 years from testing
  temp$lifetimeTB_notReInf = sapply(1:nrow(cohort), function(x) ifelse(temp$lifetimeTB[x] == 1 & cohort$tb.reinfection[x] == 0, 1, 0))
  temp$TBin2years_notReInf = sapply(1:nrow(cohort), function(x) ifelse(temp$TBin2years[x] == 1 & cohort$tb.reinfection[x] == 0, 1, 0)) # 2 years from testing
  
  # temp$TB_status_mod = sapply(1:nrow(cohort), function(x) ifelse(temp$lifetimeTB_notReInf[x], "HEALTHY" , temp$TB_status[x]))
  # temp$TB_status_mod = sapply(1:nrow(cohort), function(x) ifelse(temp$lifetimeTB_notReInf[x],temp$TB_status[x]  , "HEALTHY"))
  temp$TB_status_mod = temp$TB_status
  
  temp$rna_time.s3 = 1/12 + 1/52 
  temp$rna_time.s4 = 1/12
  
  
  # Strategy II
  temp.s2 = temp
  temp.s2$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12, "TB", temp$TB_status_mod[x]))
  temp.s2$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] < 1/12, "Dead", temp$TB_status_mod[x]))
  temp.s2$igra_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s2$healthstate_igra[x] == "Dead", NA, 
                                                                get_igra_result(temp.s2$healthstate_igra[x],
                                                                                IGRA_SPEC_HEALTHY, IGRA_SENS_LTBI, IGRA_SENS_TB, cohort$rnum.igra[x])))
  temp.s2$test.res_mod = ifelse(temp.s2$igra_pos, 1, 0)
  
  
  # Strategy III
  temp.s3 = temp.s2 %>% select(-test.res_mod)
  temp.s3$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12, "TB", temp$TB_status_mod[x]))
  temp.s3$healthstate_igra = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] < 1/12, "Dead", temp$TB_status_mod[x]))
  
  temp.s3$healthstate_rna  = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12 + 1/52, "TB", temp.s3$healthstate_igra[x]))
  temp.s3$rna_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s3$healthstate_igra[x] == "Dead", NA,
                                                               get_rna_result(health_status = temp.s3$healthstate_rna[x],
                                                                              time_to_tb = cohort$time_to_TB[x] - temp$rna_time.s3[x], 
                                                                              agent_rnum = cohort$rnum.rna[x])))
  temp.s3$rna_pos = ifelse(temp.s3$igra_pos == F, F, temp.s3$rna_pos)
  
  temp.s3$test.res_mod = ifelse(temp.s3$rna_pos, 1, 0)
  
  
  # Strategy IV
  temp.s4 = temp
  temp.s4$healthstate_rna = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_TB[x]  < 1/12, "TB", temp$TB_status_mod[x]))
  temp.s4$healthstate_rna = sapply(1:nrow(cohort), function(x) ifelse(cohort$time_to_death[x] < 1/12, "Dead", temp$TB_status_mod[x]))
  temp.s4$rna_pos = sapply(1:nrow(cohort), function(x) ifelse (temp.s4$healthstate_rna[x] == "Dead", NA,
                                                               get_rna_result(health_status = temp.s4$healthstate_rna[x],
                                                                              time_to_tb = cohort$time_to_TB[x] - temp$rna_time.s4[x], 
                                                                              agent_rnum = cohort$rnum.rna[x])))
  temp.s4$test.res_mod = ifelse(temp.s4$rna_pos, 1, 0)
  
  df_out = data.frame(matrix(NA, nrow = nrow(cohort), ncol = 11))
  df_out = data.frame(
    PARAMset_num = rep(PARAMset_num, nrow(cohort)),
    CountryName  = rep(CountryName, nrow(cohort)),
    lifetimeTB   = temp$lifetimeTB,
    TBin2years   = temp$TBin2years,
    lifetimeTB_notReInf = temp$lifetimeTB_notReInf,
    TBin2years_notReInf = temp$TBin2years_notReInf,
    TB_status = temp$TB_status,
    TB_status_mod = temp$TB_status_mod,
    testres.s2_mod   = temp.s2$test.res_mod,
    testres.s3_mod   = temp.s3$test.res_mod,
    testres.s4_mod   = temp.s4$test.res_mod)
  
  
  print(paste0(CountryName, ": Done!"))
  return(df_out)
}

# out = mclapply(CountryNames, function(x) get_ppvnpv_df(x), mc.cores= numCores)
# out = mclapply(CountryNames, function(x) get_TBcounts_df(x), mc.cores= 12)

if (SA == "PSA_noreinf"){
  out = mclapply(CountryNames, function(x) get_ppvnpv_df(x), mc.cores = numCores)
} else {
  out = mclapply(CountryNames, function(x) get_ppvnpv_mod_df(x), mc.cores = numCores)
}

out = do.call("rbind", out)
# write.csv(out, paste0(ResultsDir_PSA_current,"ppvnpv.csv"), row.names = F)
# write.csv(out, paste0(ResultsDir_PSA_current,"ppvnpv_tb.csv"), row.names = F) # tb.reinfection is incorrect
# write.csv(out, paste0(ResultsDir_PSA_current,"ppvnpv_tb2.csv"), row.names = F) # 
write.csv(out, paste0(ResultsDir_PSA_current,"ppvnpv_rev.csv"), row.names = F) 

print(paste0("paramset", PARAMset_num, "ppvnpv_rev written"))
## ============================================================================
# Description: 
# A stand alone script to read through init cohort files
# to determine NPV and PPV for Strategies II, III, IV 

# The script is built upon sections of psa.R

