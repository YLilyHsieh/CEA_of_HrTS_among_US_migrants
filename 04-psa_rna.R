#!/usr/bin/env Rscript

PARAMset_num  = commandArgs(trailingOnly=TRUE)
PARAMset_num = as.numeric(PARAMset_num)
print(paste0("PARAMset_num = ",PARAMset_num, "; class = ", class(PARAMset_num)))


## Load settings and functions ================================================
PSA = TRUE

# Read in the PSAtable  # created in psa_helpers.R
PSAtable = read.csv("Data/PSAtable1000_gam.csv")


# Subset the PSAtable to rows included in this slurm job
PARAMset = PSAtable[PSAtable$n == PARAMset_num,]; rm(PSAtable)

# # Read in the sim_preds0 table # created in data.R
sim_preds0 = read.csv("Data/sim_preds0_1000_calib.csv")
predTBcases = sim_preds0[,PARAMset_num]; rm(sim_preds0)


# Call helper functions and initiate the cohort
source("Model/01-init_cohort.R")
source("Model/02-sim_helpers.R")
source("Model/02-res_helpers.R")

# Point to the folder that will store results for this parameter set
# ResultsDir_PSA_current = paste0(ResultsDir_PSA,"paramset",PARAMset_num,"/")
# ResultsDir_PSA_current = paste0("/n/holyscratch01/menzies_lab/Users/ylh202/IncipienTB/PSA_calib_testscenar1/paramset",PARAMset_num,"/")
ResultsDir_PSA_current = paste0("~/IncipienTB/PSA_RNAgam/paramset",PARAMset_num,"/")

## Simulate and collect model output ==========================================

StrategyNums = c(0, 1, 2, 3)
STRATEGY_NAMES = STRATEGY_NAMES_complete[StrategyNums + 1]


run_model = function(input, name, strategies){
  res = NULL
  print(paste(name,"START..."))
  # print(paste("before running: ", mem_used()))
  
  res = get_sim_results(input, strategies)
  print(paste(name,"DONE..."))
  # file.out = paste0(ResultsDir,"res_",name,".RDS")
  # print(paste("after running: ", mem_used()))
  # print(paste("after gc: ", mem_used()))
  return(res)
}

# run_model2 = function(name, strategies){
#   res = NULL
#   print(paste(name,"START..."))
#   # print(paste("before running: ", mem_used()))
#   input = init_country_cohort(name, callSavedFiles = T)
#   
#   res = get_sim_results(input, strategies)
#   print(paste(name,"DONE..."))
#   # file.out = paste0(ResultsDir,"res_",name,".RDS")
#   file.out = paste0(ResultsDir_PSA_current,"res_",name,".RDS")
#   saveRDS(res, file.out)
#   
#   print(paste0(file.out, "...written."))
#   gc()
#   # print(paste("after running: ", mem_used()))
#   # print(paste("after gc: ", mem_used()))
#   # return(res)
# }

get_psa_res = function(name, strategies){
  
  # simulate the model --------------------------------------------------------
  input = init_country_cohort(name,  PSA = T)
  country_res = run_model(input, name, strategies)
  
  # # get time to TB for each strategy for the country --------------------------
  # time_to_tb = get_country_TTE(country_res = country_res, country_code = name, health_outcome = "TB")
  # time_to_tb = augmented_TTE(time_to_tb)
  # write.csv(time_to_tb, file = paste0(ResultsDir_PSA_current,"time2TB_",name,".csv"), row.names = F)
  # 
  # # get time to TB diagnosis/treatment for each strategy for the country -------
  # time_to_tbtx = get_country_TTE(country_res = country_res, country_code = name, health_outcome = c("TREAT.TB.GUARANTEE", "TREAT.TB_INIT.SUCCESS"))
  # time_to_tbtx = augmented_TTE(time_to_tbtx)
  # write.csv(time_to_tbtx, file = paste0(ResultsDir_PSA_current,"time2tbtx_",name,".csv"), row.names = F)
  # 
  # # get time to DEATH for each strategy for the country --------------------------
  # time_to_death = get_country_TTE(country_res = country_res, country_code = name, health_outcome = ALL.CAUSE.DEATH)
  # time_to_death = augmented_TTE(time_to_death)
  # write.csv(time_to_death, file = paste0(ResultsDir_PSA_current,"time2death_",name,".csv"), row.names = F)
  
  # get event tabulation for the country --------------------------------------
  event_tab = get_country_event_tab(country_res = country_res, country_code = name)
  event_tab = augmented_tab(event_tab)
  write.csv(event_tab, file = paste0(ResultsDir_PSA_current,"event_tab_",name,".csv"), row.names = F)
  
  # EFFECT =====================================================================
  # get mean eff (discounted QALY) for the country ----------------------------
  eff.dQALY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input, r = 0.03)
  names(eff.dQALY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # get mean eff (undiscounted QALY) for the country ----------------------------
  eff.QALY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input, r = 0)
  names(eff.QALY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB",  "RNA_TB")
  
  # get mean eff (discounted LY) for the country ------------------------------------
  reward_input2 = ifelse(names(reward_input) == "DEAD", 0, 1)
  names(reward_input2) = names(reward_input)
  eff.dLY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input2, r = 0.03)
  names(eff.dLY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # get mean eff (undiscounted LY) for the country ------------------------------------
  eff.LY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input2, r = 0)
  names(eff.LY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # PRODUCTIVITY ===============================================================
  prod_ctlg = PRODUCTIVITY
  prod = get_mean_prod(country_res, strategy_nums = StrategyNums, prod_ctlg, r = 0.03)
  names(prod) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # COSTS ======================================================================
  
  ## healthcare sector perspective =============================================
  # get mean cost (discounted cost) for each country ---------------------------
  costHC = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = TRUE)
  names(costHC) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # scenario analysis for RNA cost ---------------------------------------------
  TBCOSTS_HC["TEST.RNA"] = 150
  costHC150 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = TRUE)
  names(costHC150) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  TBCOSTS_HC["TEST.RNA"] = 60
  costHC60 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = TRUE)
  names(costHC60) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  TBCOSTS_HC["TEST.RNA"] = 30
  costHC30 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = TRUE)
  names(costHC30) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  TBCOSTS_HC["TEST.RNA"] = 15
  costHC15 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = TRUE)
  names(costHC15) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  ## societal perspective ======================================================
  # get mean cost (discounted cost) for each country ---------------------------
  TBCOSTS_HC["TEST.RNA"] = 300
  costSOC = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = FALSE)
  names(costSOC) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # scenario analysis for RNA cost ---------------------------------------------
  TBCOSTS_HC["TEST.RNA"] = 150
  costSOC150 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = FALSE)
  names(costSOC150) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  TBCOSTS_HC["TEST.RNA"] = 60
  costSOC60 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = FALSE)
  names(costSOC60) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  TBCOSTS_HC["TEST.RNA"] = 30
  costSOC30 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = FALSE)
  names(costSOC30) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  TBCOSTS_HC["TEST.RNA"] = 15
  costSOC15 = get_mean_cost(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03, HCperspective = FALSE)
  names(costSOC15) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  
  # # calculate NMB for each strategy -------------------------------------------
  # NMB.No_screening = eff[["No_testing"]]*WTP - cost[["No_testing"]]
  # NMB.IGRA_TB      = eff[["IGRA_TB"]]*WTP - cost[["IGRA_TB"]]
  # NMB.IGRA_RNA_TB  = eff[["IGRA_RNA_TB"]]*WTP - cost[["IGRA_RNA_TB"]]
  # iNMB.IGRA_TB = (eff[["IGRA_TB"]]*WTP - cost[["IGRA_TB"]]) - (eff[["No_testing"]]*WTP - cost[["No_testing"]])
  # iNMB.IGRA_RNA_TB =  (eff[["IGRA_RNA_TB"]]*WTP - cost[["IGRA_RNA_TB"]]) - (eff[["No_testing"]]*WTP - cost[["No_testing"]])
  
  # collect results into a list -----------------------------------------------
  cea = list(eff.dQALY = eff.dQALY,
             eff.QALY  = eff.QALY,
             eff.dLY   = eff.dLY,
             eff.LY    = eff.LY,
             prod      = prod,
             costHC = costHC,
             costHC150 = costHC150,
             costHC60  = costHC60,
             costHC30  = costHC30,
             costHC15  = costHC15,
             costSOC = costSOC,
             costSOC150 = costSOC150,
             costSOC60  = costSOC60,
             costSOC30  = costSOC30,
             costSOC15  = costSOC15)
  # NMB.No_screening = NMB.No_screening,
  # NMB.IGRA_TB = NMB.IGRA_TB,
  # NMB.IGRA_RNA_TB = NMB.IGRA_RNA_TB,
  # iNMB.IGRA_TB = iNMB.IGRA_TB,
  # iNMB.IGRA_RNA_TB = iNMB.IGRA_RNA_TB)
  
  # save cea result to RDS file -----------------------------------------------
  file.out2 = paste0(ResultsDir_PSA_current,"cea_",name,".RDS")
  saveRDS(cea, file.out2)
  
  print(paste0(file.out2, "...written."))
  print(mem_used())
  #mem_tracker = append(mem_tracker, mem_used())
  gc()
}

get_psa_res3 = function(name, strategies){
  
  # simulate the model --------------------------------------------------------
  input = init_country_cohort(name,  PSA = T)
  country_res = run_model(input, name, strategies)
  
  # # get time to TB for each strategy for the country --------------------------
  # time_to_tb = get_country_TTE(country_res = country_res, country_code = name, health_outcome = "TB")
  # time_to_tb = augmented_TTE(time_to_tb)
  # write.csv(time_to_tb, file = paste0(ResultsDir_PSA_current,"time2TB_",name,".csv"), row.names = F)
  # 
  # # get time to TB diagnosis/treatment for each strategy for the country -------
  # time_to_tbtx = get_country_TTE(country_res = country_res, country_code = name, health_outcome = c("TREAT.TB.GUARANTEE", "TREAT.TB_INIT.SUCCESS"))
  # time_to_tbtx = augmented_TTE(time_to_tbtx)
  # write.csv(time_to_tbtx, file = paste0(ResultsDir_PSA_current,"time2tbtx_",name,".csv"), row.names = F)
  # 
  # # get time to DEATH for each strategy for the country --------------------------
  # time_to_death = get_country_TTE(country_res = country_res, country_code = name, health_outcome = ALL.CAUSE.DEATH)
  # time_to_death = augmented_TTE(time_to_death)
  # write.csv(time_to_death, file = paste0(ResultsDir_PSA_current,"time2death_",name,".csv"), row.names = F)
  # 
  # # get event tabulation for the country --------------------------------------
  # event_tab = get_country_event_tab(country_res = country_res, country_code = name)
  # event_tab = augmented_tab(event_tab)
  # write.csv(event_tab, file = paste0(ResultsDir_PSA_current,"event_tab_",name,".csv"), row.names = F)
  
  # EFFECT =====================================================================
  # get mean eff (discounted QALY) for the country ----------------------------
  eff.dQALY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input, r = 0.03)
  names(eff.dQALY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # get mean eff (undiscounted QALY) for the country ----------------------------
  eff.QALY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input, r = 0)
  names(eff.QALY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB",  "RNA_TB")
  
  # get mean eff (discounted LY) for the country ------------------------------------
  reward_input2 = ifelse(names(reward_input) == "DEAD", 0, 1)
  names(reward_input2) = names(reward_input)
  eff.dLY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input2, r = 0.03)
  names(eff.dLY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # get mean eff (undiscounted LY) for the country ------------------------------------
  eff.LY = get_mean_eff(country_res, strategy_nums = StrategyNums, reward_input2, r = 0)
  names(eff.LY) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # PRODUCTIVITY ===============================================================
  prod_ctlg = PRODUCTIVITY
  prod = get_mean_prod(country_res, strategy_nums = StrategyNums, prod_ctlg, r = 0.03)
  names(prod) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  
  # COSTS ======================================================================
  
  # get mean cost (discounted cost) for each country ---------------------------
  costItems = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  # scenario analysis for RNA cost ---------------------------------------------
  TBCOSTS_HC["TEST.RNA"] = 150
  costItems150 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  TBCOSTS_HC["TEST.RNA"] = 60
  costItems60 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  TBCOSTS_HC["TEST.RNA"] = 30
  costItems30 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  TBCOSTS_HC["TEST.RNA"] = 15
  costItems15 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  
  # collect results into a list -----------------------------------------------
  cea = list(eff.dQALY = eff.dQALY,
             eff.QALY  = eff.QALY,
             eff.dLY   = eff.dLY,
             eff.LY    = eff.LY,
             prod      = prod,
             costItems = costItems,
             costItems150 = costItems150,
             costItems60  = costItems60,
             costItems30  = costItems30,
             costItems15  = costItems15)
  # NMB.No_screening = NMB.No_screening,
  # NMB.IGRA_TB = NMB.IGRA_TB,
  # NMB.IGRA_RNA_TB = NMB.IGRA_RNA_TB,
  # iNMB.IGRA_TB = iNMB.IGRA_TB,
  # iNMB.IGRA_RNA_TB = iNMB.IGRA_RNA_TB)
  
  # save cea result to RDS file -----------------------------------------------
  file.out2 = paste0(ResultsDir_PSA_current,"cea_",name,"_costItems_dft5.RDS")
  saveRDS(cea, file.out2)
  
  print(paste0(file.out2, "...written."))
  print(mem_used())
  #mem_tracker = append(mem_tracker, mem_used())
  gc()
}

get_psa_res4 = function(name, strategies){
  
  # simulate the model --------------------------------------------------------
  input = init_country_cohort(name,  PSA = T)
  country_res = run_model(input, name, strategies)
  
  # # get time to TB for each strategy for the country --------------------------
  # time_to_tb = get_country_TTE(country_res = country_res, country_code = name, health_outcome = "TB")
  # time_to_tb = augmented_TTE(time_to_tb)
  # write.csv(time_to_tb, file = paste0(ResultsDir_PSA_current,"time2TB_",name,".csv"), row.names = F)
  # 
  # # get time to TB diagnosis/treatment for each strategy for the country -------
  # time_to_tbtx = get_country_TTE(country_res = country_res, country_code = name, health_outcome = c("TREAT.TB.GUARANTEE", "TREAT.TB_INIT.SUCCESS"))
  # time_to_tbtx = augmented_TTE(time_to_tbtx)
  # write.csv(time_to_tbtx, file = paste0(ResultsDir_PSA_current,"time2tbtx_",name,".csv"), row.names = F)
  # 
  # # get time to DEATH for each strategy for the country --------------------------
  # time_to_death = get_country_TTE(country_res = country_res, country_code = name, health_outcome = ALL.CAUSE.DEATH)
  # time_to_death = augmented_TTE(time_to_death)
  # write.csv(time_to_death, file = paste0(ResultsDir_PSA_current,"time2death_",name,".csv"), row.names = F)
  # 
  # get event tabulation for the country --------------------------------------
  event_tab = get_country_event_tab(country_res = country_res, country_code = name)
  event_tab = augmented_tab(event_tab)
  write.csv(event_tab, file = paste0(ResultsDir_PSA_current,"event_tab_",name,".csv"), row.names = F)

  # EFFECT =====================================================================
  # get mean eff (discounted QALY) for the country ----------------------------
  eff.dQALY = get_mean_eff2(country_res, strategy_nums = StrategyNums, reward_input, r = 0.03)
  
  # get mean eff (undiscounted QALY) for the country ----------------------------
  eff.QALY = get_mean_eff2(country_res, strategy_nums = StrategyNums, reward_input, r = 0)
  
  # get mean eff (discounted LY) for the country ------------------------------------
  reward_input2 = ifelse(names(reward_input) == "DEAD", 0, 1)
  names(reward_input2) = names(reward_input)
  eff.dLY = get_mean_eff2(country_res, strategy_nums = StrategyNums, reward_input2, r = 0.03)
  
  # get mean eff (undiscounted LY) for the country ------------------------------------
  eff.LY = get_mean_eff2(country_res, strategy_nums = StrategyNums, reward_input2, r = 0)
  
  # PRODUCTIVITY ===============================================================
  prod_ctlg = PRODUCTIVITY
  prod = get_mean_prod2(country_res, strategy_nums = StrategyNums, prod_ctlg, r = 0.03)
  
  # COSTS ======================================================================
  
  # get mean cost (discounted cost) for each country ---------------------------
  costItems = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  # scenario analysis for RNA cost ---------------------------------------------
  TBCOSTS_HC["TEST.RNA"] = 150
  costItems150 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  TBCOSTS_HC["TEST.RNA"] = 60
  costItems60 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  TBCOSTS_HC["TEST.RNA"] = 30
  costItems30 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  TBCOSTS_HC["TEST.RNA"] = 15
  costItems15 = get_mean_cost_itemized(country_res, StrategyNums, TBCOSTS_HC, TBCOSTS_nonHC, COSTS_HC, COSTS_nonHC, r = 0.03)
  
  # collect results into a list -----------------------------------------------
  cea = list(eff.dQALY = eff.dQALY,
             eff.QALY  = eff.QALY,
             eff.dLY   = eff.dLY,
             eff.LY    = eff.LY,
             prod      = prod,
             costItems = costItems,
             costItems150 = costItems150,
             costItems60  = costItems60,
             costItems30  = costItems30,
             costItems15  = costItems15)
  # NMB.No_screening = NMB.No_screening,
  # NMB.IGRA_TB = NMB.IGRA_TB,
  # NMB.IGRA_RNA_TB = NMB.IGRA_RNA_TB,
  # iNMB.IGRA_TB = iNMB.IGRA_TB,
  # iNMB.IGRA_RNA_TB = iNMB.IGRA_RNA_TB)
  
  # save cea result to RDS file -----------------------------------------------
  file.out2 = paste0(ResultsDir_PSA_current,"cea_",name,".RDS")
  saveRDS(cea, file.out2)
  
  print(paste0(file.out2, "...written."))
  print(mem_used())
  #mem_tracker = append(mem_tracker, mem_used())
  gc()
}

# time.start = Sys.time()
# test = get_psa_res(strategies = StrategyNums, name= CountryNames[84])
# time.end = Sys.time()
# time.end - time.start
# print(paste0("Clock time:", time.end - time.start))

# file_list = list.files(path = paste0("ResultsTemp/PSA/paramset",PARAMset_num), pattern = "\\.csv$")
# ctemp = which(!CountryNames %in% str_sub(file_list, 11,13))

# file_list = list.files(path =ResultsDir_PSA_current, pattern = ".RDS")
file_list = list.files(path =ResultsDir_PSA_current, pattern = ".RDS")

ctemp = which(!CountryNames %in% str_sub(file_list, 5,7))
# ctemp = 1:length(CountryNames)


time.start.sim = Sys.time()
set.seed(02215)
# x = mclapply(c(1:length(CountryNames)), function(x) run_model(input = SIM[[x]], strategies = StrategyNums, name = CountryNames[x]), mc.cores = numCores)
# #x = future_lapply(c(1:length(CountryNames)), function(x) get_cea_res(input = SIM[[x]], strategies = StrategyNums, name = CountryNames[x]), future.seed = NULL)
# #x = future_lapply(c(1:length(CountryNames)), function(x) run_model2(strategies = StrategyNums, name = CountryNames[x]), future.seed = NULL)

# x = future_lapply(CountryNames, function(x) get_psa_res(strategies = StrategyNums, name = x), future.seed = TRUE)
# x = lapply(CountryNames, function(x) get_psa_res(strategies = StrategyNums, name = x))
# x = mclapply(CountryNames[ctemp], function(x) get_psa_res3(strategies = StrategyNums, name = x), mc.cores= numCores)]
x = mclapply(CountryNames[ctemp], function(x) get_psa_res4(strategies = StrategyNums, name = x), mc.cores = numCores)
# set.seed(NULL)
# names(x) = CountryNames[ctemp]
time.end.sim = Sys.time()
print(paste0("Clock time:", time.end.sim - time.start.sim))



 


