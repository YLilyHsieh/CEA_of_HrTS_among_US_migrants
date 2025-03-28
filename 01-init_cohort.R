# This file creates the study cohort that runs through the simulation
# This output of this file contains a csv file that is used in the main simulation file 

source("Model/00-libraries.R")
source("Model/00-gen_helpers.R")
if (SA %in%  c("PSA_main", "PSA_RNAgam")){
  source("Model/01-init_helpers.R")
} else if (SA == "PSA_reinf2"){
  source("Model/01-init_helpers_reinf2.R")
} else if (SA == "PSA_noreinf"){
  source("Model/01-init_helpers_noreinf.R")
}

## Define basic characteristics of the initial cohort

source("Data/data.R")
cohort_df = cohort2019_by_coo_age # from data.R
# rm(cohort2019_by_coo_age)
# cohort_df = read.csv("Data/cohort2019_by_coo_age.csv")


# CountryNames = sort(unique(cohort_df$place_of_birth)[which((unique(cohort_df$place_of_birth) %in% c("BGD")))])

# CountryNames = sort(unique(cohort_df$place_of_birth)[which(!(unique(cohort_df$place_of_birth) %in% c("MEX","CHN","IND", "PHL")))])
CountryNames = sort(unique(cohort_df$place_of_birth))
FLAGGED_REGIONS = CountryNames


# Update the parameters included in the PSA
if (PSA == TRUE){
  TB_Symp_to_Dx             = PARAMset[["TB_Symp_to_Dx"]]
  
  IGRA_SPEC_HEALTHY         = PARAMset[["IGRA_SPEC_HEALTHY"]]
  IGRA_SENS_LTBI            = IGRA_SENS_TB = PARAMset[["IGRA_SENS_LTBI"]]
  PROB_FAIL_TO_INIT_LTBI_TX = PARAMset[["PROB_FAIL_TO_INIT_LTBI_TX"]]
  PROB_FAIL_TO_COMP_LTBI_TX = PARAMset[["PROB_FAIL_TO_COMP_LTBI_TX"]]
  PROB_CURED_LTBI_TX        = PARAMset[["PROB_CURED_LTBI_TX"]]
  PROB_FAIL_TO_INIT_INCIPIENT_TX  = PARAMset[["PROB_FAIL_TO_INIT_INCIPIENT_TX"]]
  PROB_FAIL_TO_COMP_INCIPIENT_TX  = PARAMset[["PROB_FAIL_TO_COMP_INCIPIENT_TX"]]
  PROB_CURED_INCIPIENT_TX   = PARAMset[["PROB_CURED_INCIPIENT_TX"]]
  PROB_FAIL_TO_COMP_TB_TX_SURVIGING_TX = PARAMset[["PROB_FAIL_TO_COMP_TB_TX_SURVIGING_TX"]]
  RR_earlyTB_1 = PARAMset[["RR_earlyTB_1"]]
  RR_earlyTB_2 = PARAMset[["RR_earlyTB_2"]]
  
  RNA_SENS1 = RNA_SENS_0_6 = RNA_SENS_6_12 = RNA_SENS_12_18 = RNA_SENS_18_upper = RNA_SENS_TB = PARAMset[["RNA_SENS1"]]
  RNA_SPEC  = PARAMset[["RNA_SPEC"]]
  RNA_SENS_upper_plus = PARAMset[["RNA_SENS_upper_plus"]]
  RNA_TIME_upper = PARAMset[["RNA_TIME_upper"]]
  
  
  TBCOSTS_HC[c(TREAT.TB.GUARANTEE, TREAT.TB_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS)] = PARAMset[["COST_TB_TX_HC"]]
  TBCOSTS_HC[c(TREAT.LTBI_INIT.SUCCESS)] =  PARAMset[["COST_LTBI_TX_HC"]]
  TBCOSTS_HC[c(TREAT.INCIPIENT_INIT.SUCCESS)] = PARAMset[["COST_INCIP_TX_HC"]]
  TBCOSTS_HC[c(TEST.IGRA)] = PARAMset[["COST_IGRA"]]
  TBCOSTS_HC[c(TEST.XRAY)] = PARAMset[["COST_XRAY"]]
  

  TBCOSTS_nonHC[c(TREAT.TB.GUARANTEE, TREAT.TB_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS)] = PARAMset[["COST_TB_TX_nonHC"]]
  TBCOSTS_nonHC[c(TREAT.LTBI_INIT.SUCCESS)] =  PARAMset[["COST_LTBI_TX_nonHC"]]
  TBCOSTS_nonHC[c(TREAT.INCIPIENT_INIT.SUCCESS)] = PARAMset[["COST_INCIP_TX_nonHC"]]

  reward_input[c(LTBIposTBtx)] = PARAMset[["REWARD_POST_TB"]]
  reward_input[c(TBonTBtx, TBpostLTBItxonTBtx, TBpostINCIPtxonTBtx)] = PARAMset[["REWARD_CONTROLLED_TB"]]  
  reward_input[c(TB, TBonLTBItx, TBpostLTBItx,  TBpostTBtx, TBpostLTBItxpostTBtx, TBonINCIPtx, TBpostINCIPtx, TBpostINCIPtxpostTBtx, onTBtxDefault)] = PARAMset[["REWARD_UNCONTROLLED_TB"]]
}
