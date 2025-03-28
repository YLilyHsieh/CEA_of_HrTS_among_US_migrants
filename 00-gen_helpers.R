# This file stores functions that support functions in other helper files. 

# =============================================================================
# Directories  
# =============================================================================
ResultsDir_InitData = paste0(getwd(),"/Data/")
ResultsDir_SummRes  = paste0(getwd(),"/ResultsTemp/summary_results")
ResultsDir_PSA = paste0("~/IncipienTB/PSA_base/")
 

# =============================================================================
# Model settings 
# =============================================================================

# Strategies 
STRATEGY_NAMES_complete = c("Do Nothing", # 0
                            "IGRA_TB",    # 1
                            "IGRA_RNA_TB",# 2
                            "RNA_TB")     # 3


# =============================================================================
# CREATE FUNCTIONS 
# =============================================================================

# Functions to help create new model.matrix for each agent -------------
# the goal is to obtain the risk of developing TB over time
get_agent_tbrisk_dist = function(agentVector){
  
  sim_years = 100 - as.numeric(agentVector["entry_age"]) + 1
  
  Xnew = as.data.frame(matrix(NA, nrow = sim_years, ncol = ncol(reg.out1$model)))
  colnames(Xnew) = colnames(reg.out1$model)
  
  age = as.numeric(agentVector["entry_age"]) + seq(0,(sim_years-1))
  
  
  Xnew$tb_cases = NA
  Xnew[,c(2:14)] = 0 # dummy variables
  Xnew$entry_this_year = c(1, rep(0, sim_years-1))
  Xnew$age_over_90 = ifelse(age > 90, 1, 0) 
  Xnew$yoe_pre_1950 = 0
  Xnew$years_since = seq(0:(sim_years-1))
  Xnew$entry_age = as.numeric(agentVector["entry_age"])
  Xnew$entry_year0 = 32.56782 + (ENTRY_YEAR - 2016)   
  Xnew$`(offset)` = 0
  Xnew$origin = as.character(agentVector["place_of_birth"])
  
  
  sim.risk.annual = predict(reg.out1, Xnew, type = "response") # incidence per 100k population
   
  return(sim.risk.annual)
}

get_agent_tbrisk_dist_psa = function(agentVector){
  
  sim_years = 100 - as.numeric(agentVector["entry_age"]) + 1
  
  Xnew = as.data.frame(matrix(NA, nrow = sim_years, ncol = ncol(reg.out1$model)))
  colnames(Xnew) = colnames(reg.out1$model)
  
  age = as.numeric(agentVector["entry_age"]) + seq(0,(sim_years-1))
  
  
  Xnew$tb_cases = NA
  Xnew[,c(2:14)] = 0 # dummy variables
  Xnew$entry_this_year = c(1, rep(0, sim_years-1))
  Xnew$age_over_90 = ifelse(age > 90, 1, 0) 
  Xnew$yoe_pre_1950 = 0
  Xnew$years_since = seq(0:(sim_years-1))
  Xnew$entry_age = as.numeric(agentVector["entry_age"])
  Xnew$entry_year0 = 32.56782 + (ENTRY_YEAR - 2016)   
  Xnew$`(offset)` = 0
  Xnew$origin = as.character(agentVector["place_of_birth"])
  

  # sim.risk.annual = predict(reg.out1, Xnew, type = "response") # incidence per 100k population

  # prepare LP matrix
  lpm <- predict(fit_list5b[[1]], newdata = Xnew, type = "lpmatrix")
  for (imp in 2:5) {
    lpm <- lpm + predict(fit_list5b[[imp]], newdata = Xnew, type = "lpmatrix")
  }
  comb.lpm0  <- lpm/5

  n.sim <- 1000
  set.seed(1234)
  sim_pars <- MASS::mvrnorm(n.sim, mu = comb.est, Sigma = comb.cov)

  mx <- exp( comb.lpm0 %*% t(sim_pars) )*1.287
   
  return(mx)
}


get_agent_timeToTB = function(agentID, entryAge, sim.risk.annual, smooth){
  
  # smooth, logical: TRUE if we want time to event on a daily scale
  
  if (smooth == FALSE){
    sim.risk.annual = sim.risk.annual*1.287
    func.hazard = -log(1-sim.risk.annual[rownames(sim.risk.annual) == as.character(entryAge),]/1e5)
    
    func.hazard = func.hazard[!is.na(func.hazard)]
    func.hazard[length(func.hazard) + c(1:50)] = func.hazard[length(func.hazard)]
    func.hazard[length(func.hazard) + 1] = 1e8
    func.hazard = c(0, func.hazard)
    
    func.hazard.cum = cumsum(func.hazard[1:length(func.hazard)])
    
    func.surv = c(exp(-func.hazard.cum))
    
    cdf = 1-func.surv # the 0 at the front allows the probability of death at year zero (first year) to be what it is (i.e., it gives it its mass)
    set.seed(agentID)
    timeToTB = min(which(cdf >= runif(1, min = min(cdf), max = max(cdf))))- 1 -1 # first-1 accounts for the 0 added to cdf; second -1 accounts for zero indexing of months since entry
  } else {
    sim.risk.annual = sim.risk.annual*1.287
    func.hazard = -log(1-sim.risk.annual[rownames(sim.risk.annual) == as.character(entryAge),]/1e5)
    
    func.hazard = func.hazard[!is.na(func.hazard)]
    func.hazard[length(func.hazard) + c(1:50)] = func.hazard[length(func.hazard)]
    func.hazard[length(func.hazard) + 1] = 1e8
    func.hazard = c(0, func.hazard)
    
    func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
    func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
    
    cdf = 1-func.surv # the 0 at the front allows the probability of death at year zero (first year) to be what it is (i.e., it gives it its mass)
    set.seed(agentID)
    timeToTB = (min(which(cdf >= runif(1, min = min(cdf), max = max(cdf))))- 1 -1)/365
  }
  
  return(timeToTB)
}

get_agent_timeToTB2 = function(agentID, entryAge, sim.risk.annual, smooth){
  
  # smooth, logical: TRUE if we want time to event on a daily scale
  
  if (smooth == FALSE){
    func.hazard = -log(1-sim.risk.annual/1e5)
    
    func.hazard = func.hazard[!is.na(func.hazard)]
    func.hazard[length(func.hazard) + c(1:50)] = func.hazard[length(func.hazard)]
    func.hazard[length(func.hazard) + 1] = 1e8
    func.hazard = c(0, func.hazard)
    
    func.hazard.cum = cumsum(func.hazard[1:length(func.hazard)])
    
    func.surv = c(exp(-func.hazard.cum))
    
    cdf = 1-func.surv # the 0 at the front allows the probability of death at year zero (first year) to be what it is (i.e., it gives it its mass)
    set.seed(agentID)
    timeToTB = min(which(cdf >= runif(1, min = min(cdf), max = max(cdf))))- 1 -1 # first-1 accounts for the 0 added to cdf; second -1 accounts for zero indexing of months since entry
  } else {
    func.hazard = -log(1-sim.risk.annual/1e5)
    
    func.hazard = func.hazard[!is.na(func.hazard)]
    func.hazard[length(func.hazard) + c(1:50)] = func.hazard[length(func.hazard)]
    func.hazard[length(func.hazard) + 1] = 1e8
    func.hazard = c(0, func.hazard)
    
    func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
    func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
    
    
    cdf = 1-func.surv # the 0 at the front allows the probability of death at year zero (first year) to be what it is (i.e., it gives it its mass)
    set.seed(agentID)
    timeToTB = (min(which(cdf >= runif(1, min = min(cdf), max = max(cdf))))- 1 -1)/365
  }
  
  return(timeToTB)
}

get_agent_timeToTB_psa = function(annual.tb.risk, entryAge, agentID){
  vec = future_apply(annual.tb.risk, 2, function(x) get_agent_timeToTB2(agentID = agentID,entryAge = entryAge, sim.risk.annual = x, smooth = T))
  return(vec)
}

get_cdf_df = function(sim.risk.annual){
  func.hazard = -log(1-sim.risk.annual/1e5)
  
  func.hazard = func.hazard[!is.na(func.hazard)]
  func.hazard[length(func.hazard) + c(1:50)] = func.hazard[length(func.hazard)]
  func.hazard[length(func.hazard) + 1] = 1e8
  func.hazard = c(0, func.hazard)
  
  func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
  func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
  
  
  cdf = 1-func.surv # the 0 at the front allows the probability of death at year zero (first year) to be what it is (i.e., it gives it its mass)
  
  return(cdf)
}



get_tttb_df = function(samp, cdf){
  time_to_tb = (sapply(1:length(samp), function(x) min(which(cdf >= samp[x]))) -1 -1 )/365
  return(time_to_tb)
}

get_agent_timeToTB_psa_efficient = function(annual.tb.risk, strata_popsize, seed){
  
  cdf_df = apply(annual.tb.risk, 2, function(x) get_cdf_df(x))
  
  set.seed(seed)
  samp = runif(strata_popsize, min = 0, max = 1)
  
  timeToTB_df = apply(cdf_df, 2, function(x) get_tttb_df(samp = samp, cdf = x))
  
  return(timeToTB_df) # this is an n by 1000 matrix with time to tb for each age in this age group, across 1000 annual TB risk draw
}

# get_agent_timeToTB_efficient = function(entryAge, strata_popsize, sim.risk.annual, smooth, seed){ # sample time to TB for an age group from a COO
#   
#   func.hazard = -log(1-sim.risk.annual/1e5)
#   
#   func.hazard = func.hazard[!is.na(func.hazard)]
#   func.hazard[length(func.hazard) + c(1:50)] = func.hazard[length(func.hazard)]
#   func.hazard[length(func.hazard) + 1] = 1e8
#   func.hazard = c(0, func.hazard)
#   
#   func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
#   func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
#   
#   
#   cdf = 1-func.surv # the 0 at the front allows the probability of death at year zero (first year) to be what it is (i.e., it gives it its mass)
#  
#   set.seed(seed)
#   samp = runif(strata_popsize, min = min(cdf), max = max(cdf))
#   
#   timeToTB = (sapply(1:length(samp), function(x) min(which(cdf >= samp[x]))) -1 -1 )/365
#   
#   timeToTB = (min(which(cdf >= runif(1, min = min(cdf), max = max(cdf))))- 1 -1)/365
#   
#   return(timeToTB)
# }


# Demographis =================================================================

# Mortality rates 

# import 2017 life table
MORT = read.csv("Data/2017lt_foreign_born.csv")


get_timeToDeath = function(agentID, entry_age, current_age, life_table = MORT, update, earlyRR_1, earlyRR_2, delayedRR, smooth){
  
  # for the entry cohort, current_age = entry_age
  # update, logical, FALSE when initializing cohort, TRUE when simulation
  # smooth, logical; TRUE if daily, FALSE if annual
  # earlyRR_1 is the rate ratio in the first 6 months 
  # earlyRR_2 is the rate ratio in months 7-12 (or during treatment)
  # delayedRR is the rate ratio after first year (or after treatment)
  
  if (update == FALSE){
    
    earlyRR_1 = earlyRR_2 = delayedRR = NULL
    current_age = entry_age
    
    if (floor(entry_age) > 99){
      timeToDeath = 0
    } else {
      
      if (smooth == FALSE){
        func.hazard = life_table[which(life_table$age == floor(entry_age)) : which(life_table$age == 100), "mu_asr"]
        func.hazard[length(func.hazard)] = 1e8 # die for sure when they are 100 years old
        func.hazard = c(0, func.hazard)
        
        func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
        
        cdf = 1 - func.surv
        set.seed(agentID)
        timeToDeath = min(which(cdf >= runif(1, min = min(cdf), max = max(cdf)))) - 1 - 1 
        #first -1 accounts for the additional 0 added to cdf
        #second -1 accounts for zero indexing of years since entry
      } else {
        func.hazard = life_table[which(life_table$age == floor(entry_age)) : which(life_table$age == 100), "mu_asr"]
        func.hazard[length(func.hazard)] = 1e8 # die for sure when they are 100 years old
        # func.hazard = c(0, func.hazard)
        func.hazard = c(0,func.hazard)
        
        func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
        
        func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
        
        cdf = 1- func.surv
        set.seed(agentID)
        timeToDeath = (min(which(cdf >= runif(1, min = min(cdf), max = max(cdf)))) - 1 -1)/365
      }
    }
  } else {
    
    if (floor(entry_age) > 99){
      timeToDeath = 0
    } else {
      
      if (smooth == FALSE){
        func.hazard = life_table[which(life_table$age == floor(entry_age)) : which(life_table$age == 100), "mu_asr"]
        func.hazard[length(func.hazard)] = 1e8 # die for sure when they are 100 years old
        func.hazard = c(0, func.hazard)
        
        func.hazard[floor(current_age)] = func.hazard[floor(current_age)] * earlyRR_2 # need to be cautious about this
        func.hazard[(floor(current_age) + 1):(floor(current_age) + 7)] = func.hazard[(floor(current_age) + 1):(floor(current_age) + 7)] * delayedRR
        func.hazard[(floor(current_age) + 8):length(func.hazard)] = func.hazard[(floor(current_age) + 8):length(func.hazard)] 
        
        func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
        
        cdf = 1- func.surv
        cdf = cdf[!is.na(cdf)] # drop NA generated when a > 92 years old
        
        set.seed(agentID)
        timeToDeath = min(which(cdf >= runif(1, min = min(cdf), max = max(cdf)))) - 1 - 1 
        
      } else {
        func.hazard = life_table[which(life_table$age == floor(entry_age)) : which(life_table$age == 100), "mu_asr"]
        func.hazard[length(func.hazard)] = 1e8 # die for sure when they are 100 years old
        func.hazard = c(0, func.hazard)
        
        func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
        
        func.hazard.age = c(NA,seq(floor(entry_age-1), 99, 1/365)+1)
        
        a = which.min(abs(current_age - func.hazard.age))
        
        func.hazard[a:(a + 182)] = func.hazard[a:(a + 182)] * earlyRR_1
        func.hazard[(a + 183):(a + 364)] = func.hazard[(a + 183):(a + 364)] * earlyRR_2
        func.hazard[(a + 364*1 + 1):(a + 364*7)] = func.hazard[(a + 364*1 + 1):(a + 364*7)] * delayedRR
        func.hazard[(a + 364*7 + 1):length(func.hazard)] = func.hazard[(a + 364*7 + 1):length(func.hazard)]
        
       func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
        
        func.cdf = 1- func.surv
        func.cdf = func.cdf[!is.na(func.cdf)] # drop NA generated when a > 92 years old
        set.seed(agentID)
        timeToDeath = (min(which(func.cdf >= runif(1, min = min(func.cdf), max = max(func.cdf)))) - 1 -1)/365
        
        if(is.infinite(timeToDeath)){stop("Check current_age.")}
      }
    }
    
  }
  
  return(timeToDeath)
}

get_cdf_df_death = function(entry_age, life_table = MORT){
  
  current_age = entry_age
  
  if (floor(entry_age) > 99){
    timeToDeath = 0
  } else {

      func.hazard = life_table[which(life_table$age == floor(entry_age)) : which(life_table$age == 100), "mu_asr"]
      func.hazard[length(func.hazard)] = 1e8 # die for sure when they are 100 years old
      func.hazard = c(0,func.hazard)
      
      func.hazard = c((approx(func.hazard[1:length(func.hazard) - 1], method = "linear", n = 365*(length(func.hazard) - 2) + 1)$y)/365, 1e8)
      
      func.surv = c(exp(-cumsum(func.hazard[1:length(func.hazard)])))
      
      cdf = 1- func.surv
  }
  return(cdf)
}

get_ttdeath_df = function(samp, cdf){
  time_to_death = (future_sapply(1:length(samp), function(x) min(which(cdf >= samp[x]))) -1 -1 )/365
  return(time_to_death)
}

get_agent_timeToDeath_efficient = function(entry_age, strata_popsize, seed){
  
  cdf_df = get_cdf_df_death(entry_age, life_table = MORT)
  
  set.seed(seed)
  samp = runif(strata_popsize, min = 0, max = 1)
  
  timeToDeath_df = get_ttdeath_df(samp = samp, cdf = cdf_df)
  
  return(timeToDeath_df) # this is an n by 1000 matrix with time to tb for each age in this age group, across 1000 annual TB risk draw
}

# =============================================================================
# CREATE VARIABLES 
# =============================================================================

# A list of health states -----------------------------------------------------
HEALTHY  = "HEALTHY"
LTBI     = "LTBI"
TB       = "TB"
DEATH    = "DEATH"

LTBIonLTBItx    = "LTBIonLTBItx"
LTBIonINCIPtx   = "LTBIonINCIPtx"
LTBIpostLTBItx  = "LTBIpostLTBItx"
LTBIpostINCIPtx = "LTBIpostINCIPtx"

LTBIpostTBtx    = "LTBIpostTBtx"


HEALTHYonLTBItx    = "HEALTHYonLTBItx"
HEALTHYonINCIPtx   = "HEALTHYonINCIPtx"
HEALTHYpostLTBItx  = "HEALTHYpostLTBItx"
HEALTHYpostINCIPtx = "HEALTHYpostINCIPtx"
HEALTHYpostTBtx    = "HEALTHYpostTBtx"


TBonLTBItx    = "TBonLTBItx"
TBonINCIPtx   = "TBonINCIPtx"
TBonTBtx      = "TBonTBtx"
TBpostLTBItx  = "TBpostLTBItx"
TBpostINCIPtx = "TBpostINCIPtx"
TBpostTBtx    = "TBpostTBtx"
TBpostLTBItxpostTBtx = "TBpostLTBItxpostTBtx"
TBpostLTBItxonTBtx   = "TBpostLTBItxonTBtx"
TBpostINCIPtxonTBtx = "TBpostINCIPtxonTBtx"
TBpostINCIPtxpostTBtx = "TBpostINCIPtxpostTBtx"
LTBIposTBtx = "LTBIposTBtx"

onTBtxDefault = "onTBtxDefault"

# Mortality rates -------------------------------------------------------------

# RR           = 1.00
# RR_earlyTB   = 1.78 #1.78 #1.14 #7.29
# RR_delayedTB = 1.00 #1.78 #1.14 #1.78

RR_earlyTB_1 = 2.2658 # 1.223 # 2.2658#1.1329 #1.16 # 1.18 # calibrated to % died at diagnosis
RR_earlyTB_2 = 8.842 #9.8583 #4.92 # calibrated % died on treatment
RR_delayedTB =1.78 # Post-TB


# Test characteristics --------------------------------------------------------

IGRA_SPEC_HEALTHY = 0.98 #0.985 
IGRA_SENS_LTBI    = 0.89 #0.789  
IGRA_SENS_TB      = 0.89 #0.789
# IGRA_SENS_TB      = 1


XRAY_SPEC_HEALTHY = 1.00 # Assumption 
XRAY_SPEC_LTBI    = 1.00 # Assumption 
XRAY_SENS_TB      = 1.00 # Assumption 

RNA_SPEC          = 0.90 # WHO # 0.75 # 0.8  # Zak_2016_Lancet, Table 1
RNA_SENS1         = 0.90

RNA_SENS_0_6    = RNA_SENS1 # 0.90 # 0.712 
RNA_SENS_6_12   = RNA_SENS1 # 0.90 # 0.629
RNA_SENS_12_18  = RNA_SENS1 # 0.90 # 0.477
RNA_SENS_18_upper  = RNA_SENS1 # 0.90 # 0.291
RNA_SENS_upper_plus = 0.10 # 0.25 # 0.054

RNA_SENS_TB     = RNA_SENS_0_6

RNA_SENS_SAfactor = 1 # unless otherwise specified in sim_cohort.R when SA is performed

RNA_TIME_upper = 720/365

# Treatment care cascade ------------------------------------------------------

# Timing of post-arrival screening
POSTARRIVAL_SCREENING = 1/12 

# Lag between TB symptom development and TB diagnosis and treatment
TB_Symp_to_Dx = 0.5 # assume 6 months # 24/365
# Lag between TB diagnosis (from screening) to treatment
TB_Dx_to_Tx   = 3/365
# Lag between getting final test result to treatment
Test_to_Tx = 0
# Time between TB treatment initiation and death on treatment
TB_Tx_to_newTx = 0.5 # 6 months

if (POSTARRIVAL_SCREENING >TB_Symp_to_Dx) stop("Check assumptions on the timing of post-arrival screening and time-to-diagnosis for active TB cases.")


# A list of checkpoints 
NOT.ENTER.ALGO = "NOT.ENTER.ALGO"
ENTER.ALGO = "ENTER.ALGO"

MONITOR = "MONITOR" # corresponding timeToNextCheckpoint = NA

TEST.IGRA = "TEST.IGRA"
TEST.RNA  = "TEST.RNA"  


TEST.XRAY = "TEST.XRAY"
TEST.RNA_XRAY = "TEST.RNA_XRAY"
TESTING   = c(TEST.IGRA, TEST.RNA, TEST.XRAY, TEST.RNA_XRAY)

TREAT.TB           = "TREAT.TB"
TREAT.TB.REP       = "TREAT.TB.REP"
TREAT.TB.GUARANTEE = "TREAT.TB.GUARANTEE"
TREAT.INCIPIENT    = "TREAT.INCIPIENT"
TREAT.LTBI         = "TREAT.LTBI"

TREAT.TB_INIT.FAIL        = "TREAT.TB_INIT.FAIL"
TREAT.TB.REP_INIT.FAIL    = "TREAT.TB.REP_INIT.FAIL"
TREAT.INCIPIENT_INIT.FAIL = "TREAT.INCIPIENT_INIT.FAIL"
TREAT.LTBI_INIT.FAIL      = "TREAT.LTBI_INIT.FAIL"

TREAT.TB_INIT.SUCCESS        = "TREAT.TB_INIT.SUCCESS"
TREAT.TB.REP_INIT.SUCCESS    = "TREAT.TB.REP_INIT.SUCCESS"
TREAT.INCIPIENT_INIT.SUCCESS = "TREAT.INCIPIENT_INIT.SUCCESS"
TREAT.LTBI_INIT.SUCCESS      = "TREAT.LTBI_INIT.SUCCESS"

TREAT.TB_DEAD           = "TREAT.TB_DEAD"
TREAT.TB.REP_DEAD       = "TREAT.TB.REP_DEAD"
TREAT.TB.GUARANTEE_DEAD = "TREAT.TB.GUARANTEE_DEAD"
TREAT.INCIPIENT_DEAD    = "TREAT.INCIPIENT_DEAD"
TREAT.LTBI_DEAD         = "TREAT.LTBI_DEAD"

TREAT.TB_DEFAULT        = "TREAT.TB_DEFAULT"
TREAT.TB.REP_DEFAULT    = "TREAT.TB.REP_DEFAULT"
TREAT.TB.GUARANTEE_DEFAULT = "TREAT.TB.GUARANTEE_DEFAULT"
TREAT.INCIPIENT_DEFAULT = "TREAT.INCIPIENT_DEFAULT"
TREAT.LTBI_DEFAULT      = "TREAT.LTBI_DEFAULT"

TREAT.TB_COMPLETE           = "TREAT.TB_COMPLETE"
TREAT.TB.REP_COMPLETE       = "TREAT.TB.REP_COMPLETE"
TREAT.TB.GUARANTEE_COMPLETE = "TREAT.TB.GUARANTEE_COMPLETE"
TREAT.INCIPIENT_COMPLETE    = "TREAT.INCIPIENT_COMPLETE"
TREAT.LTBI_COMPLETE         = "TREAT.LTBI_COMPLETE"


TREAT.LTBI_COMPLETE_NOTCURED = "TREAT.LTBI_COMPLETE_NOTCURED"
TREAT.LTBI_COMPLETE_CURED = "TREAT.LTBI_COMPLETE_CURED"
TREAT.INCIPIENT_COMPLETE_NOTCURED = "TREAT.INCIPIENT_COMPLETE_NOTCURED"
TREAT.INCIPIENT_COMPLETE_CURED = "TREAT.INCIPIENT_COMPLETE_CURED"

DEAD.AT.TB.DX = "DEAD.AT.TB.DX"

TREAT.DEAD      = c(TREAT.TB_DEAD,     TREAT.TB.REP_DEAD,     TREAT.LTBI_DEAD,     TREAT.INCIPIENT_DEAD,     TREAT.TB.GUARANTEE_DEAD)
TREAT.DEFAULT   = c(TREAT.TB_DEFAULT,  TREAT.TB.REP_DEFAULT,  TREAT.LTBI_DEFAULT,  TREAT.INCIPIENT_DEFAULT,  TREAT.TB.GUARANTEE_DEFAULT)
TREAT.COMPLETE  = c(TREAT.TB_COMPLETE, TREAT.TB.REP_COMPLETE, TREAT.TB.GUARANTEE_COMPLETE,
                    TREAT.INCIPIENT_COMPLETE_NOTCURED,
                    TREAT.INCIPIENT_COMPLETE_CURED,
                    TREAT.LTBI_COMPLETE_NOTCURED,
                    TREAT.LTBI_COMPLETE_CURED)
TREAT.END       = c(TREAT.DEFAULT, TREAT.COMPLETE, TREAT.DEAD)

TREATMENT              = c(TREAT.TB,              TREAT.LTBI,              TREAT.INCIPIENT,              TREAT.TB.REP)
TREATMENT.INIT.FAIL    = c(TREAT.TB_INIT.FAIL,    TREAT.LTBI_INIT.FAIL,    TREAT.INCIPIENT_INIT.FAIL,    TREAT.TB.REP_INIT.FAIL)
TREATMENT.INIT.SUCCESS = c(TREAT.TB_INIT.SUCCESS, TREAT.LTBI_INIT.SUCCESS, TREAT.INCIPIENT_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS, TREAT.TB.GUARANTEE)

DEAD    = "DEAD"
ALL.CAUSE.DEATH   = c(DEAD, DEAD.AT.TB.DX, TREAT.DEAD)
END     = "END"

#PROB_DEAD_AT_TB_DX = 0.0151

PROB_DEAD_ON_TX_LTBI_TX = 0        
#PROB_DEAD_ON_TX_TB_TX_VIA_SYMPTOM   = 0.0555
#PROB_DEAD_ON_TX_TB_TX_VIA_SCREENING = 0.03
PROB_DEAD_ON_TX_INCIPIENT_TX = 0

PROB_FAIL_TO_INIT_LTBI_TX = 0.238     
PROB_FAIL_TO_INIT_TB_TX   = 0.000         
PROB_FAIL_TO_INIT_INCIPIENT_TX = 0.119 

PROB_FAIL_TO_COMP_LTBI_TX = 0.097       # Sterling_2011_NEJM, Table 3
#PROB_FAIL_TO_COMP_TB_TX   = 0.126       # CDC Reported TB 2019, Table 51
PROB_FAIL_TO_COMP_TB_TX_SURVIGING_TX = 0.084
PROB_FAIL_TO_COMP_INCIPIENT_TX = 0.084 # 0.126  # Assumption

PROB_CURED_TB_TX = 1
PROB_CURED_LTBI_TX = 0.64
PROB_CURED_INCIPIENT_TX = 0.64


TTE_TB_TX_DEFAULT = 3/52
TTE_TB_DEAD_ON_TX = 4/52

# CEA parameters --------------------------------------------------------------

# r = 0.03 # annual discount rate 
WTP = 150000

# reward ------------------------------------------------------------------------------------------
reward_input = NULL # rep(NA, length(unique(event_dat.aug$HealthStatus)))
reward_input[c(DEAD)] = 0
reward_input[c(HEALTHY, HEALTHYpostLTBItx, HEALTHYpostINCIPtx, 
               LTBI, LTBIpostLTBItx, LTBIpostINCIPtx)] = 1
reward_input[c(LTBIonLTBItx, LTBIonINCIPtx, HEALTHYonLTBItx, HEALTHYonINCIPtx)] = 0.999
reward_input[c(LTBIpostTBtx)] = 0.99
reward_input[c(TBonTBtx, TBpostLTBItxonTBtx, TBpostINCIPtxonTBtx)] = 0.859   # This LTBIposTBtx is treated TB, not LTBI accidentally put on TB treatment 
reward_input[c(TB, TBonLTBItx, TBpostLTBItx,  TBpostTBtx, TBpostLTBItxpostTBtx, TBonINCIPtx, TBpostINCIPtx, TBpostINCIPtxpostTBtx, onTBtxDefault)] = 0.75


# cost --------------------------------------------------------------------------------------------
# ALL IN 2021 DOLLARS

# TB-related healthcare expenditures
TBCOSTS_HC = NULL #rep(NA, length(unique(event_dat.aug$Event)))
TBCOSTS_HC[c(DEAD, "Flag", NOT.ENTER.ALGO,ENTER.ALGO, TB)] = 0
TBCOSTS_HC[c(TREAT.TB.GUARANTEE, TREAT.TB_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS)] = 20993.27
TBCOSTS_HC[c(TREAT.TB.GUARANTEE_COMPLETE, TREAT.TB.GUARANTEE_DEFAULT, TREAT.TB_DEFAULT, TREAT.TB.REP_DEFAULT, TREAT.LTBI_COMPLETE_CURED,
             TREAT.LTBI_COMPLETE_NOTCURED, TREAT.LTBI_DEFAULT, TREAT.LTBI_INIT.FAIL, TREAT.TB_COMPLETE, TREAT.TB.REP_COMPLETE, TREAT.TB.REP, TREAT.INCIPIENT_INIT.FAIL, TREAT.INCIPIENT_DEFAULT,
             TREAT.INCIPIENT_COMPLETE_CURED, TREAT.INCIPIENT_COMPLETE_NOTCURED,
             TREAT.TB.REP_COMPLETE)]= 0
TBCOSTS_HC[c(TREAT.LTBI_INIT.SUCCESS)] =  427.95
TBCOSTS_HC[c(TREAT.INCIPIENT_INIT.SUCCESS)] = 513.92
TBCOSTS_HC[c(TEST.IGRA)] = 64.38  
TBCOSTS_HC[c(TEST.RNA)]  = 300
TBCOSTS_HC[c(TEST.XRAY)] = 34.49  

# TB-related non-healthcare expenditures
TBCOSTS_nonHC = NULL #rep(NA, length(unique(event_dat.aug$Event)))
TBCOSTS_nonHC[c(DEAD, "Flag", NOT.ENTER.ALGO,ENTER.ALGO, TB)] = 0
TBCOSTS_nonHC[c(TREAT.TB.GUARANTEE, TREAT.TB_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS)] = 4622.24
TBCOSTS_nonHC[c(TREAT.TB.GUARANTEE_COMPLETE, TREAT.TB.GUARANTEE_DEFAULT, TREAT.TB_DEFAULT, TREAT.TB.REP_DEFAULT, TREAT.LTBI_COMPLETE_CURED,
             TREAT.LTBI_COMPLETE_NOTCURED, TREAT.LTBI_DEFAULT, TREAT.LTBI_INIT.FAIL, TREAT.TB_COMPLETE, TREAT.TB.REP_COMPLETE, TREAT.TB.REP, TREAT.INCIPIENT_INIT.FAIL, TREAT.INCIPIENT_DEFAULT,
             TREAT.INCIPIENT_COMPLETE_CURED, TREAT.INCIPIENT_COMPLETE_NOTCURED,
             TREAT.TB.REP_COMPLETE)]= 0
TBCOSTS_nonHC[c(TREAT.LTBI_INIT.SUCCESS)] =  112.18
TBCOSTS_nonHC[c(TREAT.INCIPIENT_INIT.SUCCESS)] = 149.57
TBCOSTS_nonHC[c(TEST.IGRA)] = 0 
TBCOSTS_nonHC[c(TEST.RNA)]  = 0
TBCOSTS_nonHC[c(TEST.XRAY)] = 0

# TB-unrelated healthcare expenditures
COSTS_HC     = read.csv("Data/econ_eval/healthcare_costs_by_age_2021usd.csv")

# non-healthcare related expenditures
COSTS_nonHC  = read.csv("Data/econ_eval/nonhealthcare_costs_by_age_2021usd.csv")

# productivity
PRODUCTIVITY = read.csv("Data/econ_eval/productivity_by_age_2021usd.csv")
