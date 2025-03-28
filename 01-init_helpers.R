# This file creates functions and variable classes that:
# (1) create the simulation study cohort 
# (2) instantiate arrivals (to doctor's appointment)

# =============================================================================
# CREATE CLASSES 
# =============================================================================

# Create AgentBio class -------------------------------------------------------
agentBio = list(eventHistory = NA , nextHealthState = NA, 
                timeToTB = NA, timeToTB.bl = NA, timeToDeath = NA, timeToDeath.bl = NA, timeToDeath.tb = NA,
                nextCheckpoint = NA, timeToNextCheckpoint = NA,
                nextEvent = NA, timeToNextEvent = NA, nextEventTime = NA,
                randomNumbers = NA)
class(agentBio) = "AgentBio"

# Create Event class ----------------------------------------------------------
currentEvent = data.frame(eventTime = NA, eventName = NA, healthStatusAtEvent = NA, ageAtEvent = NA)


# =============================================================================
# CREATE FUNCTIONS 
# =============================================================================

# A function to create the study cohort ---------------------------------------
create_cohort = function(CountryName, cohort_raw, fraction, PSA){
  
  ## get iso numeric code
  iso.num = unique(cohort_raw$iso3num)
  
  ## Expand out the rows of the cohort dataframe to the population size
  cohort_mod = cohort_raw %>% slice(rep(1:n(), pop_size_est*fraction))
  
  ## Background mortality -----------------------------------------------------
  # time_to_death
  cohort_mod$time_to_death = sapply(1:nrow(cohort_mod), function(x) get_timeToDeath(x + iso.num, entry_age = cohort_mod$entry_age[x], life_table = MORT, update = FALSE, smooth = TRUE))
  # temp1 = read.csv(paste0("Data/init_time_to_death/",CountryName,".csv"))
  # cohort_mod$time_to_death = temp1$time_to_death
  # rm(temp1)
  # age_at_death
  cohort_mod$age_at_death = cohort_mod$time_to_death + cohort_mod$entry_age
  
  ## Initial TB status --------------------------------------------------------
  # determine time_to_TBdx among all agents 
  if (PSA == FALSE){
    if (fraction < 1) {
      temp2 = read.csv(paste0(ResultsDir_InitData,"init_time_to_tb/init_time_to_TB_",CountryName,"_trunc.csv"))
    } else {
      temp2 = read.csv(paste0(ResultsDir_InitData,"init_time_to_tb/init_time_to_TB_",CountryName,".csv"))
    }
    
    cohort_mod$time_to_TBdx = temp2$time_to_TBdx
    rm(temp2)
    
  } else {
    temp2 = read.csv(paste0("Data/init_time_to_tb_psa_calib_09-30/",CountryName,".csv"))
    temp2 = temp2[,PARAMset_num]
    
    cohort_mod$time_to_TBdx = temp2
    rm(temp2)
  }
  
  
  # determine time_to_TB among all agents 
  cohort_mod$time_to_TB = cohort_mod$time_to_TBdx - TB_Symp_to_Dx
  cohort_mod$time_to_TB[cohort_mod$time_to_TB <= 0] = 0 # so that evaluate_competing_risk() doesn't return with a negative time to event value
  
  cohort_mod$age_at_TB = cohort_mod$entry_age + cohort_mod$time_to_TB
  
  cohort_mod$time_to_death_tb = sapply(1:nrow(cohort_mod), function(x) get_timeToDeath(x + iso.num, entry_age = cohort_mod$entry_age[x], current_age = cohort_mod$age_at_TB[x], 
                                                                                       earlyRR_1 = RR_earlyTB_1, earlyRR_2 = RR_earlyTB_2, delayedRR = RR_delayedTB, life_table = MORT, update = TRUE, smooth = TRUE)) 
  
  
  # determine who will ever develop & be diagnosed with TB disease 
  cohort_mod$lifetimeTB = ifelse((cohort_mod$time_to_TB) < cohort_mod$time_to_death, 1, 0) 
  # this has to be less than, without equal sign, based on how time to TB and time to death functions are set up
  
  # gather the number of people who will ever develop TB across in each age group (required for next step)
  cohort_raw$tb_caseCount = NA
  for (i in 1:nrow(cohort_raw)){
    cohort_raw$tb_caseCount[i] = sum(cohort_mod$entry_age == cohort_raw$entry_age[i] & cohort_mod$lifetimeTB == 1)
  }
  
  # find out the number of LTBI we should expect in each age group
  n.ltbi_expected = sapply(1:nrow(cohort_raw), function(x) rbinom(cohort_raw$pop_size_est[x]*fraction, 1, cohort_raw$ltbi_calib[x]/100))
  
  n.ltbi_expected = sapply(n.ltbi_expected, sum)
  
  # find out the number of LTBI we have identified in each age group based on time_to_TB data
  cohort_mod$LTBI = NA
  cohort_mod$LTBI[cohort_mod$lifetimeTB == 1] = 1
  
  n.ltbi_confirmed = cohort_mod %>% group_by(entry_age) %>% summarise(n = sum(LTBI, na.rm = T)) 
  
  # find out the number of LTBI we still need in each age group so that the LTBI prevalence matches the estimates based on NHANES data
  n.ltbi_left      = n.ltbi_expected - c(n.ltbi_confirmed$n)
  # if n.ltbi_left < 0, recode to 0 # don't need more people to enter with LTBI
  if (any(n.ltbi_left < 0) ) {write.csv(n.ltbi_left, file = paste0(ResultsDir_PSA,CountryName,"_",PARAMset_num,".csv"))}
  n.ltbi_left      = ifelse(n.ltbi_left < 0, 0, n.ltbi_left)
  
  
  # fill in the rest of LTBI 
  for (i in 1:nrow(cohort_raw)){
    
    total.na = sum(is.na(cohort_mod$LTBI) & cohort_mod$entry_age == cohort_raw$entry_age[i])
    recode_to_1 = n.ltbi_left[i]
    recode_to_0 = total.na - recode_to_1
    
    cohort_mod$LTBI[is.na(cohort_mod$LTBI) & cohort_mod$entry_age == cohort_raw$entry_age[i]] = c(rep(1, recode_to_1), rep(0, recode_to_0))
  }
  
  # determine TB status at the time of screening
  cohort_mod = cohort_mod %>% 
    mutate(., TB_status = case_when(
      (LTBI == 1 & time_to_TB > 0) ~ "LTBI",
      (LTBI == 1 & time_to_TB <= 0) ~ "TB",
      TRUE ~ "HEALTHY"
    ))
  
  
  ## Reformat and attach random numbers 
  cohort = data.frame(id                  = seq(1:nrow(cohort_mod)),
                      entry_age           = cohort_mod$entry_age,
                      year_of_entry       = ENTRY_YEAR,                   # defined in source file data.R
                      place_of_birth      = cohort_mod$place_of_birth,
                      TB_status           = cohort_mod$TB_status,
                      time_to_TB          = cohort_mod$time_to_TB, 
                      time_to_death       = cohort_mod$time_to_death,
                      time_to_death_tb    = cohort_mod$time_to_death_tb,
                      rnum.igra           = runif(nrow(cohort_mod)),
                      rnum.xray           = runif(nrow(cohort_mod)),
                      rnum.rna            = runif(nrow(cohort_mod)),
                      rnum.fail_init_tx   = runif(nrow(cohort_mod)),
                      rnum.fail_comp_tx   = runif(nrow(cohort_mod)),
                      rnum.cured          = runif(nrow(cohort_mod)))
  
  cohort$prop.preventable.tb = 0.5^(cohort$time_to_TB/20)
  cohort$prop.preventable.tb = ifelse(cohort$prop.preventable.tb > 1, 1, cohort$prop.preventable.tb)
  cohort$tb.reinfection  = sapply(cohort$prop.preventable.tb, function(x) rbinom(1, 1, prob = (1-x)))
  cohort$rnum.cured.mod = ifelse(cohort$tb.reinfection == 1, 1, cohort$rnum.cured) # won't be cured. sim_helpers.R is.cured()
  
  return(cohort)
}


# A function to determine whether an agent failed to be screened --------------
is.fail_to_be_screened = function(prob_fail_to_be_screened){
  
  is_fail_to_screen = ifelse(runif(1) < prob_fail_to_be_screened, TRUE, FALSE)
  
  return(is_fail_to_screen)
}

# A function to create each agent and flag for screening ----------------------
init_agent = function(initial_attributes, strategy){
  # initial attributes: a row entry of the output from create_cohort()
  
  
  # attach an empty event to the event list
  current_event = currentEvent
  # populate the event list
  current_event$eventTime = 0
  current_event$eventName = "Flag"
  current_event$healthStatusAtEvent = initial_attributes$TB_status
  current_event$ageAtEvent = initial_attributes$entry_age
  
  
  # attach an empty agent bio to an agent
  agentbio = agentBio
  # populate the agent bio
  agentbio$eventHistory = current_event
  
  
  agentbio$nextHealthState    = initial_attributes$TB_status
  agentbio$timeToTB           = initial_attributes$time_to_TB
  agentbio$timeToTB.bl        = initial_attributes$time_to_TB
  agentbio$timeToDeath        = initial_attributes$time_to_death
  agentbio$timeToDeath.bl     = initial_attributes$time_to_death
  agentbio$timeToDeath.tb     = initial_attributes$time_to_death_tb
  agentbio$randomNumbers      = c(igra = initial_attributes$rnum.igra,
                                  rna  = initial_attributes$rnum.rna,
                                  xray = initial_attributes$rnum.xray,
                                  fail_init_tx = initial_attributes$rnum.fail_init_tx,
                                  fail_comp_tx = initial_attributes$rnum.fail_comp_tx,
                                  cured        = initial_attributes$rnum.cured,
                                  cured.mod = initial_attributes$rnum.cured.mod)
  
  agentbio$timeToDeath.tb  = ifelse(agentbio$timeToDeath.tb < agentbio$timeToTB.bl, 
                                    agentbio$timeToTB.bl + 1/365,
                                    agentbio$timeToDeath.tb) 
  
  if(initial_attributes$place_of_birth %in% FLAGGED_REGIONS){
    # if agent is from high TB burden countries,
    # roll a dice to determine whether they fail to go onto the screening process 
    fail_to_screen = FALSE # use is.fail_to_be_screened() if we want to allow imperfect screening 
    
    if (strategy == 0 | fail_to_screen) {   
      agentbio$nextEvent = NOT.ENTER.ALGO
      agentbio$timeToNextEvent = 0  # no time lapse between agent going from init_agent to sim_agent
      agentbio$nextEventTime = current_event$eventTime + agentbio$timeToNextEvent
      
      agentbio$nextCheckpoint = MONITOR
      agentbio$timeToNextCheckpoint = NA
      
    } else {
      agentbio$nextEvent = ENTER.ALGO
      agentbio$timeToNextEvent = 0 # no time lapse between agent going from init_agent to sim_agent
      agentbio$nextEventTime = current_event$eventTime + agentbio$timeToNextEvent
      
      agentbio$timeToNextCheckpoint = POSTARRIVAL_SCREENING # defined in gen_helpers.R
      
      # nextCheckpoint depends on which strategy we are evaluating
      if (strategy == 1){
        agentbio$nextCheckpoint = TEST.IGRA
      } else if (strategy == 2){
        agentbio$nextCheckpoint = TEST.IGRA
      } else if (strategy == 3){
        agentbio$nextCheckpoint = TEST.RNA
      }
    } 
    
  } else {
    
    # if agent is not from high TB burden countries   
    agentbio$nextEvent = NOT.ENTER.ALGO
    agentbio$timeToNextEvent = 0 # no time lapse between agent going from init_agent to sim_agent
    agentbio$nextEventTime = current_event$eventTime + agentbio$timeToNextEvent
    
    agentbio$nextCheckpoint = MONITOR
    agentbio$timeToNextCheckpoint = NA
    
  }
  
  return(agentbio)
}


# Create a function to initialize cohort for each country ----------------------

init_country_cohort = function(CountryName, callSavedFiles, PSA){
  
  if (PSA == FALSE){
    
    callSavedFiles = NULL
    
    if (CountryName %in% c("MEX","CHN","IND")){
      fraction = 0.1
    } else if (CountryName == "PHL") {
      fraction = 0.5
    } else {
      fraction = 1
    }
    
    #if (callSavedFiles == FALSE){
    
    # Subset data to the country-of-origin of interest
    cohort_df_subset = NULL
    cohort_df_subset = cohort_df[cohort_df$place_of_birth == CountryName,]
    
    # Create a dataframe of agents for that country-of-origin
    cohort = create_cohort(CountryName = CountryName, cohort_raw = cohort_df_subset, fraction = fraction, PSA = PSA)
    
    
    # Initialize the cohort for simulation
    sim0 = list()
    sim0 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 0)) 
    
    sim1 = list()
    sim1 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 1))   
    
    sim2 = list()
    sim2 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 2))
    
    sim3 = list()
    sim3 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 3))
    
    
    print(paste0("sim for ",CountryName," <-- Done."))
    
    return(list(sim0 = sim0, sim1 = sim1, sim2 = sim2, sim3 = sim3))
    # }
    
  } else {
    
    callSavedFiles = NULL
    
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
    
    # Initialize the cohort for simulation
    sim0 = list()
    sim0 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 0)) 
    
    sim1 = list()
    sim1 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 1))   
    
    sim2 = list()
    sim2 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 2))
    
    sim3 = list()
    sim3 = lapply(1:nrow(cohort), function(x) init_agent(cohort[x,], strategy = 3))
    
    print(paste0("sim for ",CountryName," <-- Done."))
    
    return(list(sim0 = sim0, sim1 = sim1, sim2 = sim2, sim3 = sim3))
  }
  
}

