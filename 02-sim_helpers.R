# This file creates functions that simulates each agent through the screening and treatment processes

# =============================================================================
# CREATE SUPPORTING FUNCTIONS 
# =============================================================================

# A function to evaluate competing risks --------------------------------------
evaluate_competing_risks = function(timeToTB, timeToDeath, timeToNextCheckpoint, checkpointName){
  
  if (is.na(timeToTB) && is.na(timeToDeath) && is.na(timeToNextCheckpoint)){
    stop("All time-to-event values are null.")}
  if (any(c(timeToTB, timeToDeath, timeToNextCheckpoint) < 0, na.rm = T)){
    stop("Time-to-event values cannot be negative.")}
  
  # if (!is.na(timeToDeath) & !is.na(timeToTB)){
  #   if (abs(timeToDeath-timeToTB) < 0.0001){
  #     timeToTB = timeToTB + 1    # If time to TB and time to death is the same, then we allow agents to die 
  #   } else {
  #     timeToDeath = timeToDeath
  #     timeToTB    = timeToTB
  #   }
  # }
  # 
  # 
  # if (!is.na(timeToDeath) & !is.na(timeToNextCheckpoint)){
  #   if (abs(timeToDeath-timeToNextCheckpoint) < 0.0001){
  #     timeToDeath = timeToDeath + 1
  #   } else {
  #     timeToDeath = timeToDeath
  #     timeToNextCheckpoint    = timeToNextCheckpoint
  #   }
  # }
  
  
  dict = NULL
  dict = data.frame(Event = c("TB", "DEAD", checkpointName),
                    TimeToEvent = c(timeToTB, timeToDeath, timeToNextCheckpoint))
  
  next_event = dict[which.min(dict$TimeToEvent),"Event"]
  time_to_event = dict[which.min(dict$TimeToEvent), "TimeToEvent"]
  
  return(list(next_event = next_event, time_to_event = time_to_event))
}

# A function to append event to eventHistory-----------------------------------
append_event = function(last_agentbio_history, current_event){
  
  updated.history = rbind(last_agentbio_history, current_event) 
  
  return(updated.history)
}

# A function to determine IGRA result -----------------------------------------
get_igra_result = function(health_status, igra_spec_healthy, igra_sens_ltbi, igra_sens_tb, agent_rnum){
  
  if (health_status == HEALTHY){
    is_pos = if_else(agent_rnum > igra_spec_healthy, TRUE, FALSE)
  } else if (health_status == LTBI){
    is_pos = if_else(agent_rnum < igra_sens_ltbi, TRUE, FALSE)
  } else if (health_status == TB){
    is_pos = if_else(agent_rnum < igra_sens_tb, TRUE, FALSE)
  }
  
  return(is_pos)
}

# A function to determine XRAY result -----------------------------------------
get_xray_result = function(health_status, xray_spec_healthy, xray_spec_ltbi, xray_sens_tb, agent_rnum){
  
  if (health_status == HEALTHY){
    is_pos = if_else(agent_rnum > xray_spec_healthy, TRUE, FALSE)
  } else if (health_status == LTBI){
    is_pos = if_else(agent_rnum > xray_spec_ltbi, TRUE, FALSE)
  } else if (health_status == TB){
    is_pos = if_else(agent_rnum < xray_sens_tb, TRUE, FALSE)
  }
  
  return(is_pos)
}

# A function to determine RNA result ------------------------------------------
get_rna_result = function(health_status, time_to_tb, agent_rnum){
  
  if (health_status == HEALTHY){
    spec = RNA_SPEC                               
    is_pos = if_else(agent_rnum > spec, TRUE, FALSE)
  }
  
  if (health_status == TB){
    sens = RNA_SENS_TB * RNA_SENS_SAfactor
    is_pos = if_else(agent_rnum < sens, TRUE, FALSE)
  } 
  
  if (health_status == LTBI){
    
    if (time_to_tb <= 180/365){
      sens = RNA_SENS_0_6 * RNA_SENS_SAfactor
    } else if (time_to_tb > 180/365 && time_to_tb <= 360/365){
      sens = RNA_SENS_6_12 * RNA_SENS_SAfactor
    } else if (time_to_tb > 360/365 && time_to_tb <= 540/365){
      sens = RNA_SENS_12_18 * RNA_SENS_SAfactor
    } else if (time_to_tb > 540/365 && time_to_tb <= RNA_TIME_upper){
      sens = RNA_SENS_18_upper * RNA_SENS_SAfactor
    } else {
      sens = RNA_SENS_upper_plus * RNA_SENS_SAfactor
    }
    
    is_pos = if_else(agent_rnum < sens, TRUE, FALSE)
    
  } 
  
  return(is_pos)
}

# # A function to determine whether a TB patient is dead at diagnosis -----------
# is.dead_at_diagnosis = function(agent_rnum){
#   
#   is_dead_at_dx = if_else(agent_rnum < PROB_DEAD_AT_TB_DX, TRUE, FALSE)
#   
#   return(is_dead_at_dx)
# }

# A function to determine whether treatment is successfully initiated ---------
is.fail_to_init_treatment = function(treatment_type, agent_rnum){
  
  if(treatment_type %in% c(TREAT.TB, TREAT.TB.REP)){
    prob_fail_to_init_treatment = PROB_FAIL_TO_INIT_TB_TX
  } else if (treatment_type == TREAT.LTBI){
    prob_fail_to_init_treatment = PROB_FAIL_TO_INIT_LTBI_TX
  } else if (treatment_type == TREAT.INCIPIENT){
    prob_fail_to_init_treatment = PROB_FAIL_TO_INIT_INCIPIENT_TX
  }
  
  is_fail_to_init_tx = if_else(agent_rnum < prob_fail_to_init_treatment, TRUE, FALSE)
  
  return(is_fail_to_init_tx)
}

# A function to determine treatment outcome -----------------------------
# get.treatment_outcome = function(treatment_type, agent_rnum){
#   
#   if(treatment_type %in% c(TREAT.TB.GUARANTEE)){
#     
#     if ( agent_rnum < PROB_DEAD_ON_TX_TB_TX_VIA_SYMPTOM ){ 
#       tx_outcome = TREAT.TB.GUARANTEE_DEAD
#     } else if (agent_rnum > (1-PROB_FAIL_TO_COMP_TB_TX)){
#       tx_outcome = TREAT.TB.GUARANTEE_DEFAULT
#     } else {
#       tx_outcome = TREAT.TB.GUARANTEE_COMPLETE
#     }
#     
#   } else if (treatment_type == TREAT.TB_INIT.SUCCESS){
#     
#     if ( agent_rnum < PROB_DEAD_ON_TX_TB_TX_VIA_SCREENING ){ 
#       tx_outcome = TREAT.TB_DEAD
#     } else if (agent_rnum > (1-PROB_FAIL_TO_COMP_TB_TX)){
#       tx_outcome = TREAT.TB_DEFAULT
#     } else {
#       tx_outcome = TREAT.TB_COMPLETE
#     }
#     
#   } else if (treatment_type == TREAT.TB.REP_INIT.SUCCESS){
#     
#     if ( agent_rnum < PROB_DEAD_ON_TX_TB_TX_VIA_SYMPTOM ){ 
#       tx_outcome = TREAT.TB.REP_DEAD
#     } else if (agent_rnum > (1-PROB_FAIL_TO_COMP_TB_TX)){
#       tx_outcome = TREAT.TB.REP_DEFAULT
#     } else {
#       tx_outcome = TREAT.TB.REP_COMPLETE
#     }
#     
#     
#   } else if (treatment_type == TREAT.LTBI_INIT.SUCCESS){
#     
#     if ( agent_rnum < PROB_DEAD_ON_TX_LTBI_TX ){ 
#       tx_outcome = TREAT.LTBI_DEAD
#     } else if (agent_rnum > (1-PROB_FAIL_TO_COMP_LTBI_TX)){
#       tx_outcome = TREAT.LTBI_DEFAULT
#     } else {
#       tx_outcome = TREAT.LTBI_DEAD
#     }
#     
#   } else if (treatment_type == TREAT.INCIPIENT_INIT.SUCCESS){
#     
#     if ( agent_rnum < PROB_DEAD_ON_TX_INCIPIENT_TX ){ 
#       tx_outcome = TREAT.INCIPIENT_DEAD
#     } else if (agent_rnum > (1-PROB_FAIL_TO_COMP_INCIPIENT_TX)){
#       tx_outcome = TREAT.INCIPIENT_DEFAULT
#     } else {
#       tx_outcome = TREAT.INCIPIENT_COMPLETE
#     }
#     
#   }
#   
#   return(tx_outcome)
# }

# # A function to determine whether a person died on TB treatment ---------
# is.die_on_treatment = function(treatment_type, agent_rnum){
#   
#   if(treatment_type %in% c(TREAT.TB.REP_INIT.SUCCESS, TREAT.TB.GUARANTEE)){
#     prob_dead_on_treatment = PROB_DEAD_ON_TX_TB_TX_VIA_SYMPTOM
#   } else if (treatment_type == TREAT.TB_INIT.SUCCESS){
#     prob_dead_on_treatment = PROB_DEAD_ON_TX_TB_TX_VIA_SCREENING
#   } else if (treatment_type == TREAT.LTBI_INIT.SUCCESS){
#     prob_dead_on_treatment = PROB_DEAD_ON_TX_LTBI_TX
#   } else if (treatment_type == TREAT.INCIPIENT_INIT.SUCCESS){
#     prob_dead_on_treatment = PROB_DEAD_ON_TX_INCIPIENT_TX
#   }
#   
#   is_dead_on_tx = if_else(agent_rnum < prob_dead_on_treatment, TRUE, FALSE)
#   
#   return(is_dead_on_tx)
# }
# 
# A function to determine whether treatment is successfully completed ---------
is.fail_to_comp_treatment = function(treatment_type, agent_rnum){
  
  
  if(treatment_type %in% c(TREAT.TB_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS, TREAT.TB.GUARANTEE)){
    prob_fail_to_comp_treatment = PROB_FAIL_TO_COMP_TB_TX_SURVIGING_TX   
  } else if (treatment_type == TREAT.LTBI_INIT.SUCCESS){
    prob_fail_to_comp_treatment = PROB_FAIL_TO_COMP_LTBI_TX
  } else if (treatment_type == TREAT.INCIPIENT_INIT.SUCCESS){
    prob_fail_to_comp_treatment = PROB_FAIL_TO_COMP_INCIPIENT_TX
  }
  
  is_fail_to_comp_tx = if_else(agent_rnum < prob_fail_to_comp_treatment, TRUE, FALSE)
  
  return(is_fail_to_comp_tx)
}

# A function to sample from time to treatment default dist --------------------
get_time_to_tx_default = function(treatment_type){  #//// Placeholder /////# 
  
  # set.seed(seed)
  # if(treatment_type == TREAT.TB_INIT.SUCCESS | treatment_type == TREAT.TB.GUARANTEE){
  #   tte = 3/52
  # } else if (treatment_type == TREAT.LTBI_INIT.SUCCESS){
  #   tte = 3/52
  # } else if (treatment_type == TREAT.INCIPIENT_INIT.SUCCESS){
  #   tte = 3/52
  # }
  
  tte = TTE_TB_TX_DEFAULT
  
  return(tte)
}

# A function to sample from time to death on treatment dist -------------------
get_time_to_dead_on_tx = function(treatment_type){
  
  # if(treatment_type == TREAT.TB_INIT.SUCCESS | treatment_type == TREAT.TB.GUARANTEE){
  #   tte = TTE_TB_DEAD_ON_TX
  # } else if (treatment_type == TREAT.LTBI_INIT.SUCCESS){
  #   tte = NA
  # } else if (treatment_type == TREAT.INCIPIENT_INIT.SUCCESS){
  #   tte = NA
  # }
  
  tte = TTE_TB_DEAD_ON_TX
  
  return(tte)
}

# A function to get time to treatment completion ------------------------------
get_time_to_tx_comp = function(treatment_type){
  
  if(treatment_type %in% c(TREAT.TB_INIT.SUCCESS, TREAT.TB.REP_INIT.SUCCESS, TREAT.TB.GUARANTEE)){
    tte = 26/52
  } else if (treatment_type == TREAT.LTBI_INIT.SUCCESS) {
    tte = 12/52
  } else if (treatment_type == TREAT.INCIPIENT_INIT.SUCCESS){
    tte = 18/52
  }
  
  return(tte)
}

# A function to determine whether an agent is cured after completion of tx ----
is.cured = function(treatment_type, agent_rnum){
  
  if(treatment_type == TREAT.LTBI_COMPLETE){
    prob_cured = PROB_CURED_LTBI_TX
  } else if (treatment_type == TREAT.INCIPIENT_COMPLETE){
    prob_cured = PROB_CURED_INCIPIENT_TX
  }
  
  is.cured = if_else(agent_rnum < prob_cured, TRUE, FALSE)
  
  return(is.cured)
}


# =============================================================================
# CORE SIMULATION FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# A function to simulate agents over their trajectories 
# -----------------------------------------------------------------------------

sim_agent = function(agentbio_on_file, strategy, agentID){
  
  
  if (!(strategy %in% seq(0,3,1))){
    stop("Strategy not valid. Check function arguments.")}
  
  ## copy agentbio from last event
  agentbio = agentbio_on_file   
  
  
  if (agentbio$nextCheckpoint == MONITOR | strategy == 0){ 
    
    # Step 1. Record current event 
    current_event = currentEvent
    current_event$eventTime           = agentbio$nextEventTime
    current_event$eventName           = agentbio$nextEvent
    current_event$healthStatusAtEvent = agentbio$nextHealthState
    current_event$ageAtEvent          = agentbio$eventHistory[1,4] + current_event$eventTime
    
    # Step 2. Update agent bio on time-to-health-events based on current event
    if(current_event$eventName == TB){
      agentbio$timeToTB    = NA # already has TB
      agentbio$timeToDeath = agentbio$timeToDeath.tb - current_event$eventTime
    } else {
      agentbio$timeToTB    = agentbio$timeToTB - agentbio$timeToNextEvent     # NB: here the TimeToNextEvent is the time elapsed between last event and the current event
      agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
    }
    
    # Step 3. Determine what the next event is by evaluating competing risks 
    evaluation = evaluate_competing_risks(agentbio$timeToTB,
                                          agentbio$timeToDeath,
                                          agentbio$timeToNextCheckpoint,
                                          agentbio$nextCheckpoint)
    
    # Step 4. Update agent bio on nextEvent, timeToNextEvent, and nextEventTime
    agentbio$nextEvent       = evaluation[["next_event"]]
    agentbio$timeToNextEvent = evaluation[["time_to_event"]]
    agentbio$nextEventTime   = current_event$eventTime + agentbio$timeToNextEvent
    
    # Step 5. Update agent bio on nextHealthState based on nextEvent
    agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB,  ALL.CAUSE.DEATH),
                                       agentbio$nextEvent,
                                       current_event$healthStatusAtEvent)
    
    # Step 6. Update agent bio on nextCheckPoint and timeToNextCheckpoint based on nextEvent
    if (agentbio$nextEvent == TB){
      agentbio$nextCheckpoint       = TREAT.TB.GUARANTEE
      agentbio$timeToNextCheckpoint = TB_Symp_to_Dx
    } else if (agentbio$nextEvent %in% ALL.CAUSE.DEATH) {
      agentbio$nextCheckpoint       = END
      agentbio$timeToNextCheckpoint = NA 
    } else {
      agentbio$nextCheckpoint       = MONITOR
      agentbio$timeToNextCheckpoint = NA
    }
    
    # Step 7. Update eventHistory by appending current_event to eventHistory before exiting this checkpoint
    agentbio$eventHistory = append_event(agentbio$eventHistory, current_event)
  } 
  
  
  if (agentbio$nextCheckpoint %in% TESTING){
    
    # Step 1. Record current event
    current_event = currentEvent
    current_event$eventTime           = agentbio$nextEventTime
    current_event$eventName           = agentbio$nextEvent
    current_event$healthStatusAtEvent = agentbio$nextHealthState
    current_event$ageAtEvent          = agentbio$eventHistory[1,4] + current_event$eventTime
    
    # Step 2. Update agent bio on time-to-health-events based on current event
    if(current_event$eventName == TB){
      agentbio$timeToTB    = NA # already has TB
      agentbio$timeToDeath = agentbio$timeToDeath.tb - current_event$eventTime
    } else {
      agentbio$timeToTB    = agentbio$timeToTB - agentbio$timeToNextEvent     # Here the TimeToNextEvent is the time elapsed between last event and the current event
      agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
    }
    
    # Step 3. Determine what the next event is by evaluating competing risks 
    evaluation = evaluate_competing_risks(agentbio$timeToTB,
                                          agentbio$timeToDeath,
                                          agentbio$timeToNextCheckpoint,
                                          agentbio$nextCheckpoint)
    
    # Step 4. Update agent bio on nextEvent, timeToNextEvent, and nextEventTime
    agentbio$nextEvent       = evaluation[["next_event"]]
    agentbio$timeToNextEvent = evaluation[["time_to_event"]]
    agentbio$nextEventTime   = current_event$eventTime + agentbio$timeToNextEvent
    
    # Step 5. Update agent bio on nextHealthState for nextEvent
    agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),
                                       agentbio$nextEvent,
                                       current_event$healthStatusAtEvent)
    
    # Step 6. Update agent bio on nextCheckPoint and timeToNextCheckpoint based on nextEvent
    # if (agentbio$nextEvent == TB){
    # 
    #   if (agentbio$nextEventTime < POSTARRIVAL_SCREENING){
    #     agentbio$nextCheckpoint       = agentbio$nextCheckpoint        # TEST.RNA or TEST.IGRA
    #     agentbio$timeToNextCheckpoint = POSTARRIVAL_SCREENING - agentbio$nextEventTime # agentbio$timeToNextCheckpoint  # This should be the same as POSTARRIVAL_SCREENING
    #   }
    #   else {
    #     agentbio$nextCheckpoint       = TREAT.TB.GUARANTEE
    #     agentbio$timeToNextCheckpoint = TB_Symp_to_Dx
    #   }
    #    
    if (agentbio$nextEvent == TB & agentbio$nextEventTime < POSTARRIVAL_SCREENING){
      
      agentbio$nextCheckpoint       = agentbio$nextCheckpoint        # TEST.RNA or TEST.IGRA
      agentbio$timeToNextCheckpoint = POSTARRIVAL_SCREENING - agentbio$nextEventTime # agentbio$timeToNextCheckpoint  # This should be the same as POSTARRIVAL_SCREENING
      
    } else if (agentbio$nextEvent %in% ALL.CAUSE.DEATH){
      agentbio$nextCheckpoint = END
      agentbio$timeToNextCheckpoint = NA  
    } else if (agentbio$nextEvent == MONITOR){
      agentbio$timeToNextCheckpoint = NA
    } else {
      
      # if nextEvent is part of a testing algorithm
      if (strategy == 1){
        
        if (agentbio$nextEvent == TEST.IGRA){
          
          IGRA_pos = get_igra_result(current_event$healthStatusAtEvent,
                                     IGRA_SPEC_HEALTHY, IGRA_SENS_LTBI, IGRA_SENS_TB, agentbio$randomNumbers[["igra"]])
          
          if(IGRA_pos){
            agentbio$nextCheckpoint = TEST.XRAY  
            agentbio$timeToNextCheckpoint = 1/52 
          } else {
            if(current_event$healthStatusAtEvent == TB){
              agentbio$nextCheckpoint = TREAT.TB.GUARANTEE
              agentbio$timeToNextCheckpoint = ifelse ((agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime) < 0, 
                                                      0, agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime)
            } else {
              agentbio$nextCheckpoint = MONITOR
              agentbio$timeToNextCheckpoint = NA
            }
            # agentbio$nextCheckpoint = MONITOR
            # agentbio$timeToNextCheckpoint = NA
          }
          
        } else if (agentbio$nextEvent == TEST.XRAY){
          
          XRAY_pos = get_xray_result(current_event$healthStatusAtEvent, 
                                     XRAY_SPEC_HEALTHY, XRAY_SPEC_LTBI, XRAY_SENS_TB, agentbio$randomNumbers[["xray"]])
          
          if(XRAY_pos){
            agentbio$nextCheckpoint = TREAT.TB
            agentbio$timeToNextCheckpoint = TB_Dx_to_Tx  
          } else {
            agentbio$nextCheckpoint = TREAT.LTBI
            agentbio$timeToNextCheckpoint = Test_to_Tx
          }
        }
        
      } else if (strategy == 2) {
        
        if(agentbio$nextEvent == TEST.IGRA){
          
          IGRA_pos = get_igra_result(current_event$healthStatusAtEvent,
                                     IGRA_SPEC_HEALTHY, IGRA_SENS_LTBI, IGRA_SENS_TB, agentbio$randomNumbers[["igra"]])
          
          if(IGRA_pos){
            agentbio$nextCheckpoint = TEST.RNA  
            agentbio$timeToNextCheckpoint = 1/52  
          } else {
            if(current_event$healthStatusAtEvent == TB){
              agentbio$nextCheckpoint = TREAT.TB.GUARANTEE
              agentbio$timeToNextCheckpoint = ifelse ((agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime) < 0, 
                                                      0, agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime)
            } else {
              agentbio$nextCheckpoint = MONITOR
              agentbio$timeToNextCheckpoint = NA
            }
            # agentbio$nextCheckpoint = MONITOR
            # agentbio$timeToNextCheckpoint = NA
          }
          
        } else if (agentbio$nextEvent == TEST.RNA){
          
          RNA_pos = get_rna_result(current_event$healthStatusAtEvent, 
                                   agentbio$timeToTB - agentbio$timeToNextEvent, agentbio$randomNumbers[["rna"]])
          
          if(RNA_pos){
            agentbio$nextCheckpoint = TEST.XRAY
            agentbio$timeToNextCheckpoint = 1/52  
          } else{
            if(current_event$healthStatusAtEvent == TB){
              agentbio$nextCheckpoint = TREAT.TB.GUARANTEE
              agentbio$timeToNextCheckpoint = ifelse ((agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime) < 0, 
                                                      0, agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime)
            } else {
              agentbio$nextCheckpoint = MONITOR
              agentbio$timeToNextCheckpoint = NA
            }
            # agentbio$nextCheckpoint = MONITOR
            # agentbio$timeToNextCheckpoint = NA 
          }
          
        } else if (agentbio$nextEvent == TEST.XRAY){
          
          XRAY_pos = get_xray_result(current_event$healthStatusAtEvent, 
                                     XRAY_SPEC_HEALTHY, XRAY_SPEC_LTBI, XRAY_SENS_TB, agentbio$randomNumbers[["xray"]])
          
          if(XRAY_pos){
            agentbio$nextCheckpoint = TREAT.TB
            agentbio$timeToNextCheckpoint =  TB_Dx_to_Tx 
          } else {
            agentbio$nextCheckpoint = TREAT.INCIPIENT
            agentbio$timeToNextCheckpoint = Test_to_Tx 
          }
        }
      } else if (strategy == 3){
        
        if (agentbio$nextEvent == TEST.RNA){
          
          RNA_pos = get_rna_result(current_event$healthStatusAtEvent,
                                   agentbio$timeToTB - agentbio$timeToNextEvent, agentbio$randomNumbers[["rna"]])
          
          if(RNA_pos){
            agentbio$nextCheckpoint = TEST.XRAY  
            agentbio$timeToNextCheckpoint = 1/52 
          } else {
            if(current_event$healthStatusAtEvent == TB){
              agentbio$nextCheckpoint = TREAT.TB.GUARANTEE
              agentbio$timeToNextCheckpoint = ifelse ((agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime) < 0, 
                                                      0, agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime)
            } else {
              agentbio$nextCheckpoint = MONITOR
              agentbio$timeToNextCheckpoint = NA
            }
            # agentbio$nextCheckpoint = MONITOR
            # agentbio$timeToNextCheckpoint = NA
          }
          
        } else if (agentbio$nextEvent == TEST.XRAY){
          
          XRAY_pos = get_xray_result(current_event$healthStatusAtEvent, 
                                     XRAY_SPEC_HEALTHY, XRAY_SPEC_LTBI, XRAY_SENS_TB, agentbio$randomNumbers[["xray"]])
          
          if(XRAY_pos){
            agentbio$nextCheckpoint = TREAT.TB
            agentbio$timeToNextCheckpoint = TB_Dx_to_Tx  
          } else {
            agentbio$nextCheckpoint = TREAT.INCIPIENT
            agentbio$timeToNextCheckpoint = Test_to_Tx
          }
        } 
      }
      # Step 7. update eventHistory by appending current_event to eventHistory before exiting this checkpoint
      agentbio$eventHistory = append_event(agentbio$eventHistory, current_event)
    }
  }
  
  
  if (agentbio$nextCheckpoint %in% TREATMENT){
    # Determine whether there's failure to initiate treatment
    fail_to_treat = is.fail_to_init_treatment(treatment_type = agentbio$nextCheckpoint, 
                                              agentbio$randomNumbers[["fail_init_tx"]])
    
    if (fail_to_treat){
      agentbio$nextCheckpoint = paste0(agentbio$nextCheckpoint,"_INIT.FAIL")
    } else {
      agentbio$nextCheckpoint = paste0(agentbio$nextCheckpoint,"_INIT.SUCCESS")
    }
  } 
  
  
  if (agentbio$nextCheckpoint %in% c(TREATMENT.INIT.FAIL)){
    # Step 1. Record current event
    current_event = currentEvent
    current_event$eventTime           = agentbio$nextEventTime
    current_event$eventName           = agentbio$nextEvent
    current_event$healthStatusAtEvent = agentbio$nextHealthState
    current_event$ageAtEvent          = agentbio$eventHistory[1,4] + current_event$eventTime
    
    # Step 2. Update agent bio on time-to-health-events based on current event
    if(current_event$eventName == TB){
      agentbio$timeToTB    = NA # already has TB
      agentbio$timeToDeath = agentbio$timeToDeath.tb - current_event$eventTime
    } else {
      agentbio$timeToTB    = agentbio$timeToTB - agentbio$timeToNextEvent     # NB: here the TimeToNextEvent is the time elapsed between last event and the current event                
      agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
    }
    
    # Step 3. Determine what the next event is by evaluating competing risks 
    evaluation = evaluate_competing_risks(agentbio$timeToTB,
                                          agentbio$timeToDeath,
                                          agentbio$timeToNextCheckpoint,
                                          agentbio$nextCheckpoint)
    
    # Step 4. Update agent bio on nextEvent, timeToNextEvent, and nextEventTime
    agentbio$nextEvent       = evaluation[["next_event"]]
    agentbio$timeToNextEvent = evaluation[["time_to_event"]]
    agentbio$nextEventTime   = current_event$eventTime + agentbio$timeToNextEvent
    
    # Step 5. Update agent bio on nextHealthState based on nextEvent
    agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),
                                       agentbio$nextEvent,
                                       current_event$healthStatusAtEvent)
    
    # Step 6. Update agent bio on nextCheckPoint and timeToNextCheckpoint based on nextEvent
    if (agentbio$nextEvent == TB){
      agentbio$nextCheckpoint       = TREAT.TB.GUARANTEE
      agentbio$timeToNextCheckpoint = TB_Symp_to_Dx
    } else if (agentbio$nextEvent %in% c(ALL.CAUSE.DEATH)){
      agentbio$nextCheckpoint       = END
      agentbio$timeToNextCheckpoint = NA 
    } else if (agentbio$nextEvent %in% TREATMENT.INIT.FAIL){
      agentbio$nextCheckpoint       = MONITOR
      agentbio$timeToNextCheckpoint = NA
    }
    
    # Step 7. update eventHistory by appending current_event to eventHistory before exiting this checkpoint
    agentbio$eventHistory = append_event(agentbio$eventHistory, current_event)
  }
  
  
  if (agentbio$nextCheckpoint %in% TREATMENT.INIT.SUCCESS){
    # Step 1. Record current event
    current_event = currentEvent
    current_event$eventTime           = agentbio$nextEventTime 
    current_event$eventName           = agentbio$nextEvent
    current_event$healthStatusAtEvent = agentbio$nextHealthState
    current_event$ageAtEvent          = agentbio$eventHistory[1,4] + current_event$eventTime
    
    # Step 2. Update agent bio on time-to-health-events based on current event
    if(current_event$eventName == TB){
      agentbio$timeToTB    = NA # already has TB
      agentbio$timeToDeath = agentbio$timeToDeath.tb - current_event$eventTime
    } else {
      agentbio$timeToTB    = agentbio$timeToTB - agentbio$timeToNextEvent     # NB: here the TimeToNextEvent is the time elapsed between last event and the current event
      agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
    }
    
    # Step 3. Determine what the next event is by evaluating competing risks 
    evaluation = evaluate_competing_risks(agentbio$timeToTB,
                                          agentbio$timeToDeath,
                                          agentbio$timeToNextCheckpoint,
                                          agentbio$nextCheckpoint)
    
    # Step 4. Update agent bio on nextEvent, timeToNextEvent, and nextEventTime
    agentbio$nextEvent       = evaluation[["next_event"]]            # if dead here, then death prior to treatment initiation
    agentbio$timeToNextEvent = evaluation[["time_to_event"]]
    agentbio$nextEventTime   = current_event$eventTime + agentbio$timeToNextEvent
    
    # Step 5. Update agent bio on nextHealthState based on nextEvent
    if (agentbio$nextCheckpoint == TREAT.LTBI_INIT.SUCCESS){
      agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),
                                         agentbio$nextEvent,
                                         paste0(current_event$healthStatusAtEvent, "onLTBItx"))
    } else if (agentbio$nextCheckpoint == TREAT.INCIPIENT_INIT.SUCCESS) {
      agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),
                                         agentbio$nextEvent,
                                         paste0(current_event$healthStatusAtEvent, "onINCIPtx"))
    } else if (agentbio$nextCheckpoint == TREAT.TB.REP_INIT.SUCCESS) {
      agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),
                                         agentbio$nextEvent,
                                         paste0("TBonTBtx"))
    } else {
      agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),
                                         agentbio$nextEvent,
                                         paste0(current_event$healthStatusAtEvent, "onTBtx"))
    } 
    
    # Step 6. Update agent bio on nextCheckPoint and timeToNextCheckpoint based on nextEvent
    if (agentbio$nextEvent == TB){
      agentbio$nextCheckpoint       = TREAT.TB.GUARANTEE
      agentbio$timeToNextCheckpoint = TB_Symp_to_Dx 
    } else if (agentbio$nextEvent == DEAD){
      agentbio$nextCheckpoint       = END
      agentbio$timeToNextCheckpoint = NA 
    } else if (agentbio$nextEvent %in% TREATMENT.INIT.SUCCESS){
      
      
      # conditional on surviving, determine if agent will successfully complete treatment
      treatment_default = is.fail_to_comp_treatment(agentbio$nextEvent, agentbio$randomNumbers[["fail_comp_tx"]])
      
      if (treatment_default){ # conditional on surviving treatment, did not finish the full course of treatment
        
        # update nextCheckpoint
        if (agentbio$nextEvent == TREAT.TB_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.TB_DEFAULT
        } else if (agentbio$nextEvent == TREAT.TB.REP_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.TB.REP_DEFAULT
        } else if (agentbio$nextEvent == TREAT.TB.GUARANTEE){
          agentbio$nextCheckpoint = TREAT.TB.GUARANTEE_DEFAULT
        } else if (agentbio$nextEvent == TREAT.LTBI_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.LTBI_DEFAULT
        } else if (agentbio$nextEvent == TREAT.INCIPIENT_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.INCIPIENT_DEFAULT
        }
        # determine timeToNextCheckpoint
        agentbio$timeToNextCheckpoint = get_time_to_tx_default(agentbio$nextEvent)
        
      } else {  # conditional on surviving treatment, completed the full course of treatment
        
        # update nextCheckpoint
        if (agentbio$nextEvent == TREAT.TB_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.TB_COMPLETE
        } else if (agentbio$nextEvent == TREAT.TB.REP_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.TB.REP_COMPLETE
        } else if (agentbio$nextEvent == TREAT.TB.GUARANTEE){
          agentbio$nextCheckpoint = TREAT.TB.GUARANTEE_COMPLETE
        } else if (agentbio$nextEvent == TREAT.LTBI_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.LTBI_COMPLETE
          cured = is.cured(agentbio$nextCheckpoint, agentbio$randomNumbers[["cured.mod"]])  
          if (cured) {
            agentbio$nextCheckpoint = paste0(agentbio$nextCheckpoint,"_CURED")      # rename nextCheckpoint to record whether the agent is cured
          } else {
            agentbio$nextCheckpoint = paste0(agentbio$nextCheckpoint,"_NOTCURED")  # rename nextCheckpoint to record whether the agent is not cured
          }
        } else if (agentbio$nextEvent == TREAT.INCIPIENT_INIT.SUCCESS){
          agentbio$nextCheckpoint = TREAT.INCIPIENT_COMPLETE
          cured = is.cured(agentbio$nextCheckpoint, agentbio$randomNumbers[["cured.mod"]])  
          if (cured) {
            agentbio$nextCheckpoint = paste0(agentbio$nextCheckpoint,"_CURED")      # rename nextCheckpoint to record whether the agent is cured
          } else {
            agentbio$nextCheckpoint = paste0(agentbio$nextCheckpoint,"_NOTCURED")  # rename nextCheckpoint to record whether the agent is not cured
          }
        }
        # determine timeToNextCheckpoint
        agentbio$timeToNextCheckpoint = get_time_to_tx_comp(agentbio$nextEvent)
      }
    }
    #}
    
    
    # Step 7. update eventHistory by appending current_event to eventHistory before exiting this checkpoint
    agentbio$eventHistory = append_event(agentbio$eventHistory, current_event)
  }
  
  
  if (agentbio$nextCheckpoint %in% TREAT.END){
    
    # Step 1. Record current event
    current_event = currentEvent
    current_event$eventTime           = agentbio$nextEventTime
    current_event$eventName           = agentbio$nextEvent
    current_event$healthStatusAtEvent = agentbio$nextHealthState
    current_event$ageAtEvent          = agentbio$eventHistory[1,4] + current_event$eventTime
    
    # Step 2. Update agent bio on time-to-health-events based on current event AND effect of treatment outcomes
    if(current_event$eventName == TB){
      agentbio$timeToTB    = NA # already has TB
      agentbio$timeToDeath = agentbio$timeToDeath.tb - current_event$eventTime
      
    } else {
      if (agentbio$nextCheckpoint %in% TREAT.DEFAULT){
        agentbio$timeToTB    = agentbio$timeToTB - agentbio$timeToNextEvent 
        agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
        
      } else if (agentbio$nextCheckpoint %in% c(TREAT.LTBI_COMPLETE_CURED, TREAT.INCIPIENT_COMPLETE_CURED)){
        agentbio$timeToTB    = NA # won't progress to TB
        agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
        
      } else if (agentbio$nextCheckpoint %in% c(TREAT.LTBI_COMPLETE_NOTCURED, TREAT.INCIPIENT_COMPLETE_NOTCURED)){
        agentbio$timeToTB    = agentbio$timeToTB - agentbio$timeToNextEvent 
        agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
        
      } else {
        agentbio$timeToTB    = NA  # had TB already, cured of TB 
        agentbio$timeToDeath = agentbio$timeToDeath - agentbio$timeToNextEvent
      }
    }
    
    # Step 3. Determine what the next event is by evaluating competing risks 
    # Step 4. Update agent bio on nextEvent, timeToNextEvent, and nextEventTime
    if (agentbio$nextCheckpoint %in% c(TREAT.TB.GUARANTEE_DEFAULT, TREAT.TB.GUARANTEE_COMPLETE) & agentbio$timeToDeath <= 0.5){  # Those who were set to die on treatment would die 
      agentbio$nextEvent       = DEAD
      agentbio$timeToNextEvent = agentbio$timeToDeath
      agentbio$nextEventTime   = current_event$eventTime + agentbio$timeToNextEvent
    } else {
      evaluation = evaluate_competing_risks(agentbio$timeToTB,
                                            agentbio$timeToDeath,
                                            agentbio$timeToNextCheckpoint,
                                            agentbio$nextCheckpoint)
      
      agentbio$nextEvent       = evaluation[["next_event"]]
      agentbio$timeToNextEvent = evaluation[["time_to_event"]]
      agentbio$nextEventTime   = current_event$eventTime + agentbio$timeToNextEvent
    }
    
    
    # Step 5. Update agent bio on nextHealthState based on nextEvent
    if (agentbio$nextCheckpoint %in% c(TREAT.DEFAULT)){
      agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),   
                                         agentbio$nextEvent,
                                         str_replace(current_event$healthStatusAtEvent, "on", "post"))
      
    } else if (agentbio$nextCheckpoint == TREAT.LTBI_COMPLETE_NOTCURED){
      if (agentbio$nextEvent %in% ALL.CAUSE.DEATH){
        agentbio$nextHealthState = agentbio$nextEvent
      } else if (agentbio$nextEvent == TB){
        agentbio$nextHealthState = "TBonLTBItx"
      } else {
        agentbio$nextHealthState = str_replace(current_event$healthStatusAtEvent, "on", "post")
      }
      
    } else if (agentbio$nextCheckpoint == TREAT.INCIPIENT_COMPLETE_NOTCURED){
      if (agentbio$nextEvent %in% ALL.CAUSE.DEATH){
        agentbio$nextHealthState = agentbio$nextEvent
      } else if (agentbio$nextEvent == TB){
        agentbio$nextHealthState = "TBonINCIPtx"
      } else {
        agentbio$nextHealthState = str_replace(current_event$healthStatusAtEvent, "on", "post")
      }
      
    } else if (agentbio$nextCheckpoint %in% c(TREAT.LTBI_COMPLETE_CURED, TREAT.INCIPIENT_COMPLETE_CURED)){
      
      if (current_event$healthStatusAtEvent == TBonLTBItx){
        agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),   
                                           agentbio$nextEvent,
                                           "LTBIpostLTBItx")
      } else if(current_event$healthStatusAtEvent == TBonINCIPtx){
        agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),   
                                           agentbio$nextEvent,
                                           "LTBIpostINCIPtx")
      } else {
        agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),   
                                           agentbio$nextEvent,
                                           str_replace(current_event$healthStatusAtEvent, "on", "post"))
      }
    } else if (agentbio$nextCheckpoint %in% c(TREAT.TB_COMPLETE, TREAT.TB.REP_COMPLETE, TREAT.TB.GUARANTEE_COMPLETE)){
      
      if (current_event$healthStatusAtEvent == TBonTBtx){
        agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),   
                                           agentbio$nextEvent,
                                           "LTBIpostTBtx")
      } else {
        agentbio$nextHealthState = if_else(agentbio$nextEvent %in% c(TB, ALL.CAUSE.DEATH),   
                                           agentbio$nextEvent,
                                           str_replace(current_event$healthStatusAtEvent, "on", "post"))
      }
    }  
    
    # Step 6. Update agent bio on nextCheckPoint and timeToNextCheckpoint based on nextEvent
    if (agentbio$nextEvent == TB){
      agentbio$nextCheckpoint       = agentbio$nextCheckpoint
      agentbio$timeToNextCheckpoint = agentbio$timeToNextCheckpoint - agentbio$timeToNextEvent
    } else if (agentbio$nextEvent %in% ALL.CAUSE.DEATH){
      agentbio$nextCheckpoint       = END
      agentbio$timeToNextCheckpoint = NA 
    } else if (agentbio$nextEvent %in% c(TREAT.TB_DEFAULT, TREAT.TB.REP_DEFAULT, TREAT.TB.GUARANTEE_DEFAULT)){
      agentbio$nextCheckpoint       = TREAT.TB.REP
      agentbio$timeToNextCheckpoint = TB_Tx_to_newTx
      agentbio$randomNumbers[["fail_comp_tx"]] = runif(1, 0, 1) # resample, otherwise will keep failing the treatment
    } else if (agentbio$nextEvent %in% c(TREAT.LTBI_COMPLETE_NOTCURED, TREAT.INCIPIENT_COMPLETE_NOTCURED) & 
               current_event$healthStatusAtEvent %in% c(TBonLTBItx, TBonINCIPtx)){
      agentbio$nextCheckpoint       = TREAT.TB.GUARANTEE
      agentbio$timeToNextCheckpoint = ifelse ((agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime) < 0, 
                                              0, agentbio$timeToTB.bl + TB_Symp_to_Dx - agentbio$nextEventTime)
    } else {
      agentbio$nextCheckpoint       = MONITOR
      agentbio$timeToNextCheckpoint = NA
    }
    
    
    # Step 7. update eventHistory by appending current_event to eventHistory before exiting this checkpoint
    agentbio$eventHistory = append_event(agentbio$eventHistory, current_event)
  }
  
  
  if (agentbio$nextEvent %in% ALL.CAUSE.DEATH){
    # Record current event
    current_event = currentEvent
    current_event$eventTime           = agentbio$nextEventTime
    current_event$eventName           = agentbio$nextEvent
    current_event$healthStatusAtEvent = agentbio$nextHealthState
    current_event$ageAtEvent          = agentbio$eventHistory[1,4] + current_event$eventTime
    
    # update eventHistory by appending current_event to eventHistory
    agentbio$eventHistory = append_event(agentbio$eventHistory, current_event)
  }
  
  
  return(agentbio)
}


# -----------------------------------------------------------------------------
# A function to extract eventHistory from sim_agent output 
# -----------------------------------------------------------------------------


get_agent_record = function(agent, strategy, agentID){
  
  # simulate agent
  while(agent$nextCheckpoint != END){
    agent = agent %>% sim_agent(., strategy, agentID = agentID)
  }
  
  # check simulation output
  print(paste("Event history for agent",agentID,"under strategy",strategy,":"))
  
  # extract and format agent event history 
  
  agent_record = agent$eventHistory
  
  colnames(agent_record) = c("Time", "Event", "HealthStatus","Age")
  agent_record[,"Time"]= as.numeric(agent_record[,"Time"])
  agent_record[,"Age"] = as.numeric(agent_record[,"Age"])
  
  agent_record$AgentID = agentID
  agent_record$Strategy = strategy
  
  agent_record = agent_record[,c("Strategy", "AgentID", "Event", "Time", "HealthStatus","Age")]
  
  return(agent_record)
}


# -----------------------------------------------------------------------------
# A function to simulate all agents from the same country-of-origin 
# -----------------------------------------------------------------------------

sim_cohort = function(initialized.cohort, strategy){
  
  agent_records = NULL
  
  agent_records = mclapply(1:length(initialized.cohort), function(x) get_agent_record(agent = initialized.cohort[[x]], strategy = strategy, agentID = x), mc.cores = numCores)
  # agent_records = future_lapply(1:length(initialized.cohort), function(x) get_agent_record(agent = initialized.cohort[[x]], strategy = strategy, agentID = x), future.seed = TRUE)
  # agent_records = lapply(1:length(initialized.cohort), function(x) get_agent_record(agent = initialized.cohort[[x]], strategy = strategy, agentID = x))
  
  # agent_records = rlist::list.rbind(agent_records)
  agent_records = do.call("rbind", agent_records)
  
  
  return(agent_records)
}


# -----------------------------------------------------------------------------
# A function to simulate all countries 
# -----------------------------------------------------------------------------

get_sim_results = function(initialized.cohort.list, strategies){
  
  if(0 %in% strategies){
    full_record_0 = sim_cohort(initialized.cohort = initialized.cohort.list[[1]], strategy = 0) # No testing
  } else {full_record_0 = NA}
  
  if(1 %in% strategies){
    full_record_1 = sim_cohort(initialized.cohort = initialized.cohort.list[[2]], strategy = 1) # IGRA_TB
  } else {full_record_1 = NA}
  
  if(2 %in% strategies){
    full_record_2 = sim_cohort(initialized.cohort = initialized.cohort.list[[3]], strategy = 2) # IGRA_RNA_TB
  } else {full_record_2 = NA}
  
  if(3 %in% strategies){
    full_record_3 = sim_cohort(initialized.cohort = initialized.cohort.list[[4]], strategy = 3) # RNA_TB
  } else {full_record_3 = NA}
  
  return(list(full_record_0 = full_record_0, full_record_1 = full_record_1, full_record_2 = full_record_2, full_record_3 = full_record_3))
}


