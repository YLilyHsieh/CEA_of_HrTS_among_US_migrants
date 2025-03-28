# This file contains functions that assist with results compilation

# NUMBER OF EVENTS ============================================================

get_country_event_tab = function(country_res, country_code){
  
  full_record = do.call("rbind", country_res)
  
  event.tab_by_strategy = full_record %>% 
    count(Strategy, Event, HealthStatus) %>% 
    mutate(Country = country_code)
  
  return(event.tab_by_strategy)
}

get_full_event_tab = function(country_res_all, writeCSV){
  
  full_tab = lapply(1:length(CountryNames), function(x) get_country_event_tab(country_res = country_res_all[[x]],
                                                                              country_code = CountryNames[x]))
  full_tab = do.call("rbind", full_tab)
  
  if (writeCSV == TRUE){
    write.csv(full_tab, paste0(ResultsDir_SummRes,"/full_event_tab_Event.csv"))
    print(paste0("writing ",ResultsDir_SummRes,"/full_event_tab_Event.csv "," <-- done."))
  } else {
    return(full_tab)
  }
}


augmented_tab = function(tab_data){
  temp1 = tab_data %>% filter(Country %in% c("MEX", "IND", "CHN")) %>%
    mutate(n = n*10)
  
  temp2 = tab_data %>% filter(Country %in% c("PHL")) %>%
    mutate(n = n*2)
  
  temp3 = tab_data %>% filter(!(Country %in% c("MEX", "IND", "CHN", "PHL")))
  
  final = rbind(temp1, temp2, temp3)
  
  final = final %>%
    filter(Strategy %in% c(0, 1, 2, 3))
  
  return(final)
}

# TIME-TO-EVENT DISTRIBUTIONS =================================================
get_country_TTE = function(country_res, country_code, health_outcome){
  
  full_record = do.call("rbind", country_res)
  
  time2event_by_strategy = full_record %>% 
    filter(Event %in% health_outcome) %>%
    mutate(Country = country_code)
  
  return(time2event_by_strategy)
}

get_full_TTE = function(country_res_all, health_outcome, writeCSV){
  
  full_TTE = lapply(1:length(CountryNames), function(x) get_country_TTE(country_res = country_res_all[[x]],
                                                                        country_code = CountryNames[x],
                                                                        health_outcome = health_outcome))
  
  full_TTE = do.call("rbind", full_TTE)
  
  if (writeCSV == TRUE){
    write.csv(full_TTE, paste0(ResultsDir_SummRes,"/full_model_TTE_",health_outcome,".csv"), row.names = F)
    print(paste0("writing ",ResultsDir_SummRes,"/full_model_TTE_",health_outcome,".csv "," <-- done."))
  } else {
    return(full_TTE)
  }
}

augmented_TTE = function(TTE_data){
  temp1 = TTE_data %>% filter(Country %in% c("MEX", "IND", "CHN")) %>% 
    slice(rep(1:n(), each = 10))
  
  temp2 = TTE_data %>% filter(Country %in% c("PHL")) %>%
    slice(rep(1:n(), each = 2))
  
  temp3 = TTE_data %>% filter(!(Country %in% c("MEX", "IND", "CHN", "PHL")))
  
  final = rbind(temp1, temp2, temp3)
  
  # final = final %>% filter(!(Time == 0))
  
  return(final)
}
# COST-EFFECTIVENESS ANALYSIS =================================================

# Time Rewards 
get_effectiveness = function(strategy_record, reward, r){
  
  # r = annual discount rate
  
  r_inst = log(1+r)
  
  dat = strategy_record %>%
    group_by(AgentID) %>% 
    mutate(change = case_when(
      HealthStatus != lag(HealthStatus) ~ TRUE,
      TRUE ~ FALSE),
      n_change = cumsum(change))
  
  dat2 = dat %>% 
    group_by(AgentID, n_change) %>%
    mutate(start_date = min(Time))
  
  dat3 = dat2 %>% 
    group_by(AgentID, HealthStatus) %>%
    distinct(start_date)
  
  dat4 = dat3 %>% 
    group_by(AgentID) %>%
    mutate(end_date = lead(start_date),
           length = end_date - start_date,
           QALY = length*reward[HealthStatus],
           disQALY = reward[HealthStatus]*(exp(-r_inst*end_date) - exp(-r_inst*start_date))/-r_inst) %>%
    summarise(QALY_sum = sum(QALY, na.rm = T),
              disQALY_sum = sum(disQALY, na.rm = T))
  
  out = ifelse(r == 0, mean(dat4$QALY_sum, na.rm = T), mean(dat4$disQALY_sum, na.rm = T))
  
  return(out)
}

get_effectiveness2 = function(strategy_record, reward, r){
  
  # r = annual discount rate
  
  r_inst = log(1+r)
  
  dat = strategy_record %>%
    group_by(AgentID) %>% 
    mutate(change = case_when(
      HealthStatus != lag(HealthStatus) ~ TRUE,
      TRUE ~ FALSE),
      n_change = cumsum(change))
  
  dat2 = dat %>% 
    group_by(AgentID, n_change) %>%
    mutate(start_date = min(Time))
  
  dat3 = dat2 %>% 
    group_by(AgentID, HealthStatus) %>%
    distinct(start_date)
  
  dat4 = dat3 %>% 
    group_by(AgentID) %>%
    mutate(end_date = lead(start_date),
           length = end_date - start_date) %>%
    mutate(HealthStatus = case_when(HealthStatus %in% c(TBonTBtx, TBpostLTBItxonTBtx, TBpostINCIPtxonTBtx) & length == TTE_TB_TX_DEFAULT ~ onTBtxDefault,
            TRUE ~ HealthStatus))
  
  dat5 = dat4 %>% 
    group_by(AgentID) %>%
    mutate(QALY = length*reward[HealthStatus],
           disQALY = reward[HealthStatus]*(exp(-r_inst*end_date) - exp(-r_inst*start_date))/-r_inst) %>%
    summarise(QALY_sum = sum(QALY, na.rm = T),
              disQALY_sum = sum(disQALY, na.rm = T))
  
  if (r == 0){
    out = dat5[,c("AgentID","QALY_sum")]
  } else {
    out = dat5[,c("AgentID", "disQALY_sum")]
  }
  
  #out = ifelse(r == 0, dat4[,c("AgentID", "QALY_sum")], dat4[,c("AgentID", "disQALY_sum")]) # return a data frame 
  
  return(out)
}

get_costs_agent = function(strategy_record, agentID, HC_costs_ctlg, nonHC_costs_ctlg, r){
  
  print(agentID)
  dat = strategy_record %>% filter(AgentID == agentID)
  
  ageb4deadage = ifelse ( floor(dat$Age[dat$Event == "DEAD"]) < 1, 0, floor(dat$Age[dat$Event == "DEAD"]) - 1 )

  # extract values from catalogues
  HC_costs.vec    = c(HC_costs_ctlg[which(HC_costs_ctlg$age == dat$Age[dat$Event == "Flag"]):
                              which(HC_costs_ctlg$age == ageb4deadage),]$costs, 
                    HC_costs_ctlg[which(HC_costs_ctlg$age == floor(dat$Age[dat$Event == "DEAD"])),]$costs *  
                      (dat$Age[dat$Event == "DEAD"] - floor(dat$Age[dat$Event == "DEAD"])))
              
  nonHC_costs.vec = c(nonHC_costs_ctlg[which(nonHC_costs_ctlg$age == dat$Age[dat$Event == "Flag"]):
                            which(nonHC_costs_ctlg$age == ageb4deadage),]$costs,  
                      nonHC_costs_ctlg[which(nonHC_costs_ctlg$age == floor(dat$Age[dat$Event == "DEAD"])),]$costs *  
                        (dat$Age[dat$Event == "DEAD"] - floor(dat$Age[dat$Event == "DEAD"])))
  
   
  # apply discounting
  HC_costs.vec.dist    = HC_costs.vec / (1+r)^c(min(dat$Time):floor(max(dat$Time))) 
  nonHC_costs.vec.dist = nonHC_costs.vec / (1+r)^c(min(dat$Time):floor(max(dat$Time))) 
   
  # summing
  HC_costs      = sum(HC_costs.vec)
  HC_costs.dist = sum(HC_costs.vec.dist)
  
  nonHC_costs      = sum(nonHC_costs.vec)
  nonHC_costs.dist = sum(nonHC_costs.vec.dist)
  
  return(c(HC_costs    = HC_costs,    HC_costs.dist    = HC_costs.dist, 
           nonHC_costs = nonHC_costs, nonHC_costs.dist = nonHC_costs.dist))
}

get_costs = function(strategy_record, TBHCcost_profile, TBnonHCcost_profile, HC_costs_ctlg, nonHC_costs_ctlg, r, HCperspective){
  
  # r = annual discount rate 
  
  dat1 = strategy_record %>%
    mutate(TBHC_costs         = TBHCcost_profile[Event],
           TBHC_costs.dist    = TBHC_costs/(1+r)^Time,
           TBnonHC_costs      = TBnonHCcost_profile[Event],
           TBnonHC_costs.dist = TBnonHC_costs/(1+r)^Time,
           ) %>%
    group_by(AgentID) %>% 
    summarise(TBHC_costs         = sum(TBHC_costs),
              TBHC_costs.dist    = sum(TBHC_costs.dist),
              TBnonHC_costs      = sum(TBnonHC_costs),
              TBnonHC_costs.dist = sum(TBnonHC_costs.dist))
  
  dat2 = lapply(1:max(strategy_record$AgentID), function(x) get_costs_agent(strategy_record = strategy_record, 
                                                                                     agentID = x, 
                                                                                     HC_costs_ctlg = HC_costs_ctlg,
                                                                                     nonHC_costs_ctlg = nonHC_costs_ctlg, 
                                                                                     r = r))
  dat2 = do.call("rbind", dat2)
  
  cost_df = cbind(dat1, dat2)
  
  out.healthcare = mean(rowSums(cost_df[,c("TBHC_costs.dist", "HC_costs.dist")]))
  out.societal   = mean(rowSums(cost_df[,c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist")]))
  
  out = ifelse(HCperspective, out.healthcare, out.societal)
 
  return(out)
  
}

get_costs_itemized = function(strategy_record, strategy_no, TBHCcost_profile, TBnonHCcost_profile, HC_costs_ctlg, nonHC_costs_ctlg, r, HCperspective){
  
  # r = annual discount rate 
  
  dat1 = strategy_record %>%
    mutate(TBHC_costs         = TBHCcost_profile[Event],
           TBHC_costs.dist    = TBHC_costs/(1+r)^Time,
           TBnonHC_costs      = TBnonHCcost_profile[Event],
           TBnonHC_costs.dist = TBnonHC_costs/(1+r)^Time,
    ) %>%
    group_by(AgentID) %>% 
    summarise(TBHC_costs         = sum(TBHC_costs),
              TBHC_costs.dist    = sum(TBHC_costs.dist),
              TBnonHC_costs      = sum(TBnonHC_costs),
              TBnonHC_costs.dist = sum(TBnonHC_costs.dist))
  
  dat2 = lapply(1:max(strategy_record$AgentID), function(x) get_costs_agent(strategy_record = strategy_record, 
                                                                            agentID = x, 
                                                                            HC_costs_ctlg = HC_costs_ctlg,
                                                                            nonHC_costs_ctlg = nonHC_costs_ctlg, 
                                                                            r = r))
  dat2 = do.call("rbind", dat2)
  
  cost_df = cbind(strategy = strategy_no, dat1, dat2)
  
  # out = colMeans(cost_df[2:9])
  
  # out.healthcare = mean(rowSums(cost_df[,c("TBHC_costs.dist", "HC_costs.dist")]))
  # out.societal   = mean(rowSums(cost_df[,c("TBHC_costs.dist", "HC_costs.dist", "TBnonHC_costs.dist", "nonHC_costs.dist")]))
  #out = ifelse(HCperspective, out.healthcare, out.societal)
  
  return(cost_df)
  
}

get_prod_agent = function(strategy_record, agentID, prod_ctlg, r){
  
  dat = strategy_record %>% filter(AgentID == agentID)
  
  ageb4deadage = ifelse ( floor(dat$Age[dat$Event == "DEAD"]) < 1, 0, floor(dat$Age[dat$Event == "DEAD"]) - 1 )
  
  # extract values from catalogues
  prod.vec        = c(prod_ctlg[which(prod_ctlg$age == dat$Age[dat$Event == "Flag"]):
                                  which(prod_ctlg$age == ageb4deadage),]$productivity, 
                      prod_ctlg[which(prod_ctlg$age == floor(dat$Age[dat$Event == "DEAD"])),]$productivity * 
                        (dat$Age[dat$Event == "DEAD"] - floor(dat$Age[dat$Event == "DEAD"])))
  
  # apply discounting
  prod.vec.dist        = prod.vec / (1+r)^c(min(dat$Time):floor(max(dat$Time))) 
  
  # summing
  prod      = sum(prod.vec)
  prod.dist = sum(prod.vec.dist)
  
  return(c(prod = prod, prod.dist = prod.dist))
}

get_prod = function(strategy_record, prod_ctlg, r){
  
  # r = annual discount rate 
  
  dat = lapply(1:max(strategy_record$AgentID), function(x) get_prod_agent(strategy_record = strategy_record, 
                                                                                     agentID = x, 
                                                                                     prod_ctlg = prod_ctlg, 
                                                                                     r = 0.03))
  dat = do.call("rbind", dat)
  out = mean(dat[,"prod.dist"])
  
  return(out)
}

get_prod2 = function(strategy_record, prod_ctlg, r){
  
  # r = annual discount rate 
  
  dat = lapply(1:max(strategy_record$AgentID), function(x) get_prod_agent(strategy_record = strategy_record, 
                                                                          agentID = x, 
                                                                          prod_ctlg = prod_ctlg, 
                                                                          r = 0.03))
  dat = do.call("rbind", dat)

  if (r == 0){
    out = data.frame(AgentID = c(1:max(strategy_record$AgentID)), 
                     prod = dat[,"prod"])
  } else {
    out = data.frame(AgentID = c(1:max(strategy_record$AgentID)),
                     prod.dist = dat[,"prod.dist"])
  }
  
  return(out)
}

get_mean_eff = function(country_record, strategy_nums, reward, r){
  eff = sapply(strategy_nums + 1, function(x) get_effectiveness(country_record[[x]], reward, r))
  
  return(eff)
}

get_mean_eff2 = function(country_record, strategy_nums, reward, r){
  eff = lapply(strategy_nums + 1, function(x) get_effectiveness2(country_record[[x]], reward, r))
  
  
  eff = do.call("rbind", eff)
  out = cbind(strategy = rep(c(0:3), each = nrow(eff)/4), eff)
  return(out)
}

get_mean_cost = function(country_record, strategy_nums, TBHCcost_profile, TBnonHCcost_profile, HC_costs_ctlg, nonHC_costs_ctlg, r, HCperspective){
  cost = sapply(strategy_nums + 1, 
                function(x) get_costs(country_record[[x]], x, TBHCcost_profile, TBnonHCcost_profile, 
                                      HC_costs_ctlg, nonHC_costs_ctlg, r, HCperspective))
  
  return(cost)
}

get_mean_cost_itemized = function(country_record, strategy_nums, TBHCcost_profile, TBnonHCcost_profile, HC_costs_ctlg, nonHC_costs_ctlg, r){
  cost = lapply(strategy_nums + 1, 
                function(x) get_costs_itemized(country_record[[x]], x, TBHCcost_profile, TBnonHCcost_profile, 
                                      HC_costs_ctlg, nonHC_costs_ctlg, r, HCperspective))
  
  cost = do.call("rbind", cost)
  
  # colnames(cost) = c("No_testing", "IGRA_TB", "IGRA_RNA_TB", "RNA_TB")
  return(cost)
}

get_mean_prod = function(country_record, strategy_nums, prod_ctlg, r){
  prod = sapply(strategy_nums + 1, function(x) get_prod(country_record[[x]], prod_ctlg, r))
  
  return(prod)
}

get_mean_prod2 = function(country_record, strategy_nums, prod_ctlg, r){
  prod = lapply(strategy_nums + 1, function(x) get_prod2(country_record[[x]], prod_ctlg, r))
  
  prod = do.call("rbind", prod)
  out = cbind(strategy = rep(c(0:3), each = nrow(prod)/4), 
              prod)
  return(out)
}
# TREATMENT CASCADE ===========================================================

get_cascade = function(data, strategy, disease_state){
  # disease_state, character, "LTBI", "Incipient TB", "TB disease"
  # prop, logical, default is T. If FALSE, then absolute values are reported
  
  tab = data %>% 
    mutate( prop = .$n/max(.$n))
  
  if (disease_state == LTBI){
    
    if (strategy == 1){
      
      df = tab %>% 
        filter(Strategy == 1 & Event %in% c("Flag", "TEST.IGRA", "TEST.RNA", 
                                            "TREAT.LTBI_INIT.FAIL", "TREAT.LTBI_INIT.SUCCESS", "TREAT.LTBI_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>%
        add_row(Event = "TEST.RNA.neg", .after = which(.$Event == "TEST.RNA")) %>% 
        fill(Strategy, Country) 
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.RNA",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.RNA",]$prop
      
      if ( nrow(df[df$Event == "TREAT.LTBI_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.RNA.neg",]$n = NA
        df[df$Event == "TEST.RNA.neg",]$prop = NA
      } else {
        df[df$Event == "TEST.RNA.neg",]$n  = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$n + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$n
        df[df$Event == "TEST.RNA.neg",]$prop  = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$prop + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.LTBI_INIT.FAIL") %>%
        arrange(desc(prop)) 
      
    } else if (strategy == 2){
      
      df = tab %>%
        filter(Strategy == 2 & Event %in% c("Flag", "TEST.IGRA", "TEST.XRAY", "TEST.RNA", 
                                            "TREAT.LTBI_INIT.FAIL", "TREAT.LTBI_INIT.SUCCESS", "TREAT.LTBI_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.neg", .after = which(.$Event == "TEST.XRAY")) %>%
        add_row(Event = "TEST.RNA.neg",  .after = which(.$Event == "TEST.RNA")) %>%
        fill(Strategy, Country)
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      df[df$Event == "TEST.XRAY.neg",]$n = df[df$Event == "TEST.RNA",]$n
      df[df$Event == "TEST.XRAY.neg",]$prop = df[df$Event == "TEST.RNA",]$prop
      
      if ( nrow(df[df$Event == "TREAT.LTBI_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.RNA.neg",]$n = NA
        df[df$Event == "TEST.RNA.neg",]$prop = NA
      } else {
        df[df$Event == "TEST.RNA.neg",]$n  = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$n + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$n
        df[df$Event == "TEST.RNA.neg",]$prop  = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$prop + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.LTBI_INIT.FAIL") %>%
        arrange(desc(prop))
      
    } else if (strategy == 3){
      
      df = tab %>%
        filter(Strategy == 3 & Event %in% c("Flag", "TEST.IGRA", "TEST.RNA_XRAY",
                                            "TREAT.LTBI_INIT.FAIL", "TREAT.LTBI_INIT.SUCCESS", "TREAT.LTBI_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.neg_RNA.neg", .after = which(.$Event == "TEST.RNA_XRAY")) %>%
        fill(Strategy, Country)
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.RNA_XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.RNA_XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.LTBI_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.neg_RNA.neg",]$n = NA
        df[df$Event == "TEST.XRAY.neg_RNA.neg",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.neg_RNA.neg",]$n = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$n + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.neg_RNA.neg",]$prop = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$prop + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.LTBI_INIT.FAIL") %>%
        arrange(desc(prop))
      
    } else if (strategy == 4){
      
      df = NA
      
    } else if (strategy == 5) {
      
      df = tab %>%
        filter(Strategy == 5 & Event %in% c("Flag", "TEST.IGRA", "TEST.XRAY",
                                            "TREAT.LTBI_INIT.FAIL", "TREAT.LTBI_INIT.SUCCESS", "TREAT.LTBI_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.neg", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country)
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.LTBI_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.neg",]$n = NA
        df[df$Event == "TEST.XRAY.neg",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.neg",]$n = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$n + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.neg",]$prop = df[df$Event == "TREAT.LTBI_INIT.FAIL",]$prop + df[df$Event == "TREAT.LTBI_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.LTBI_INIT.FAIL") %>%
        arrange(desc(prop))
    }
    
    
  } else if (disease_state == "Incipient TB") {
    
    if (strategy == 1){
      
      df = tab %>% 
        filter(Strategy == 1 & Event %in% c("Flag", "TEST.IGRA", "TEST.RNA", "TEST.XRAY",
                                            "TREAT.INCIPIENT_INIT.FAIL", "TREAT.INCIPIENT_INIT.SUCCESS", "TREAT.INCIPIENT_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>%
        add_row(Event = "TEST.RNA.pos", .after = which(.$Event == "TEST.RNA")) %>% 
        add_row(Event = "TEST.XRAY.neg", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country) 
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.RNA",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.RNA",]$prop
      df[df$Event == "TEST.RNA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.RNA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.neg",]$n = NA
        df[df$Event == "TEST.XRAY.neg",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.neg",]$n  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$n + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.neg",]$prop  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$prop + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.INCIPIENT_INIT.FAIL") %>%
        arrange(desc(prop))
      
    } else if (strategy == 2){
      
      df = tab %>%
        filter(Strategy == 2 & Event %in% c("Flag", "TEST.IGRA", "TEST.XRAY", "TEST.RNA", 
                                            "TREAT.INCIPIENT_INIT.FAIL", "TREAT.INCIPIENT_INIT.SUCCESS", "TREAT.INCIPIENT_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.neg", .after = which(.$Event == "TEST.XRAY")) %>%
        add_row(Event = "TEST.RNA.pos",  .after = which(.$Event == "TEST.RNA")) %>%
        fill(Strategy, Country)
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      df[df$Event == "TEST.XRAY.neg",]$n = df[df$Event == "TEST.RNA",]$n
      df[df$Event == "TEST.XRAY.neg",]$prop = df[df$Event == "TEST.RNA",]$prop
      
      if ( nrow(df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.RNA.pos",]$n = NA
        df[df$Event == "TEST.RNA.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.RNA.pos",]$n  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$n + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$n
        df[df$Event == "TEST.RNA.pos",]$prop  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$prop + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.INCIPIENT_INIT.FAIL") %>%
        arrange(desc(prop))
      
    } else if (strategy == 3){
      
      df = tab %>%
        filter(Strategy == 3 & Event %in% c("Flag", "TEST.IGRA", "TEST.RNA_XRAY",
                                            "TREAT.INCIPIENT_INIT.FAIL", "TREAT.INCIPIENT_INIT.SUCCESS", "TREAT.INCIPIENT_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.neg_RNA.pos", .after = which(.$Event == "TEST.RNA_XRAY")) %>%
        fill(Strategy, Country)
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.RNA_XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.RNA_XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.neg_RNA.pos",]$n = NA
        df[df$Event == "TEST.XRAY.neg_RNA.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.neg_RNA.pos",]$n  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$n + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.neg_RNA.pos",]$prop  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$prop + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.INCIPIENT_INIT.FAIL") %>%
        arrange(desc(prop))
      
    } else if (strategy == 4){
      
      df = tab %>%
        filter(Strategy == 4 & Event %in% c("Flag", "TEST.RNA", "TEST.XRAY",
                                            "TREAT.INCIPIENT_INIT.FAIL", "TREAT.INCIPIENT_INIT.SUCCESS", "TREAT.INCIPIENT_COMPLETE")) %>%
        add_row(Event = "TEST.RNA.pos", .after = which(.$Event == "TEST.RNA")) %>% 
        add_row(Event = "TEST.XRAY.neg", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country)
      
      df[df$Event == "TEST.RNA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.RNA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.neg",]$n = NA
        df[df$Event == "TEST.XRAY.neg",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.neg",]$n  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$n + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.neg",]$prop  = df[df$Event == "TREAT.INCIPIENT_INIT.FAIL",]$prop + df[df$Event == "TREAT.INCIPIENT_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.INCIPIENT_INIT.FAIL") %>%
        arrange(desc(prop))
      
    } else if (strategy == 5) {
      
      df = NA
    }
    
  } else if (disease_state == "TB disease"){
    
    if (strategy == 1){
      
      df = tab %>% 
        filter(Strategy == 1 & Event %in% c("Flag", "TEST.IGRA", "TEST.RNA", "TEST.XRAY",
                                            "TREAT.TB_INIT.FAIL", "TREAT.TB_INIT.SUCCESS", "TREAT.TB_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>%
        add_row(Event = "TEST.RNA.pos", .after = which(.$Event == "TEST.RNA")) %>% 
        add_row(Event = "TEST.XRAY.pos", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country) 
      
      df2 = tab %>%
        filter(Strategy == 1 & Event %in% c("TREAT.TB.GUARANTEE", "TREAT.TB.GUARANTEE_COMPLETE"))
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.RNA",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.RNA",]$prop
      df[df$Event == "TEST.RNA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.RNA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.TB_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.TB_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.pos",]$n = NA
        df[df$Event == "TEST.XRAY.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.pos",]$n  = df[df$Event == "TREAT.TB_INIT.FAIL",]$n + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.pos",]$prop  = df[df$Event == "TREAT.TB_INIT.FAIL",]$prop + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$prop
      }
      
      
      df = df %>% filter(!Event == "TREAT.TB_INIT.FAIL") %>%
        replace_na(list(n = 0, prop = 0)) %>%
        arrange(desc(prop)) %>%
        rbind(., df2)
      
    } else if (strategy == 2){
      
      df = tab %>%
        filter(Strategy == 2 & Event %in% c("Flag", "TEST.IGRA", "TEST.XRAY", 
                                            "TREAT.TB_INIT.FAIL", "TREAT.TB_INIT.SUCCESS", "TREAT.TB_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.pos", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country)
      
      df2 = tab %>%
        filter(Strategy == 2 & Event %in% c("TREAT.TB.GUARANTEE", "TREAT.TB.GUARANTEE_COMPLETE"))
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.TB_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.TB_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.pos",]$n = NA
        df[df$Event == "TEST.XRAY.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.pos",]$n  = df[df$Event == "TREAT.TB_INIT.FAIL",]$n + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.pos",]$prop  = df[df$Event == "TREAT.TB_INIT.FAIL",]$prop + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$prop
      }
      
      df = df %>% filter(!Event == "TREAT.TB_INIT.FAIL") %>%
        replace_na(list(n = 0, prop = 0)) %>%
        arrange(desc(prop)) %>%
        rbind(., df2)
      
    } else if (strategy == 3){
      
      df = tab %>%
        filter(Strategy == 3 & Event %in% c("Flag", "TEST.IGRA", "TEST.RNA_XRAY",
                                            "TREAT.TB_INIT.FAIL", "TREAT.TB_INIT.SUCCESS", "TREAT.TB_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.pos", .after = which(.$Event == "TEST.RNA_XRAY")) %>%
        fill(Strategy, Country)
      
      df2 = tab %>%
        filter(Strategy == 3 & Event %in% c("TREAT.TB.GUARANTEE", "TREAT.TB.GUARANTEE_COMPLETE"))
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.RNA_XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.RNA_XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.TB_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.TB_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.pos",]$n = NA
        df[df$Event == "TEST.XRAY.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.pos",]$n  = df[df$Event == "TREAT.TB_INIT.FAIL",]$n + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.pos",]$prop  = df[df$Event == "TREAT.TB_INIT.FAIL",]$prop + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$prop
      }
      
      
      df = df %>% filter(!Event == "TREAT.TB_INIT.FAIL") %>%
        replace_na(list(n = 0, prop = 0)) %>%
        arrange(desc(prop)) %>%
        rbind(., df2)
      
    } else if (strategy == 4){
      
      df = tab %>%
        filter(Strategy == 4 & Event %in% c("Flag", "TEST.RNA", "TEST.XRAY",
                                            "TREAT.TB_INIT.FAIL", "TREAT.TB_INIT.SUCCESS", "TREAT.TB_COMPLETE")) %>%
        add_row(Event = "TEST.RNA.pos", .after = which(.$Event == "TEST.RNA")) %>% 
        add_row(Event = "TEST.XRAY.pos", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country)
      
      df2 = tab %>%
        filter(Strategy == 4 & Event %in% c("TREAT.TB.GUARANTEE", "TREAT.TB.GUARANTEE_COMPLETE"))
      
      df[df$Event == "TEST.RNA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.RNA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.TB_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.TB_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.pos",]$n = NA
        df[df$Event == "TEST.XRAY.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.pos",]$n  = df[df$Event == "TREAT.TB_INIT.FAIL",]$n + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.pos",]$prop  = df[df$Event == "TREAT.TB_INIT.FAIL",]$prop + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$prop
      }
      
      
      df = df %>% filter(!Event == "TREAT.TB_INIT.FAIL") %>%
        replace_na(list(n = 0, prop = 0)) %>%
        arrange(desc(prop)) %>%
        rbind(., df2)
      
    } else if (strategy == 5) {
      
      df = tab %>%
        filter(Strategy == 5 & Event %in% c("Flag", "TEST.IGRA", "TEST.XRAY",
                                            "TREAT.TB_INIT.FAIL", "TREAT.TB_INIT.SUCCESS", "TREAT.TB_COMPLETE")) %>%
        add_row(Event = "TEST.IGRA.pos", .after = which(.$Event == "TEST.IGRA")) %>% 
        add_row(Event = "TEST.XRAY.pos", .after = which(.$Event == "TEST.XRAY")) %>%
        fill(Strategy, Country)
      
      df2 = tab %>%
        filter(Strategy == 5 & Event %in% c("TREAT.TB.GUARANTEE", "TREAT.TB.GUARANTEE_COMPLETE"))
      
      df[df$Event == "TEST.IGRA.pos",]$n = df[df$Event == "TEST.XRAY",]$n
      df[df$Event == "TEST.IGRA.pos",]$prop = df[df$Event == "TEST.XRAY",]$prop
      
      if ( nrow(df[df$Event == "TREAT.TB_INIT.FAIL",]) == 0 | nrow(df[df$Event == "TREAT.TB_INIT.SUCCESS",]) == 0 ){
        df[df$Event == "TEST.XRAY.pos",]$n = NA
        df[df$Event == "TEST.XRAY.pos",]$prop = NA
      } else {
        df[df$Event == "TEST.XRAY.pos",]$n  = df[df$Event == "TREAT.TB_INIT.FAIL",]$n + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$n
        df[df$Event == "TEST.XRAY.pos",]$prop  = df[df$Event == "TREAT.TB_INIT.FAIL",]$prop + df[df$Event == "TREAT.TB_INIT.SUCCESS",]$prop
      }
      
      
      df = df %>% filter(!Event == "TREAT.TB_INIT.FAIL") %>%
        replace_na(list(n = 0, prop = 0)) %>%
        arrange(desc(prop)) %>%
        rbind(., df2)
    }
  }
  
  return(df)
}

