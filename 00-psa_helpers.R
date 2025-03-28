# This file generates the table of parameter values to run PSA 

library(randtoolbox)


# Generate the PSA table using Sobol low-discrepancy quasi-random sequences -----------------------

makePSAtable = function(n, vals_cont, vals_disc){
  # n    = numeric, number of samples to draw
  # dist = list, distibutions (variable name, 'model.param', 
  #                            distribution name, 'dist.name',
  #                            distribution params, 'dist.params')
  
  # Generate a Sobol' sequence
  if (is.null(vals_disc)){
    sob = sobol(n, length(vals_cont))
  } else {
    sob = sobol(n, length(vals_cont) + length(vals_disc))
  }
  
  # Fill a matrix with the values inverted from uniform values to distributions of choice for continuous params
  samp = matrix(NA, nrow = n, ncol = length(vals_cont) + 1)
  samp[,1] = 1:n

  for (i in 1:length(vals_cont)) {
    l = vals_cont[[i]]
    dist   = l$dist.name
    params = l$dist.params
    fname  = paste("q",l$dist.name,sep="")
    samp[,i+1] = do.call(fname,c(list(p=sob[,i]),params))
  }
  
  # Convert matrix to data frame and add labels
  samp        = as.data.frame(samp)
  names(samp) = c("n", lapply(vals_cont, function(l) l$model.param))
  
  # Fill a matrix with the values inverted from uniform values to distributions of choice for discrete params
  if (!is.null(vals_disc)){
    samp2 = matrix(NA, nrow = n, ncol = length(vals_disc))
    for (i in 1:length(vals_disc)) {
      for (j in 1:n){
        if (sob[j,i] <= 0.25){samp2[j,i] = 1
        } else if (sob[j,i] > 0.25 & sob[j,i] <= 0.5) {samp2[j,i] = 2
        } else if (sob[j,i] > 0.5  & sob[j,i] <= 0.75) {samp2[j,i] = 3
        } else {samp2[j,i] = 4}
      }
    }
    samp2        = as.data.frame(samp2)
    names(samp2) = vals_disc
    
    samp = cbind(samp, samp2)
  }
  
  
  return(samp)
}

vals_cont = list(list(model.param="TB_Symp_to_Dx",
                 dist.name="lnorm",
                 dist.params=list(meanlog = 0.5, sdlog = 0.125)),
            list(model.param="IGRA_SPEC_HEALTHY",
                 dist.name="beta",
                 dist.params=list(shape1 = 191.1, shape2 = 3.9)),
            list(model.param="IGRA_SENS_LTBI",   # same as IGRA_SENS_TB
                 dist.name="beta",
                 dist.params=list(shape1 = 138.5196, shape2 = 17.1204)),
            list(model.param="PROB_FAIL_TO_INIT_LTBI_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 1325, shape2 = 4235)),
            list(model.param="PROB_FAIL_TO_COMP_LTBI_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 126, shape2 = 1168)),
            list(model.param="PROB_CURED_LTBI_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 7.16, shape2 = 4.03)),
            list(model.param="PROB_FAIL_TO_INIT_INCIPIENT_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 5.41, shape2 = 40.06)),
            list(model.param="PROB_FAIL_TO_COMP_INCIPIENT_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 9.68, shape2 = 122.99)),
            list(model.param="PROB_CURED_INCIPIENT_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 7.16, shape2 = 4.03)),
            list(model.param="PROB_FAIL_TO_COMP_TB_TX_SURVIGING_TX",
                 dist.name="beta",
                 dist.params=list(shape1 = 41.28, shape2 = 450.16)),
            list(model.param="RR_earlyTB_1",
                 dist.name="gamma",
                 dist.params=list(shape = 232.96, rate = 105.08)),
            list(model.param="RR_earlyTB_2",
                 dist.name="gamma",
                 dist.params=list(shape = 21.11, rate = 2.33)),
            list(model.param="COST_TB_TX_HC",
                 dist.name="gamma",
                 dist.params=list(shape = 70.51, scale = 297.71)),
            list(model.param="COST_TB_TX_nonHC",
                 dist.name="gamma",
                 dist.params=list(shape = 145.30, scale = 31.79)),
            list(model.param="COST_LTBI_TX_HC",
                 dist.name="gamma",
                 dist.params=list(shape = 13.46, scale = 31.81)),
            list(model.param="COST_LTBI_TX_nonHC",
                 dist.name="gamma",
                 dist.params=list(shape = 125.84, scale = 0.89)),
            list(model.param="COST_INCIP_TX_HC",
                 dist.name="gamma",
                 dist.params=list(shape = 14.86, scale = 34.59)),
            list(model.param="COST_INCIP_TX_nonHC",
                 dist.name="gamma",
                 dist.params=list(shape = 125.84, scale = 1.19)),
            list(model.param="COST_IGRA",
                 dist.name="gamma",
                 dist.params=list(shape = 59.68, scale = 1.08)),
            list(model.param="COST_XRAY",
                 dist.name="gamma",
                 dist.params=list(shape = 47.58, scale = 0.72)),
            list(model.param="REWARD_UNCONTROLLED_TB",
                 dist.name="gamma",
                 dist.params=list(shape = 34.60208, scale = 0.007225)),            
            list(model.param="REWARD_CONTROLLED_TB",
                 dist.name="gamma",
                 dist.params=list(shape = 7.9524, scale = 0.0177305)),
            list(model.param="REWARD_POST_TB",
                 dist.name="gamma",
                 dist.params=list(shape = 16, scale = 0.00075)),
            list(model.param="RNA_SENS1",
                 dist.name="unif",
                 dist.params=list(min = 0.75, max = 1.0)),
            list(model.param="RNA_SPEC",
                 dist.name="unif",
                 dist.params=list(min = 0.75, max = 1.0)),
            list(model.param="RNA_TIME_upper",
                 dist.name="unif",
                 dist.params=list(min = 2, max = 10))
)

# vals_disc = c("RNA_SENS1","RNA_SPEC")
# vals_disc = c("RNA_TIME_upper")

# n = length(vals)*10 

# Generate the PSA params table
set.seed(08122022)
PSAtable = makePSAtable(1000,vals_cont, vals_disc = NULL)
PSAtable[,"TB_Symp_to_Dx"] = log(PSAtable[,"TB_Symp_to_Dx"])
PSAtable[,"REWARD_POST_TB"]   = 1 - PSAtable[,"REWARD_POST_TB"] 
PSAtable[,"REWARD_CONTROLLED_TB"]   = 1 - PSAtable[,"REWARD_CONTROLLED_TB"] 
PSAtable[,"REWARD_UNCONTROLLED_TB"] = 1 - PSAtable[,"REWARD_UNCONTROLLED_TB"]

PSAtable$RNA_SENS1 = 0.9
PSAtable$RNA_SPEC  = 0.9
PSAtable$RNA_SENS_upper_plus = 1 -  PSAtable$RNA_SPEC
PSAtable$RNA_TIME_upper = 720/365
write.csv(PSAtable, paste0("Data/PSAtable1000.csv"), row.names = F) # base case


# # Adding test characteristics scenarios1 
# PSAtable$RNA_SENS1 = ifelse(PSAtable$RNA_SENS1 %in% c(1, 2), 0.9, 1.0)
# PSAtable$RNA_SPEC  = ifelse(PSAtable$RNA_SPEC  %in% c(1, 3), 0.9, 1.0)
# 
# PSAtable$RNA_SENS_upper_plus = 1 -  PSAtable$RNA_SPEC
# PSAtable$RNA_TIME_upper = 720/365
# 
# write.csv(PSAtable, paste0("Data/PSAtable1000_testscenar1.csv"), row.names = F) # varying test characteristics for t<=2

# Create PSAtable for test characteristics GAM
set.seed(08122022)
PSAtable = makePSAtable(1000,vals_cont, vals_disc = NULL)
PSAtable[,"TB_Symp_to_Dx"] = log(PSAtable[,"TB_Symp_to_Dx"])
PSAtable[,"REWARD_CONTROLLED_TB"]   = 1 - PSAtable[,"REWARD_CONTROLLED_TB"] 
PSAtable[,"REWARD_UNCONTROLLED_TB"] = 1 - PSAtable[,"REWARD_UNCONTROLLED_TB"]
 
PSAtable$RNA_SENS_upper_plus = 1 -  PSAtable$RNA_SPEC
write.csv(PSAtable, paste0("Data/PSAtable1000_gam.csv"), row.names = F) 

