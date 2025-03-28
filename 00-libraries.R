
# This R script will always be called first.
library(tidyverse)
library(hash)

library(pryr)

library(dampack)
# library(kableExtra)
# 
# library(future)
# library(future.apply)
# plan(multisession)
# plan(multicore)

# library(ggplot2)
# library(ggpubr)

library(parallel)
numCores =  8



# install.packages("pacman", repos="http://cran.r-project.org")

# pacman::p_load(tidyverse, hash,
#                proftools, profvis,
#                dampack, kableExtra,
#                ggplot2, ggpubr,
#                parallel)

