# Thie file contains functions that help plot results 

# library(ggpubr)
 library(grid)
 
# Theme settings 
theme_bar = function(){
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           text = element_text(size = 12))
}

theme_line = function(){
  theme_bw()  + theme(plot.title = element_text(hjust = 0.5),
                      text = element_text(size = 16),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14))
}



