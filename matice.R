library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(vegan)
library(iNEXT)
library(bipartiteD3)
library(grid)
library(gridExtra)
library(betapart)
library(usedist)
library(bipartite)
library(readxl)
library(writexl)
library(openxlsx)



#loading data and prepare data -all parasitoids and caterpillars without singletons and doubletons
main_table <- read_excel("data/input/Beta_August22.xlsx")
main_table <- rename (main_table, par_morph = NEW_M_code, cat_morph = cat_whole_name, plant_sp = sample_code) 

testik <- main_table %>% 
  filter(Par_remove == "0", Remove_cats == "0", singlecat == "0", singlepara == "0") %>% 
  select(LOC2, par_morph, cat_morph, plant_sp)

Distance <- 
  read.csv2(
    here::here("~/ownCloud/000_/000_R_stat/Beta_diversity/data/input/Distance.csv"), row.names=1)  