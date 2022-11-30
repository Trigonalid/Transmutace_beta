library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(vegan)
library(iNEXT)
library(bipartiteD3)
library(bipartite)
library(betapart)
library(usedist)
library(readxl)
library(writexl)
library(openxlsx)



#loading data and prepare data -all parasitoids and caterpillars without singletons and doubletons
main_table <- read_excel("data/input/Beta_August22.xlsx")
main_table <- rename (main_table, par_morph = NEW_M_code, cat_morph = cat_whole_name, plant_sp = sample_code) 

#Main dataset - reduced
testik <- main_table %>% 
  filter(Par_remove == "0", Remove_cats == "0", singlecat == "0", singlepara == "0") %>% 
  select(LOC2, par_morph, cat_morph, plant_sp)

# distance matrix
Distance <- 
  read.csv2(
    here::here("data/input/Distance.csv"), row.names=1)


#######   Food-web interactions between caterpillars and parasitoids

#Preparing contingency table - caterpillars vs parasitoids

Elem_PC <- testik %>% 
  filter(LOC2 =="Elem") %>%
  select(cat_morph, par_morph)
Elem_PC <- as.matrix(table(Elem_PC$cat_morph, Elem_PC$par_morph))

Morox_PC <- testik %>% 
  filter(LOC2 =="Morox") %>%
  select(cat_morph, par_morph)
Morox_PC <- as.matrix(table(Morox_PC$cat_morph, Morox_PC$par_morph))

Niksek_PC <- testik %>% 
  filter(LOC2 =="Niksek") %>%
  select(cat_morph, par_morph)
Niksek_PC <- as.matrix(table(Niksek_PC$cat_morph, Niksek_PC$par_morph))

Ohu1_PC <- testik %>% 
  filter(LOC2 =="Ohu1") %>%
  select(cat_morph, par_morph)
Ohu1_PC <- as.matrix(table(Ohu1_PC$cat_morph, Ohu1_PC$par_morph))

Ohu2_PC <- testik %>% 
  filter(LOC2 =="Ohu2") %>%
  select(cat_morph, par_morph)
Ohu2_PC <- as.matrix(table(Ohu2_PC$cat_morph, Ohu2_PC$par_morph))

Utai_PC <- testik %>% 
  filter(LOC2 =="Utai") %>%
  select(cat_morph, par_morph)
Utai_PC <- as.matrix(table(Utai_PC$cat_morph, Utai_PC$par_morph))

Wamangu_PC <- testik %>% 
  filter(LOC2 =="Wamangu") %>%
  select(cat_morph, par_morph)
Wamangu_PC <- as.matrix(table(Wamangu_PC$cat_morph, Wamangu_PC$par_morph))

Wanang_PC <- testik %>% 
  filter(LOC2 =="Wanang") %>%
  select(cat_morph, par_morph)
Wanang_PC <- as.matrix(table(Wanang_PC$cat_morph, Wanang_PC$par_morph))

Yapsiei_PC <- testik %>% 
  filter(LOC2 =="Yapsiei") %>%
  select(cat_morph, par_morph)
Yapsiei_PC <- as.matrix(table(Yapsiei_PC$cat_morph, Yapsiei_PC$par_morph))

################################# run the analyses

FW_PC_res <-betalinkr_multi(webs2array(Elem_PC,Morox_PC,Niksek_PC,Ohu1_PC,Ohu2_PC,Utai_PC, Wamangu_PC, Wanang_PC, Yapsiei_PC), partitioning="commondenom", partition.st=TRUE)
FW_PC_res<-FW_PC_res %>% add_column(dataset = "Cat_Par")   # adding colums for excel comparison


######### Prepare tables for Mantle test

#function, which prepare the matrix suitable for the Mantle test

# FW_S_all_parasitoids
xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$S
xm[upper.tri(xm,diag=F)] <-FW_PC_res$S
FW_S_CP<-xm

# FW_OS_all_parasitoids
xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$OS
xm[upper.tri(xm,diag=F)] <-FW_PC_res$OS
FW_OS_CP<-xm



# FW_ST.lh_all_parasitoids

xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$ST.lh
xm[upper.tri(xm,diag=F)] <-FW_PC_res$ST.lh
FW_ST.lh_CP<-xm

# FW_WN_all_parasitoids

xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$WN
xm[upper.tri(xm,diag=F)] <-FW_PC_res$WN
FW_WN_CP<-xm

# FW_ST_all_parasitoids

xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$ST
xm[upper.tri(xm,diag=F)] <-FW_PC_res$ST
FW_ST_CP<-xm

# FW_ST.l_all_parasitoids

xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$ST.l
xm[upper.tri(xm,diag=F)] <-FW_PC_res$ST.l
FW_ST.l_CP<-xm

# FW_ST.h_all_parasitoids

xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$ST.h
xm[upper.tri(xm,diag=F)] <-FW_PC_res$ST.h
FW_ST.h_CP<-xm

#FW_ST.lh_all_parasitoids

xmnames<-unique(c(FW_PC_res$i,FW_PC_res$j))
xm <- matrix( rep(0,81),nrow=9,ncol=9)
colnames(xm) <- xmnames; rownames(xm) <- xmnames
xm[lower.tri(xm,diag=F)] <-FW_PC_res$ST.lh
xm[upper.tri(xm,diag=F)] <-FW_PC_res$ST.lh
FW_ST.lh_CP<-xm

sink("data/output/results.txt", append = T)
#Mantel test
FW_S_CP_mantel <-  mantel(FW_S_CP, Distance , method="spear")
FW_OS_CP_mantel <- mantel(FW_OS_CP, Distance , method="spear")
FW_WN_CP_mantel <- mantel(FW_WN_CP, Distance , method="spear")
FW_ST_CP_mantel <- mantel(FW_ST_CP, Distance , method="spear")
FW_ST.l_CP_mantel <- mantel(FW_ST.l_CP, Distance , method="spear")
FW_ST.h_CP_mantel <- mantel(FW_ST.h_CP, Distance , method="spear")
FW_ST.lh_CP_mantel <- mantel(FW_ST.lh_CP, Distance , method="spear")

FW_S_CP_mantel
FW_OS_CP_mantel
FW_WN_CP_mantel
FW_ST_CP_mantel
FW_ST.l_CP_mantel
FW_ST.h_CP_mantel
FW_ST.lh_CP_mantel
sink()


