# script used to analyse data reported in https://doi.org/10.1101/2022.08.25.505296
# data and analysis scripts can be found at www.github.com/BriVandrey/lec-to-mec-analysis
# to use, download files and set working directory (***) to relevant path for your machine on line 14 
# outputs: statistical results are saved as .txt files, plots are saved as .eps files
# Contact: bvandrey@ed.uk

# import libraries
library(rstatix)
library(tidyverse)
library(gridExtra)
library(ggpubr)

#set working directory - add relevant path for your machine
wd = "~/Desktop/analysis" 
#wd = " " # *** Set working directory here
setwd(wd)

#function to create and move into subdirectory
create_dir <- function(main_dir, sub_dir){
  new_dir <- file.path(main_dir, sub_dir) 
  dir.create(new_dir, showWarnings=FALSE)  #ignore if already exists
  setwd(new_dir)
}


#-----------------------------------------------------------------------------------------------------------
#PLOTTING UTILITY-------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#function for producing lineplot with individual observations - for pharmacology data
pharma_lp <- function(df, col, filter, xlabel1, xlabel2, ylabel, filename)
{
  # handle col variables
  col2 <- enquo(col) 
  
  #filter for cell type
  df2 <- filter(df, type==filter) 
  
  #Turn into a factor with the levels in the correct order
  df2$drug <- as.character(df2$drug)
  df2$drug <- factor(df2$drug, levels=unique(df2$drug))
  
  #summarise data 
  sum <-group_by(df2, drug, cell) %>% 
  summarise(mean_col = mean(!!col2))
  sum_mean_col <-group_by(sum, drug) %>% 
  get_summary_stats(type = "mean_se")

  #lineplot
  epsp_lp <- ggplot(sum_mean_col, aes(x=drug, y=mean))+
  theme_classic()+
  geom_line(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean_col, group = cell), size =0.6, colour="grey74")+ #individual traces
  geom_point(data = sum, inherit.aes = FALSE, aes(x = drug, y = mean_col, group = cell), size = 2, shape=111, colour="grey74") + #individual points
  geom_point(size = 1) + #mean points
  geom_line(aes(group=1), size = 0.6) + #mean trace
  geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width = 0.2, size = 0.6) + #error bars 
  scale_x_discrete(labels=c("Baseline", xlabel1, xlabel2)) +
  xlab("") + ylab(ylabel)+
  ylim(0,NA)
  epsp_lp + theme(axis.text.x = element_text(size = 10, angle=-50, color='black'), 
                  axis.text.y = element_text(size = 10, color='black'),
                  axis.title.y = element_text(size=12, margin= margin(t=0,r=0,b=0,l=0))) 

  ggsave(file=filename, width =35, height =45, units="mm") #save plot
}

#-----------------------------------------------------------------------------------------------------------
#ANALYSIS UTILITY-------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

# function to perform Friedmans repeated measures test on average data 
avg_fried <- function(df, cell_type, filename){
  
  df <- filter(df, type == cell_type) #filter by cell-type
  df<- as.matrix(df); df<- as.data.frame(df) # reformat data to correct type
  
  f_epsp <- friedman_test(data=df, epsp~drug|cell)
  eff_epsp<- friedman_effsize(data=df, epsp ~ drug|cell)
  f_ipsp <- friedman_test(data=df, ipsp_abs~drug|cell)
  eff_ipsp <- friedman_effsize(data=df, ipsp_abs ~ drug|cell)
  df$epsp <- as.numeric(as.character(df$epsp)); df$ipsp_abs <- as.numeric(as.character(df$ipsp_abs)); # re-make numeric for pairwise tests
  pw_epsp <- wilcox_test(data=df, epsp~drug, paired = TRUE, p.adjust.method = "bonferroni")
  pw_ipsp <- wilcox_test(data=df, ipsp_abs~drug, paired = TRUE, p.adjust.method = "bonferroni")
  
  df <- na.omit(df) # omit any empty rows for epsp halfwidth
  df<- as.matrix(df); df<- as.data.frame(df) #reformat again
  f_epsp_hw <- friedman_test(data=df, epsp_halfwidth~drug|cell)
  eff_epsp_hw<- friedman_effsize(data=df, epsp_halfwidth ~ drug|cell)
  df$epsp_halfwidth <- as.numeric(as.character(df$epsp_halfwidth))
  pw_epsp_hw <- wilcox_test(data=df, epsp_halfwidth~drug, paired = TRUE, p.adjust.method = "bonferroni")
  
  results <- list(f_epsp, eff_epsp, pw_epsp,
                  f_ipsp, eff_ipsp, pw_ipsp,
                  f_epsp_hw,eff_epsp_hw, pw_epsp_hw)
  capture.output(results, file=filename)
  
}

# function to perform Friedmans repeated measures test on sweep data - epsp & ipsp amplitude
cell_by_cell_fried <- function (df, cell_name, filename){
  
  df<- filter(df, cell==cell_name)
  
  #friedmans and posthoc tests
  f_epsp <- friedman_test(data=df, epsp~drug|pulse)
  pw_epsp <- wilcox_test(data=df, epsp~drug, paired = TRUE, p.adjust.method = "bonferroni")
  eff_epsp<- friedman_effsize(data=df, epsp ~ drug|pulse)
  f_ipsp <- friedman_test(data=df, ipsp_abs~drug|pulse)
  eff_ipsp <- friedman_effsize(data=df, ipsp_abs ~ drug|pulse)
  pw_ipsp <- wilcox_test(data=df, ipsp_abs~drug, paired = TRUE, p.adjust.method = "bonferroni")
  
  results <- list(f_epsp, eff_epsp, pw_epsp,
                  f_ipsp, eff_ipsp, pw_ipsp)
  capture.output(results, file=filename)
}

# function to perform Friedmans repeated measures test on sweep data - halfwidth
cell_by_cell_fried_hw <- function (df, cell_name, filename){
  
  df<- filter(df, cell==cell_name)
  df<- na.omit(df)
  df <- as.matrix(df); df <- as.data.frame(df)
  
  #friedmans and posthoc tests
  f_epsp_hw <- friedman_test(data=df, epsp_halfwidth~drug|pulse)
  eff_epsp_hw<- friedman_effsize(data=df, epsp_halfwidth ~ drug|pulse)
  df$epsp_halfwidth <- as.numeric(as.character(df$epsp_halfwidth)) # re-make numeric for pairwise tests
  pw_epsp_hw <- wilcox_test(data=df, epsp_halfwidth~drug, paired = TRUE, p.adjust.method = "bonferroni")
  
  results <- list(f_epsp_hw, eff_epsp_hw, pw_epsp_hw)
  capture.output(results, file=filename)
}

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------

# import data from .csv files and subset for plotting
glut_df <- read.csv(paste(wd, "data/glut_pharma.csv", sep="/"), head=TRUE, sep=",") #NBQX and APV - AMPA receptor antags
gaba_df <- read.csv(paste(wd, "data/gaba_pharma.csv", sep="/"), head=TRUE, sep=",") #Gabazine and CGP - GABA receptor antags
gaba_df_pulses<- read.csv(paste(wd, "data/gaba_pharma_individual_traces.csv", sep="/"), head=TRUE, sep=",") # individual traces - GABA receptor antags
ttx_df <- read.csv(paste(wd, "data/ttx_pharma.csv", sep="/"), head=TRUE, sep=",") #ttx and ap4 - monosynaptic expts

#add absolute values for ipsp
glut_df$ipsp_abs <- abs(glut_df$ipsp)
gaba_df$ipsp_abs <- abs(gaba_df$ipsp)
gaba_df_pulses$ipsp_abs <- abs(gaba_df_pulses$ipsp)
ttx_df$ipsp_abs <- abs(ttx_df$ipsp)

#subset based on order of drug application
cgp <- subset(gaba_df, pharma=="gaba_cgpgz")#CGP -> Gabazine
gz <- subset(gaba_df, pharma=="gaba_gzcgp")#Gabazine -> CGP

#create directory for outputs
create_dir(wd, "pharmacology")


#-----------------------------------------------------------------------------------------------------------
#FIGURE 5 --------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

# NBQX + APV - L2 SCs & L2 PCs - Panels A & B
pharma_lp(glut_df, epsp, "l2_sc", "NBQX", "APV", "EPSP (mV)",  "l2_sc_nbqx_apv_epsp.eps")
pharma_lp(glut_df, epsp, "l2_pc", "NBQX", "APV", "EPSP (mV)", "l2_pc_nbqx_apv_epsp.eps")
pharma_lp(glut_df, ipsp_abs, "l2_sc", "NBQX", "APV", "IPSP (mV)", "l2_sc_nbqx_apv_ipsp.eps")
pharma_lp(glut_df, ipsp_abs, "l2_pc", "NBQX", "APV", "IPSP (mV)", "l2_pc_nbqx_apv_ipsp.eps")

# TTX & 4AP - L2 SCs & L2 PCs - Panels E & F
pharma_lp(ttx_df, epsp, "l2_sc", "TTX", "4-AP", "EPSP (mV)", "l2_sc_ttx_4ap_epsp.eps")
pharma_lp(ttx_df, epsp, "l2_pc", "TTX", "4-AP", "EPSP (mV)", "l2_pc_ttx_4ap_epsp.eps")


#-----------------------------------------------------------------------------------------------------------
#FIGURE 6 --------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# statistics using averages from multiple traces------------------------------------------------------------
cgp_nas_removed <- na.omit(cgp) # manually remove NAs - required for Friedmans

# friedmans tests
# stellate cells - gabazine, CGP
avg_fried(gz, 'l2_sc', 'friedmans_l2_sc_gzcgp.txt')
avg_fried(gz, 'l2_pc', 'friedmans_l2_pc_gzcgp.txt')
avg_fried(cgp_nas_removed, 'l2_sc', 'friedmans_l2_sc_cgpgz.txt')

# statistics using individual traces----------------------------------------------------------------------------
# cell by cell analysis ----------------------------------------------------------------------------------------
create_dir(getwd(), "pulses_stats")

#Friedmans repeated measures
# gabazine -> cgp, stellates, cell by cell
cell_by_cell_fried(gaba_df_pulses, "1237_S2C2", "1237_S2C2_fried.txt")
cell_by_cell_fried(gaba_df_pulses, "1241_S2C2", "1241_S2C2_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1241_S2C2", "1241_S2C2_fried_halfwidth.txt")
cell_by_cell_fried(gaba_df_pulses, "1356_S5C1", "1356_S5C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1356_S5C1", "1356_S5C1_fried_halfwidth.txt")
cell_by_cell_fried(gaba_df_pulses, "1355_S2C1", "1355_S2C1_fried.txt")
cell_by_cell_fried(gaba_df_pulses, "1354_S4C1", "1354_S4C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1354_S4C1", "1354_S4C1_fried_halfwidth.txt")

# cgp -> gabazine, stellates, cell by cell
cell_by_cell_fried(gaba_df_pulses, "1309_S1C2", "1309_S1C2_fried.txt")
cell_by_cell_fried(gaba_df_pulses, "1309_S4C1", "1309_S4C1_fried.txt")
cell_by_cell_fried(gaba_df_pulses, "180321_S1C1", "180321_S1C1_fried.txt")
cell_by_cell_fried(gaba_df_pulses, "180321_S2C1", "180321_S2C1_fried.txt")

# gabazine -> cgp, pyramidal cells, cell by cell
cell_by_cell_fried(gaba_df_pulses, "1312_S2C1", "1312_S2C1_fried.txt")
cell_by_cell_fried(gaba_df_pulses, "1354_S3C1", "1354_S3C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1354_S3C1", "1354_S3C1_fried_halfwidth.txt")
cell_by_cell_fried(gaba_df_pulses, "1355_S1C1", "1355_S1C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1355_S1C1", "1355_S1C1_fried_halfwidth.txt")
cell_by_cell_fried(gaba_df_pulses, "1355_S4C1", "1355_S4C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1355_S4C1", "1355_S4C1_fried_halfwidth.txt")
cell_by_cell_fried(gaba_df_pulses, "1356_S1C1", "1356_S1C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1356_S1C1", "1356_S1C1_fried_halfwidth.txt")

# gabazine -> cgp, pyramidal cells, *** only one cell
cell_by_cell_fried(gaba_df_pulses, "1309_S2C1", "1309_S2C1_fried.txt")
cell_by_cell_fried_hw(gaba_df_pulses, "1309_S2C1", "1309_S2C1_fried_halfwidth.txt")


# --------------------------------------------------------------------------------------------------------------    
# nested ANOVAs with cell as a factor, includes all pulses instead of avg per cell -----------------------------
#nested_aov(gaba_df_pulses, 'l2_sc', 'gaba_gzcgp', "l2_sc_gzcgp_nested_aov.txt")
#nested_aov(gaba_df_pulses, 'l2_sc', 'gaba_cgpgz', "l2_sc_cgpgz_nested_aov.txt")
#nested_aov(gaba_df_pulses, 'l2_pc', 'gaba_gzcgp', "l2_pc_gzcgp_nested_aov.txt")

#--------------------------------------------------------------------------------------------------------------
setwd("..") # return to main pharmacology folder for plots

# lineplots
# Gabazine + CGP 55485 - L2 SCs & L2 PCs - Panels A & B
pharma_lp(gz, epsp, "l2_sc", "Gabazine", "CGP", "EPSP (mV)",  "l2_sc_gz_cgp_epsp.eps")
pharma_lp(gz, epsp, "l2_pc", "Gabazine", "CGP", "EPSP (mV)", "l2_pc_gz_cgp_epsp.eps")
pharma_lp(gz, epsp_halfwidth, "l2_sc", "Gabazine", "CGP", "EPSP Hafwidth (ms)",  "l2_sc_gz_cgp_epsp_hw.eps")
pharma_lp(gz, epsp_halfwidth, "l2_pc", "Gabazine", "CGP", "EPSP Hafwidth (ms)", "l2_pc_gz_cgp_epsp_hw.eps")
pharma_lp(gz, ipsp_abs, "l2_sc", "Gabazine", "CGP", "IPSP (mV)", "l2_sc_gz_cgp_ipsp.eps")
pharma_lp(gz, ipsp_abs, "l2_pc", "Gabazine", "CGP", "IPSP (mV)", "l2_pc_gz_cgp_ipsp.eps")

# CGP 55485  + Gabazine - L2 SCs & L2 PCs - Panels C & D
pharma_lp(cgp, epsp, "l2_sc", "CGP", "Gabazine", "EPSP (mV)",  "l2_sc_cgp_gz_epsp.eps")
pharma_lp(cgp, epsp, "l2_pc", "CGP", "Gabazine", "EPSP (mV)", "l2_pc_cgp_gz_epsp.eps")
pharma_lp(cgp, epsp_halfwidth, "l2_sc", "CGP", "Gabazine", "EPSP Hafwidth (ms)",  "l2_sc_cgp_gz_epsp_hw.eps")
pharma_lp(cgp, epsp_halfwidth, "l2_pc", "CGP", "Gabazine", "EPSP Hafwidth (ms)", "l2_pc_cgp_gz_epsp_hw.eps")
pharma_lp(cgp, ipsp_abs, "l2_sc", "CGP", "Gabazine", "IPSP (mV)", "l2_sc_cgp_gz_ipsp.eps")
pharma_lp(cgp, ipsp_abs, "l2_pc", "CGP", "Gabazine", "IPSP (mV)", "l2_pc_cgp_gz_ipsp.eps")


#-----------------------------------------------------------------------------------------------------------
#FIGURE 7 --------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

# NBQX + APV - L1 & L2 INS - Panel E
pharma_lp(glut_df, epsp, "l1_in", "NBQX", "APV", "EPSP (mV)",  "l1_in_nbqx_apv_epsp.eps")
pharma_lp(glut_df, epsp, "l2_in", "NBQX", "APV", "EPSP (mV)", "l2_in_nbqx_apv_epsp.eps")

# TTX & 4AP - L1 & L2 INS - Panel F
pharma_lp(ttx_df, epsp, "l1_in", "TTX", "4-AP", "EPSP (mV)", "l1_in_ttx_4ap_epsp.eps")
pharma_lp(ttx_df, epsp, "l2_in", "TTX", "4-AP", "EPSP (mV)", "l2_in_ttx_4ap_epsp.eps")
