# script used to analyse data reported in https://doi.org/10.1101/2022.08.25.505296
# data and analysis scripts can be found at www.github.com/BriVandrey/lec-to-mec-analysis
# to use, download files and set working directory (***) to relevant path for your machine on line 13 
# outputs: statistical results are saved as .txt files, plots are saved as .eps files

#import required libraries
library (rstatix)
library(dplyr)
library(tidyverse)
library(FSA)

#set working directory 
wd = "~/Documents/lecmec manuscript/analysis"  
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

lp_fluor <- function(df, exp, ymax){
  
  # subset normalised data for experiment type
  df2 <- filter(df, experiment==exp)
  df2 <- filter(df2, region== 'mec_l1_normalised' | region == 'dg_oml_normalised' )
  
  #turn into a factor with the levels in the correct order
  df2$region <- as.character(df2$region)
  df2$region<- factor(df2$region, levels=unique(df2$region))
  
  #summarise data
  sum <-group_by(df2, region, mouse) %>% 
    summarise(mean=mean(value))
  sum_mean<-group_by(sum,region) %>% 
    get_summary_stats(type = "mean_se")
  
  #lineplot 
  lp <- ggplot(sum_mean, aes(x=region, y=mean))+
    theme_classic()+
    geom_line(data = sum, inherit.aes = FALSE, aes(x = region, y = mean, group = mouse), size =0.6, colour="grey74")+ #individual traces
    geom_point(data = sum, inherit.aes = FALSE, aes(x = region, y = mean, group = mouse), size = 2, shape=111, colour="grey74") + #individual points
    geom_point(size = 1) + #mean points
    geom_line(aes(group=1), size = 0.6) + #mean trace
    geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width = 0.2, size = 0.6) + #error bars 
    scale_x_discrete(labels=c("MEC", "DG")) +
    xlab("") + ylab("Gray Value (pixels)") +
    ylim(0, ymax)
  lp + theme(axis.text.x = element_text(size = 10, color='black'), 
                  axis.text.y = element_text(size = 10, color='black'),
                  axis.title.y = element_text(size=12, margin= margin(t=0,r=0,b=0,l=0))) 
  
  ggsave(file=paste(exp, "fluorescence_intensity.eps", sep="_"), width =30, height =30, units="mm") #save plot
}


#function to plot fluorescence intensity x distance 
lp_fluor_dist <- function(df, label, xmax){
  
  #summarise data 
  sum <-group_by(df, distance_microns, mouse) %>% 
    summarise(mean=mean(grey_value_normalised))
  sum_mean<-group_by(sum,distance_microns) %>% 
    get_summary_stats(type = "mean_se")
  
  lp <- ggplot(df)+
    theme_classic()+
    geom_line(aes(x=distance_microns, y=grey_value_normalised, group=mouse), size=0.3, colour="grey74")+
    stat_smooth(aes(y=grey_value_normalised, x=distance_microns), method = lm, formula = y ~ poly(x, 10), se = FALSE, size=0.8, colour='black')+

    xlab("Distance") + ylab("Gray Value (pixels)") +
    ylim(0, 1)+
    xlim(0, xmax)
  lp + theme(axis.text.x = element_text(size = 10, color='black'), 
             axis.text.y = element_text(size = 10, color='black'),
             axis.title.y = element_text(size=12, margin= margin(t=0,r=0,b=0,l=0))) 

  ggsave(file=paste(label, "intensity_distance.eps", sep="_"), width =35, height =30, units="mm") #save plot
  
}


#function to plot boxplot of psp characteristics for neurons in multiple layers
ephys_bp <- function(df, x_col, y_col, y_label, filename){
  
  # handle col variables
  x_var <- enquo(x_col) 
  y_var <- enquo(y_col)
  
  # plot data
  bp<-ggplot(data=df, aes(!!x_var, !!y_var, fill=!!x_var)) +
    geom_boxplot(outlier.color=NA, show.legend=FALSE) + #remove outliers and legend
    geom_point(fill=NA, shape=1, size= 2, na.rm=TRUE) + #point fills
    labs(title="", x="", y=y_label) + #axis labels
    theme_classic()+
    scale_x_discrete(labels=c("L2 SC", "L2 PC", "L3 PC", "L5a PC", "L5b PC")) + #x tick labels
    scale_fill_brewer(palette="PRGn") #colour palette
  bp + theme(axis.text.x = element_text(size = 10, angle=-45, color='black'), 
             axis.text.y = element_text(size = 10, color='black'),
             axis.title.y = element_text(size=12, margin= margin(t=0,r=12,b=0,l=0)))
  
  ggsave(file=filename, width =48, height =60, units='mm') #save plot
}


# function to generate cumulative frequency plot, one measurement
cumulative_plot <- function(df, cell_type, x_col, x_label, y_label, filename){
  
  x_col <- enquo(x_col) # handle col variable
  df2 <- subset(df, type == cell_type) #subset df by cell type
  
  cum_plt <- ggplot(df2, aes(!!x_col)) + stat_ecdf(geom = "step", pad = FALSE) +
    labs(title="", x=x_label, y=y_label) + #axis labels
    theme_classic()
  cum_plt + theme(axis.text.x = element_text(size = 9, angle=0, color='black'), 
                  axis.text.y = element_text(size = 9, color='black'),
                  axis.title.y = element_text(size=9, margin= margin(t=0,r=8,b=0,l=0)),
                  axis.title.x = element_text(size=9, margin= margin(t=8,r=0,b=0,l=0)))
  ggsave(file=filename, width =40, height =40, units='mm') #save plot 
}


# function to generate cumulative frequency plot with two measurements overlaid
cumulative_plot_overlay <- function(df, x_col, y_col, x_label, filename){
  
  x_col <- enquo(x_col) # handle col variables
  y_col <- enquo(y_col)
  
  cum_plt <- ggplot(df, aes(x=(!!x_col), color=(!!y_col))) + geom_step(aes(y=..y..),stat="ecdf") +
    labs(title="", x=x_label, y="") + #axis labels
    xlim(0,1)+
    theme_classic()
  
  cum_plt + theme(axis.text.x = element_text(size = 9, angle=0, color='black'), 
                  axis.text.y = element_text(size = 9, color='black'),
                  axis.title.y = element_text(size=9, margin= margin(t=0,r=8,b=0,l=0)),
                  axis.title.x = element_text(size=9, margin= margin(t=8,r=0,b=0,l=0)),
                  legend.position='none')
  ggsave(file=filename, width =45, height =50, units='mm') #save plot 
}
# function to plot boxplot for two groups
bp_two_groups <- function(df, x_col, y_col, x_label1, x_label2, y_label, filename){
  
  #handle col variables
  x_var <- enquo(x_col) 
  y_var <- enquo(y_col)
  
  bp <- ggplot(data=df, aes(!!x_var, !!y_var, fill=!!x_var)) +
    geom_boxplot(outlier.color=NA, show.legend=FALSE) + #remove outliers and legend
    geom_point(fill=NA, shape=1, size= 2, na.rm=TRUE) + #point fills
    labs(title="", x="", y=y_label) + #axis labels
    theme_classic()+
    scale_x_discrete(labels=c(x_label1, x_label2)) + #x tick labels
    scale_fill_brewer(palette="PRGn") #colour palette
  bp + theme(axis.text.x = element_text(size = 10, angle=-45, color='black'), 
             axis.text.y = element_text(size = 10, color='black'),
             axis.title.y = element_text(size=12, margin= margin(t=0,r=12,b=0,l=0)))
  
  ggsave(file=filename, width =40, height =60, units='mm') #save plot
}


# function to generate scatterplot with simple linear model
corr_plot <- function(df, cell_type, x_col, y_col, x_label, y_label, filename){
  
  # handle col variables
  x_var <- enquo(x_col)
  y_var <- enquo(y_col)
  
  df <- filter(df, type==cell_type)  #filter by cell type
  
  # plot data - option to colour code by cell type, uncommented
  plt <- ggplot(df, aes(x=!!x_var, y=!!y_var)) +
    geom_point(fill=NA, shape=1, size=2, na.rm=TRUE, show.legend=FALSE)+ #aes(colour=factor(type)) 
    labs(title="", x=x_label, y=y_label) +
    scale_color_brewer(palette="BuPu")+  
    #xlim(0, 25)+ #limits set for figures used in paper
    #ylim(0,1.05)+ 
    theme_classic() 
  
  plt + theme(axis.text.x = element_text(size = 12, angle=-40, color='black'), 
              axis.text.y = element_text(size = 12, color='black'),
              axis.title.y = element_text(size=12, margin= margin(t=0,r=14,b=0,l=0)), 
              axis.title.x = element_text(size=12, margin= margin(t=14,r=0,b=0,l=0)),
              legend.position = 'none') +
    
    geom_smooth(method='glm') 
  
  ggsave(file=filename, width =40, height =45, units='mm') #save plot
}


# function to create line plot - requires frequency input to filter (10 or 20 Hz)
# automatically labels 1, 10 on x axis so need to change manually for Figure 10 - Supplement 2
lp_hf <- function(df, col, type_filter, frequency_filter, pulse_to_remove, y_label, name){
  
  #filter for cell type
  df2 <- filter(df, type==type_filter & frequency==frequency_filter & pulse!=pulse_to_remove) 
  
  # handle col variable
  col <- enquo(col)
  
  #turn into a factor with the levels in the correct order
  df2$pulse <- as.character(df2$pulse)
  df2$pulse<- factor(df2$pulse, levels=unique(df2$pulse))
  
  #summarise data - epsp width
  sum <-group_by(df2, pulse, cell) %>% 
    summarise(mean=mean(!!col))
  sum_mean<-group_by(sum, pulse) %>% 
    get_summary_stats(type = "mean_se")
  
  #lineplot -epsp width
  epsp_lp <- ggplot(sum_mean, aes(x=pulse, y=mean))+
    theme_classic()+
    geom_line(data = sum, inherit.aes = FALSE, aes(x = pulse, y = mean, group = cell), size =0.6, colour="grey74")+ #individual traces
    geom_point(data = sum, inherit.aes = FALSE, aes(x = pulse, y = mean, group = cell), size = 2, shape=111, colour="grey74") + #individual points
    geom_point(size = 1) + #mean points
    geom_line(aes(group=1), size = 0.6) + #mean trace
    geom_errorbar(aes(ymax = mean+se, ymin = mean-se), width = 0.2, size = 0.6) + #error bars 
    scale_x_discrete(labels=c("1", "10")) +
    xlab("Pulse") + ylab(y_label) +
    ylim(0,NA)
  epsp_lp + theme(axis.text.x = element_text(size = 10, color='black'), 
                  axis.text.y = element_text(size = 10, color='black'),
                  axis.title.y = element_text(size=12, margin= margin(t=0,r=0,b=0,l=0))) 
  
  ggsave(file=paste(type_filter, frequency_filter, name, "high_frequency_stim.eps", sep="_"), width =30, height =30, units="mm") #save plot
}


#-----------------------------------------------------------------------------------------------------------
#FIGURE 1---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#analyse puncta density across layers & mediolaterial axis

#import data
df_syn_layers <- read.csv(file= paste(wd, "data/syn_layers.csv", sep="/"), head=TRUE, sep=",") 
df_syn_medlat <- read.csv(file= paste(wd, "data/syn_medlat.csv", sep="/"), head=TRUE, sep=",") 

#create directory for outputs
create_dir(wd, "puncta")

#summary statistics 
df.summary_layers <- df_syn_layers %>% 
  group_by(layer) %>% 
  summarise(count = n(), 
            sd = sd(density, na.rm = TRUE),
            se = sd/(sqrt(count)),
            density = mean(density))

df.summary_medlat <- df_syn_medlat %>% 
  group_by(distance) %>% 
  summarise(count = n (),
            sd = sd(density, na.rm = TRUE),
            se = sd/(sqrt(count)),
            density = mean(density))

#friedmans test, medlat density
fried_medlat <- df_syn_medlat %>% friedman_test(density~distance|id)
eff_size_medlat<- df_syn_medlat %>% friedman_effsize(density~distance|id)
pw_fried_medlat <- df_syn_medlat %>% #pairwise comparisons
  wilcox_test(density ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

#friedmans test, laminar density
fried_layers <- df_syn_layers %>% friedman_test(density~layer|id)
eff_size_layers<- df_syn_layers %>% friedman_effsize(density~layer|id)
pw_fried_layers <- df_syn_layers %>% #pairwise comparisons
  wilcox_test(density ~ layer, paired = TRUE, p.adjust.method = "bonferroni")

#line plots of puncta data - Figure 1, Panel D
#line plot for layers
layers <- ggplot(df_syn_layers, aes(x=layer, y=density)) +
  geom_jitter(position = position_jitter(0.2), size=1, color = "darkgray") + #jitterplot
  geom_line(aes(group = 1), data = df.summary_layers, size=0.6) + #line plot
  geom_errorbar(aes(ymin = density-se, ymax = density+se),  #error bar
  data = df.summary_layers, width = 0.2, size=0.5) +
  geom_point(data = df.summary_layers, size = 1) +
  ylim(0,10000)+
  scale_x_discrete(labels=c("L1", "L2", "L3", "L5a", "L5b")) + #x tick labels
  labs(title="", x="", y="Density (mm )") + #axis labels
  theme_classic()

  layers + theme(axis.text.x = element_text(size = 10, angle=-40, color='black'), 
            axis.text.y = element_text(size =10, color='black'),
            axis.title.y = element_text(size=12, margin= margin(t=0,r=0,b=0,l=0)),
            legend.title = element_blank()) 
  
  ggsave(filename="syn_layers.eps", width=39, height=45, units="mm") #save plot

#line plots for mediolateral distance
medlat <- ggplot(df_syn_medlat, aes(x=distance, y=density)) +
    geom_jitter(position = position_jitter(0.2), size=1, color = "darkgray") + #jitterplot
    geom_line(aes(group = 1), data = df.summary_medlat, size=0.6) + #line plot
    geom_errorbar(aes(ymin = density-se, ymax = density+se), #error bar
    data = df.summary_medlat, width = 0.2, size=0.5) +
    geom_point(data = df.summary_medlat, size = 1) +
    ylim(0,10000)+
    scale_x_discrete(labels=c("150", "300", "450", "600", "750", "900")) + #x tick labels
    labs(title="", x="", y="Density (mm )") + #axis labels
    theme_classic()
  
  medlat + theme(axis.text.x = element_text(size = 10, angle=-40, color='black'), 
                 axis.text.y = element_text(size = 10, color='black'),
                 axis.title.y = element_text(size=12, margin= margin(t=0,r=0,b=0,l=0)),
                 legend.title = element_blank()) 
  
  ggsave(filename= 'syn_medlat.eps', width = 39, height= 45, units="mm") #save plot
  
  #save analysis outputs to text file
  results <- list(df.summary_layers, df.summary_medlat, 
                  fried_layers, eff_size_layers, pw_fried_layers,
                  fried_medlat, eff_size_medlat, pw_fried_medlat)
  capture.output(results, file='puncta_stats.txt')

  
#-----------------------------------------------------------------------------------------------------------  
#FIGURE 2---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

setwd(wd)
#create directory for outputs
create_dir(wd, "fluorescence_intensity")

# import data
fluor_mec_dg <- read.csv(file= paste(wd, "data/fluorescence_intensity/mec_dg.csv", sep="/"), head=TRUE, sep=",") 
fluor_dg <- read.csv(file= paste(wd, "data/fluorescence_intensity/retrograde_dg.csv", sep="/"), head=TRUE, sep=",") 
fluor_mec <- read.csv(file= paste(wd, "data/fluorescence_intensity/retrograde_mec.csv", sep="/"), head=TRUE, sep=",") 
distance_dg <- read.csv(file= paste(wd, "data/fluorescence_intensity/distance_dg.csv", sep="/"), head=TRUE, sep=",") 
distance_mec <- read.csv(file= paste(wd, "data/fluorescence_intensity/distance_mec.csv", sep="/"), head=TRUE, sep=",") 

normalised = filter(fluor_mec_dg, region == 'mec_l1_normalised'| region == 'dg_oml_normalised') # normalised data only
df_retro = filter(normalised, experiment == 'retro_cre')

#summary statistics - mec vs dg
df.retro_cre <- df_retro %>% 
  group_by(region) %>% 
  summarise(count = n(), 
            sd = sd(value, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(value))

#compare dg & mec
retro_cre_wilcox <- wilcox.test(value~region, data=df_retro) 

#summary statistics - intensity x distance
df.summary_distance_dg <- distance_dg %>% 
  group_by(distance) %>% 
  summarise(count = n(), 
            sd = sd(average, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(average))

df.summary_distance_mec <- distance_mec %>% 
  group_by(distance) %>% 
  summarise(count = n(), 
            sd = sd(average, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(average))

#friedmans test, intensity x distance, 50 uM steps (averaged) - dg
fried_dg <- distance_dg %>% friedman_test(average~distance|mouse)
eff_size_dg <- distance_dg %>% friedman_effsize(average~distance|mouse)
pw_dg <- distance_dg %>% #pairwise comparisons
  wilcox_test(average ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

#friedmans test, intensity x distance, 50 uM steps (averaged) - mec
fried_mec <- distance_mec %>% friedman_test(average~distance|mouse)
eff_size_mec <- distance_mec %>% friedman_effsize(average~distance|mouse)
pw_mec <- distance_mec %>% #pairwise comparisons
  wilcox_test(average ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

# line plot comparing intensity in DG & MEC - Figure 2, Panel D
lp_fluor(fluor_mec_dg, 'retro_cre', 40)

# line plot of intensity x distance - Figure 2, Panel E
lp_fluor_dist(fluor_dg, 'retro_dg', 200)
lp_fluor_dist(fluor_mec, 'retro_mec', 400)

#save analysis outputs to text file
results <- list(df.retro_cre, retro_cre_wilcox, 
                df.summary_distance_dg, fried_dg, eff_size_dg, pw_dg,
                df.summary_distance_mec, fried_mec, eff_size_mec, pw_mec)
capture.output(results, file='retro_cre_intensity.txt')


#-----------------------------------------------------------------------------------------------------------  
#FIGURE 3---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# import data
fluor_dg_pir <- read.csv(file= paste(wd, "data/fluorescence_intensity/anterograde_pir_dg.csv", sep="/"), head=TRUE, sep=",") 
fluor_mec_pir <- read.csv(file= paste(wd, "data/fluorescence_intensity/anterograde_pir_mec.csv", sep="/"), head=TRUE, sep=",") 
fluor_dg_mpfc <- read.csv(file= paste(wd, "data/fluorescence_intensity/anterograde_mpfc_dg.csv", sep="/"), head=TRUE, sep=",") 
fluor_mec_mpfc <- read.csv(file= paste(wd, "data/fluorescence_intensity/anterograde_mpfc_mec.csv", sep="/"), head=TRUE, sep=",") 

distance_dg_pir <- read.csv(file= paste(wd, "data/fluorescence_intensity/distance_dg_pir.csv", sep="/"), head=TRUE, sep=",") 
distance_mec_pir <- read.csv(file= paste(wd, "data/fluorescence_intensity/distance_mec_pir.csv", sep="/"), head=TRUE, sep=",") 
distance_dg_mpfc <- read.csv(file= paste(wd, "data/fluorescence_intensity/distance_dg_mpfc.csv", sep="/"), head=TRUE, sep=",") 
distance_mec_mpfc <- read.csv(file= paste(wd, "data/fluorescence_intensity/distance_mec_mpfc.csv", sep="/"), head=TRUE, sep=",") 

df_pir <- filter(normalised, experiment == 'anter_cre_pir')
df_mpfc <- filter(normalised, experiment == 'anter_cre_mpfc')

#summary statistics - mec vs dg, piriform cortex
df.antero_pir <- df_pir %>% 
  group_by(region) %>% 
  summarise(count = n(), 
            sd = sd(value, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(value))

#compare dg & mec, piriform
anter_cre_wilcox_pir <- wilcox.test(value~region, data=df_pir)

#summary statistics - mec vs dg, medial prefrontal cortex
df.antero_mpfc <- df_mpfc %>% 
  group_by(region) %>% 
  summarise(count = n(), 
            sd = sd(value, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(value))

#compare dg & mec, medial prefrontal
anter_cre_wilcox_mpfc <- wilcox.test(value~region, data=df_mpfc)

#summary statistics - intensity x distance, piriform
df.summary_distance_dg_pir <- distance_dg_pir %>% 
  group_by(distance) %>% 
  summarise(count = n(), 
            sd = sd(average, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(average))

df.summary_distance_mec_pir <- distance_mec_pir %>% 
  group_by(distance) %>% 
  summarise(count = n(), 
            sd = sd(average, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(average))

#friedmans test, intensity x distance, 50 uM steps (averaged) - dg
fried_dg_pir <- distance_dg_pir %>% friedman_test(average~distance|mouse)
eff_size_dg_pir <- distance_dg_pir %>% friedman_effsize(average~distance|mouse)
pw_dg_pir <- distance_dg_pir %>% #pairwise comparisons
  wilcox_test(average ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

#friedmans test, intensity x distance, 50 uM steps (averaged) - mec
fried_mec_pir <- distance_mec_pir %>% friedman_test(average~distance|mouse)
eff_size_mec_pir <- distance_mec_pir %>% friedman_effsize(average~distance|mouse)
pw_mec_pir <- distance_mec %>% #pairwise comparisons
  wilcox_test(average ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

#summary statistics - intensity x distance, medial prefrontal cortex
df.summary_distance_dg_mpfc <- distance_dg_mpfc %>% 
  group_by(distance) %>% 
  summarise(count = n(), 
            sd = sd(average, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(average))

df.summary_distance_mec_mpfc <- distance_mec_mpfc %>% 
  group_by(distance) %>% 
  summarise(count = n(), 
            sd = sd(average, na.rm = TRUE),
            se = sd/(sqrt(count)),
            average = mean(average))

#friedmans test, intensity x distance, 50 uM steps (averaged) - dg
fried_dg_mpfc <- distance_dg_mpfc %>% friedman_test(average~distance|mouse)
eff_size_dg_mpfc <- distance_dg_mpfc %>% friedman_effsize(average~distance|mouse)
pw_dg_mpfc <- distance_dg_mpfc %>% #pairwise comparisons
  wilcox_test(average ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

#friedmans test, intensity x distance, 50 uM steps (averaged) - mec
fried_mec_mpfc <- distance_mec_mpfc %>% friedman_test(average~distance|mouse)
eff_size_mec_mpfc <- distance_mec_mpfc %>% friedman_effsize(average~distance|mouse)
pw_mec_mpfc <- distance_mec %>% #pairwise comparisons
  wilcox_test(average ~ distance, paired = TRUE, p.adjust.method = "bonferroni")

# line plots comparing intensity in DG & MEC - Figure 3, Panels C, G
lp_fluor(fluor_mec_dg, 'anter_cre_pir', 150)  
lp_fluor(fluor_mec_dg, 'anter_cre_mpfc', 40)  

# line plots of intensity x distance - Figure 3, Panels G, H
lp_fluor_dist(fluor_dg_pir, 'antero_pir_dg', 200)
lp_fluor_dist(fluor_mec_pir, 'antero_pir_mec', 400)
lp_fluor_dist(fluor_dg_mpfc, 'antero_mpfc_dg', 200)
lp_fluor_dist(fluor_mec_mpfc, 'antero_mpfc_mec', 400)

#save analysis outputs to text file
results <- list(df.antero_pir, anter_cre_wilcox_pir, 
                df.antero_mpfc, anter_cre_wilcox_mpfc,
                df.summary_distance_dg_pir, fried_dg_pir, eff_size_dg_pir, pw_dg_pir,
                df.summary_distance_mec_pir, fried_mec_pir, eff_size_mec_pir, pw_mec_pir,
                df.summary_distance_dg_mpfc, fried_dg_mpfc, eff_size_dg_mpfc, pw_dg_mpfc,
                df.summary_distance_mec_mpfc, fried_mec_mpfc, eff_size_mec_mpfc, pw_mec_mpfc)
capture.output(results, file='antero_cre_intensity.txt')


#-----------------------------------------------------------------------------------------------------------  
#FIGURE 4---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#comparison of response amplitudes and latency across L1-L5b

#create new directory for outputs
setwd(wd)
create_dir(wd, "psp_characteristics")

#import data
df <- read.csv(file= paste(wd, "data/main.csv", sep="/"), head=TRUE, sep=",") 
  
#subset electrophysiology data 
df <- filter(df, response=='Y') #filter out non-responsive cells
df_pc <- subset(df, type == 'l2_sc' | type == 'l2_pc' | type == 'l3_pc' | type == 'l5a_pc'|type == 'l5b_pc') #principal neurons
df_l2 <- subset(df, type == 'l2_sc' | type == 'l2_pc')

# summary statistics
sum_stats_epsp <- df_pc %>% #summary statistics - epsp amplitudes
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(epsp, na.rm=TRUE), 
            sd = sd(epsp, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp, na.rm=TRUE),
            iqr = IQR(epsp, na.rm=TRUE))

#convert ipsps to absolute values
df_pc$ipsp_abs <- abs(df_pc$ipsp) #average values
df_pc$ipsp_pk_abs <- abs(df_pc$ipsp_pk) #peak values

sum_stats_ipsp <- df_pc %>% #summary statistics - ipsp amplitudes
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(ipsp_abs, na.rm=TRUE), 
            sd = sd(ipsp_abs, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ipsp_abs, na.rm=TRUE),
            iqr = IQR(ipsp_abs, na.rm=TRUE))

sum_stats_epsp_latency <- df_pc %>% #summary statistics - epsp latency
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(epsp_latency, na.rm=TRUE), 
            sd = sd(epsp_latency, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp_latency, na.rm=TRUE),
            iqr = IQR(epsp_latency, na.rm=TRUE))

sum_stats_ipsp_latency <- df_pc %>% #summary statistics - ipsp latency
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(ipsp_latency, na.rm=TRUE), 
            sd = sd(ipsp_latency, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ipsp_latency, na.rm=TRUE),
            iqr = IQR(ipsp_latency, na.rm=TRUE))

#comparative statistics - kruskal-wallis with multiple pairwise comparisons
epsp_kw <- kruskal.test(epsp~type, data=df_pc) # epsp x cell type
epsp_pw_dunn = dunnTest(epsp~type, data=df_pc, method="bonferroni") # full results in Figure 4, Supplement 2
ipsp_kw <- kruskal.test(ipsp_abs~type, data=df_pc) # ipsp x cell type
epsp_latency_kw <- kruskal.test(epsp_latency~type, data=df_pc) # epsp latency x cell type
ipsp_latency_kw <- kruskal.test(ipsp_latency~type, data=df_pc) # ipsp latency x cell type

# reorder type factor for plotting
df_pc$type <- factor(df_pc$type, levels=c("l2_sc", "l2_pc", "l3_pc", "l5a_pc", "l5b_pc")) 

# generate boxplots for Figure 4, panels E & F
ephys_bp(df_pc, type, epsp, "EPSP (mV)", "epsp_pc.eps") 
ephys_bp(df_pc, type, ipsp_abs, "IPSP (mV)", "ipsp_pc.eps") 
ephys_bp(df_pc, type, epsp_latency, "EPSP Latency (ms)", "epsp_latency_pc.eps") 
ephys_bp(df_pc, type, ipsp_latency, "IPSP Latency (ms)", "ipsp_latency_pc.eps") 

# for ww Neurise talk
bp_two_groups(df_l2, 'type', 'epsp', "L2 SC", "L2 PC", "EPSP (mV)", "epsp_wwneurise.eps")
bp_two_groups(df_l2, 'type', 'ipsp_abs', "L2 SC", "L2 PC", "IPSP (mV)", "ipsp_wwneurise.eps")

#-----------------------------------------------------------------------------------------------------------
#FIGURE 5---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#cumulative frequency plots for sd of latency

# summary statistics
sum_stats_epsp_latency_sd <- df_pc %>% #summary statistics - epsp latency sd
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(epsp_latency, na.rm=TRUE), 
            sd = sd(epsp_latency, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp_latency, na.rm=TRUE),
            iqr = IQR(epsp_latency, na.rm=TRUE))

#generate cum. frequency plots for Figure 5, panels C & D
cumulative_plot(df_pc, "l2_sc", epsp_latency_sd, "SD Latency (ms)", "Cum. Probability", "cum_freq_latency_sd_l2_sc.eps")
cumulative_plot(df_pc, "l2_pc", epsp_latency_sd, "SD Latency (ms)", "Cum. Probability", "cum_freq_latency_sd_l2_pc.eps")

# save results
results <- list(sum_stats_epsp, sum_stats_ipsp, sum_stats_epsp_latency, sum_stats_ipsp_latency,
                epsp_kw, epsp_pw_dunn,
                ipsp_kw, epsp_latency_kw, ipsp_latency_kw,
                sum_stats_epsp_latency_sd)
capture.output(results, file="psp_characteristics_stats_principal_cells.txt")


#-----------------------------------------------------------------------------------------------------------
#FIGURE 7---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#comparison of epsp amplitude and latency across L1 & L2 IN

df_in <- subset(df, type == 'l1_in' | type == 'l2_in') #interneurons

sum_stats_epsp_in <- df_in %>% # summary statistics - epsp amplitude
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(epsp, na.rm=TRUE), 
            sd = sd(epsp, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp_latency, na.rm=TRUE),
            iqr = IQR(epsp_latency, na.rm=TRUE))

sum_stats_epsp_latency_in <- df_in %>% #summary statistics - epsp latency
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(epsp_latency, na.rm=TRUE), 
            sd = sd(epsp_latency, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp_latency, na.rm=TRUE),
            iqr = IQR(epsp_latency, na.rm=TRUE))

#compare epsp amplitude and latency 
epsp_wilcox <- wilcox.test(epsp~type, data=df_in)
epsp_latency_wilcox <- wilcox.test(epsp_latency~type, data=df_in)

# generate boxplots for Figure 7, panel B
df_in$type <- factor(df_in$type, levels=c("l1_in", "l2_in")) #reorder type factor for plotting
bp_two_groups(df_in, type, epsp, "L1 IN", "L2 IN", y_label="EPSP (mV)", "epsp_in.eps")
bp_two_groups(df_in, type, epsp_latency, "L1 IN", "L2 IN", "EPSP Latency (ms)", "epsp_latency_in.eps")

# generate cumulative frequency plots --- not included in manuscript
cumulative_plot(df_in, "l1_in", epsp_latency_sd, "SD Latency (ms)", "Cum. Probability", "cum_freq_latency_sd_l1_in.eps")
cumulative_plot(df_in, "l2_in", epsp_latency_sd, "SD Latency (ms)", "Cum. Probability", "cum_freq_latency_sd_l2_in.eps")

# save results
results <- list(sum_stats_epsp_in, sum_stats_epsp_latency_in,
                epsp_wilcox, epsp_latency_wilcox)
capture.output(results, file="psp_characteristics_stats_interneurons.txt")


#-----------------------------------------------------------------------------------------------------------
#FIGURE 8---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# paired recordings between layer 1 interneurons and L2 stellate/pyramidal cells

#create new directory
setwd(wd)
create_dir(wd, "paired_recordings")

#import data
df_paired <- read.csv(file = paste(wd, "data/paired.csv", sep="/"), head=TRUE, sep=",") 

#calculate absolute values for ipsp
df_paired$abs_amplitude <- abs(df_paired$amplitude)
df_paired$abs_pk_amplitude <- abs(df_paired$pk_amplitude)

#convert to factor and order by pulses
df_paired$type <- as.character(df_paired$type)
df_paired$type <- factor(df_paired$type, levels=unique(df_paired$type))

sum_stats_paired_amplitude <- df_paired %>% #summary statistics - ipsp amplitude
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(abs_amplitude, na.rm=TRUE), 
            sd = sd(abs_amplitude, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(abs_amplitude, na.rm=TRUE),
            iqr = IQR(abs_amplitude, na.rm=TRUE))

sum_stats_paired_width <- df_paired %>% #summary statistics - ipsp amplitude
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(ipsp_width_from_avg, na.rm=TRUE), 
            sd = sd(ipsp_width_from_avg, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ipsp_width_from_avg, na.rm=TRUE),
            iqr = IQR(ipsp_width_from_avg, na.rm=TRUE))

# generate boxplots for Figure 8, panel D
bp_two_groups(df_paired, type, abs_amplitude, "L2 SC", "L2 PC", "IPSP (mV)", "ipsp_amplitude.eps") 
bp_two_groups(df_paired, type, ipsp_width_from_avg, "L2 SC", "L2 PC", "Halfwidth (ms)", "ipsp_halfwidth.eps") 

# save results
results <- list(sum_stats_paired_amplitude, sum_stats_paired_width)
capture.output(results, file="paired_recording_stats.txt")


#-----------------------------------------------------------------------------------------------------------
#FIGURE 9---------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#compare relative strength of excitation and inhibition 

setwd(wd)
create_dir(wd, "ei")

#summary statistics
sum_stats_ei <- df_pc %>%
  group_by(type) %>% 
  summarise(count = n(),
            mean = mean(ei, na.rm=TRUE), 
            sd = sd(ei, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ei, na.rm=TRUE),
            iqr = IQR(ei, na.rm=TRUE))

# subset l2 stellates and pyramidal cells
df_l2 <- subset(df_pc, type=='l2_sc' | type =='l2_pc')

#log transform e/ipsp
df_l2$log_epsp <- log(df_l2$epsp)
df_l2$log_ipsp <- log(df_l2$ipsp_abs)

#subset again into sc/pc
df_l2_sc <- subset(df_l2, type=='l2_sc')
df_l2_pc <- subset(df_l2, type=='l2_pc')

#compare ei ratio - Wilcoxon/ Mann Whitney U
ei_wilcox <- wilcox.test(ei~type, data=df_l2)

# generate boxplot for Figure 9, Panel A
bp_two_groups(df_l2, type, ei, 'L2 SC', 'L2 PC', 'E-I Ratio', 'ei.eps')
bp_two_groups(df_l2, type, epsp, 'L2 SC', 'L2 PC', 'EPSP (mV)', 'epsp_neurise.eps')
bp_two_groups(df_l2, type, ipsp_abs, 'L2 SC', 'L2 PC', 'IPSP (mV)', 'ipsp_neurise.eps')

#correlational stats - excitation vs inhibition
ei_sc_log <- summary(lm(df_l2_sc$log_epsp~df_l2_sc$log_ipsp)) #log transformed
ei_pc_log <- summary(lm(df_l2_pc$log_epsp~df_l2_pc$log_ipsp))

# plot excitation v inhibition, Figure 9, Panel B
corr_plot(df_l2, 'l2_sc', log_ipsp, log_epsp, "Log (EPSP)", "Log (IPSP)", 'ei_corr_log_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', log_ipsp, log_epsp, "Log (EPSP)", "Log (IPSP)", 'ei_corr_log_l2_pc.eps')

#correlational stats - ei ratio vs epsp & position
ei_epsp_sc <- summary(lm(df_l2_sc$epsp~df_l2_sc$ei)) # epsp
ei_epsp_pc <- summary(lm(df_l2_pc$epsp~df_l2_pc$ei))
ei_medlat_sc <- summary(lm(df_l2_sc$dist_psub~df_l2_sc$ei)) # medlat distance
ei_medlat_pc <- summary(lm(df_l2_pc$dist_psub~df_l2_pc$ei)) 
ei_dv_sc <- summary(lm(df_l2_sc$dv~df_l2_sc$ei)) # dv 
ei_dv_pc <- summary(lm(df_l2_pc$dv~df_l2_pc$ei)) 

# plot e-i against EPSP amplitude (Panel C), dorsoventral & medio-lateral position (Panels E-F)
corr_plot(df_l2, 'l2_sc', epsp, ei, x_label=Delta~EPSP, y_label="E-I Ratio", 'ei_epsp_l2_sc.eps') # Panel C
corr_plot(df_l2, 'l2_pc', epsp, ei, x_label=Delta~EPSP, y_label="E-I Ratio", 'ei_epsp_l2_pc.eps') 
corr_plot(df_l2, 'l2_sc', dist_psub, ei, 'Distance (mm)', 'E-I Ratio', 'ei_medlat_l2_sc.eps') # Panel E
corr_plot(df_l2, 'l2_pc', dist_psub, ei, x_label='Distance (mm)', y_label='E-I Ratio', 'ei_medlat_l2_pc.eps')
corr_plot(df_l2, 'l2_sc', dv, ei, x_label='DV (mm)', y_label='E-I Ratio', 'ei_dv_l2_sc.eps') #Panel F
corr_plot(df_l2, 'l2_pc', dv, ei, x_label='DV (mm)', y_label='E-I Ratio', 'ei_dv_l2_pc.eps') 

# save results
results <- list(sum_stats_ei, ei_wilcox, 
                ei_sc_log, ei_pc_log, 
                ei_epsp_sc, ei_epsp_pc,
                ei_medlat_sc, ei_medlat_pc, ei_dv_sc, ei_dv_pc)
capture.output(results, file="ei_stats.txt")


#-----------------------------------------------------------------------------------------------------------
#Figure 9, Supplement 1-------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#create directory for supplemental plots/data
dir <- getwd()
create_dir(dir, 'Fig9_supplemental')

#compare ei ratio x resting membrane potential 
ei_rmp_sc <- summary(lm(df_l2_sc$vm~df_l2_sc$ei))
ei_rmp_pc <- summary(lm(df_l2_pc$vm~df_l2_pc$ei))

#compare excitation vs inhibition
ei_sc <- summary(lm(df_l2_sc$epsp~df_l2_sc$ipsp_abs)) #raw data
ei_pc <- summary(lm(df_l2_pc$epsp~df_l2_pc$ipsp_abs))

#compare e/ipsp amplitude x neuron position - principal neurons 
epsp_medlat_sc <- summary(lm(df_l2_sc$dist_psub~df_l2_sc$epsp)) 
epsp_medlat_pc <- summary(lm(df_l2_pc$dist_psub~df_l2_pc$epsp)) 
epsp_dv_sc <- summary(lm(df_l2_sc$dv~df_l2_sc$epsp)) 
epsp_dv_pc <- summary(lm(df_l2_pc$dv~df_l2_pc$epsp))
ipsp_medlat_sc <- summary(lm(df_l2_sc$dist_psub~df_l2_sc$ipsp_abs)) 
ipsp_medlat_pc <- summary(lm(df_l2_pc$dist_psub~df_l2_pc$ipsp_abs)) 
ipsp_dv_sc <- summary(lm(df_l2_sc$dv~df_l2_sc$ipsp_abs)) 
ipsp_dv_pc <- summary(lm(df_l2_pc$dv~df_l2_pc$ipsp_abs))

#plot e-i ratio against resting membrane potential - Panel A
corr_plot(df_l2, 'l2_sc', vm, ei, x_label='RMP', y_label='E-I Ratio', 'ei_rmp_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', vm, ei, x_label='RMP', y_label='E-I Ratio', 'ei_rmp_l2_pc.eps')

# plot non-transformed EPSP x IPSP - Panel B
corr_plot(df_l2, 'l2_sc', ipsp_abs, epsp, "EPSP (mV)", "IPSP (mV)", 'ei_corr_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', ipsp_abs, epsp, "EPSP (mV)", "IPSP (mV)", 'ei_corr_l2_pc.eps')

# plot e/ipsp x location - Panels C-D
corr_plot(df_l2, 'l2_sc', dist_psub, epsp, 'Distance (mm)', "EPSP (mV)", 'epsp_medlat_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', dist_psub, epsp, 'Distance (mm)', "EPSP (mV)", 'epsp_medlat_l2_pc.eps')
corr_plot(df_l2, 'l2_sc', dv, epsp, 'DV', "EPSP (mV)", 'epsp_dv_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', dv, epsp, 'DV', "EPSP (mV)", 'epsp_dv_l2_pc.eps')
corr_plot(df_l2, 'l2_sc', dist_psub, ipsp_abs, 'Distance (mm)', "IPSP (mV)", 'ipsp_medlat_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', dist_psub, ipsp_abs, 'Distance (mm)', "IPSP (mV)", 'ipsp_medlat_l2_pc.eps')
corr_plot(df_l2, 'l2_sc', dv, ipsp_abs, 'DV', "IPSP (mV)", 'ipsp_dv_l2_sc.eps') 
corr_plot(df_l2, 'l2_pc', dv, ipsp_abs, 'DV', "IPSP (mV)", 'ipsp_dv_l2_pc.eps')

# save results
results <- list(ei_rmp_sc, ei_rmp_pc,
                ei_sc, ei_pc,
                epsp_medlat_sc, epsp_medlat_pc, 
                ipsp_medlat_sc, ipsp_medlat_pc,
                epsp_dv_sc, epsp_dv_pc,
                ipsp_dv_sc, ipsp_dv_pc)
capture.output(results, file="supplemental_stats_pc.txt")


#-----------------------------------------------------------------------------------------------------------
#Figure 9, Supplement 2-------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#compare e/ipsp amplitude x neuron position - l1 & l2 interneuons

# subset l1 & l2 interneurons
df_l1_in <- subset(df_in, type == 'l1_in')
df_l2_in <- subset(df_in, type == 'l2_in')

#correlational stats 
epsp_medlat_l1_in <- summary(lm(df_l1_in$dist_psub~df_l1_in$epsp)) 
epsp_medlat_l2_in <- summary(lm(df_l2_in$dist_psub~df_l2_in$epsp)) 
epsp_dv_l1_in <- summary(lm(df_l1_in$dv~df_l1_in$epsp)) 
epsp_dv_l2_in <- summary(lm(df_l2_in$dv~df_l2_in$epsp))

#epsp
corr_plot(df_in, 'l1_in', dist_psub, epsp, 'Distance (mm)', y_label=Delta~EPSP, 'epsp_medlat_l1_in.eps') 
corr_plot(df_in, 'l2_in', dist_psub, epsp, 'Distance (mm)', y_label=Delta~EPSP, 'epsp_medlat_l2_in.eps')
corr_plot(df_in, 'l1_in', dv, epsp, 'DV', "EPSP (mV)", 'epsp_dv_l1_in.eps') 
corr_plot(df_in, 'l2_in', dv, epsp, 'DV', "EPSP (mV)", 'epsp_dv_l2_in.eps')

# save results
results <- list(epsp_medlat_l1_in, epsp_medlat_l2_in, 
                epsp_dv_l1_in, epsp_dv_l2_in)
capture.output(results, file="supplemental_stats_in.txt")


#-----------------------------------------------------------------------------------------------------------
#Figure 10 -------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# analyse psp amplitudes across high frequency stimulation

#create directory
setwd(wd)
create_dir(wd, 'high_frequency')

#import data
df_hf <- read.csv(file=paste(wd, "data/high_frequency.csv", sep="/"), head=TRUE, sep=",") 

#summary statistics - 10 Hz, stellates
sum_stats_hf_epsp_l2_sc <- df_hf %>%
  filter(frequency=='10hz' & type=='l2_sc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(epsp, na.rm=TRUE), 
            sd = sd(epsp, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp, na.rm=TRUE),
            iqr = IQR(epsp, na.rm=TRUE))
  
sum_stats_hf_ipsp_l2_sc <- df_hf %>%
  filter(frequency=='10hz' & type=='l2_sc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
              mean = mean(ipsp_abs, na.rm=TRUE), 
              sd = sd(ipsp_abs, na.rm=TRUE),
              se = sd/(sqrt(count)),
              median = median(ipsp_abs, na.rm=TRUE),
              iqr = IQR(ipsp_abs, na.rm=TRUE))

sum_stats_hf_ei_l2_sc <- df_hf %>%
  filter(frequency=='10hz' & type=='l2_sc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
              mean = mean(ei, na.rm=TRUE), 
              sd = sd(ei, na.rm=TRUE),
              se = sd/(sqrt(count)),
              median = median(ei, na.rm=TRUE),
              iqr = IQR(ei, na.rm=TRUE))
  
#summary statistics - 10 Hz, pyramidal cells
  sum_stats_hf_epsp_l2_pc <- df_hf %>%
  filter(frequency=='10hz' & type=='l2_pc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
              mean = mean(epsp, na.rm=TRUE), 
              sd = sd(epsp, na.rm=TRUE),
              se = sd/(sqrt(count)),
              median = median(epsp, na.rm=TRUE),
              iqr = IQR(epsp, na.rm=TRUE))
  
  sum_stats_hf_ipsp_l2_pc <- df_hf %>%
  filter(frequency=='10hz' & type=='l2_pc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
              mean = mean(ipsp_abs, na.rm=TRUE), 
              sd = sd(ipsp_abs, na.rm=TRUE),
              se = sd/(sqrt(count)),
              median = median(ipsp_abs, na.rm=TRUE),
              iqr = IQR(ipsp_abs, na.rm=TRUE))
  
  sum_stats_hf_ei_l2_pc <- df_hf %>%
  filter(frequency=='10hz' & type=='l2_pc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
              mean = mean(ei, na.rm=TRUE), 
              sd = sd(ei, na.rm=TRUE),
              se = sd/(sqrt(count)),
              median = median(ei, na.rm=TRUE),
              iqr = IQR(ei, na.rm=TRUE))  

# subset for ks test & cumulative plots
hf_10 <- df_hf[df_hf$frequency=='10hz',]
hf_20 <- df_hf[df_hf$frequency=='20hz',] 
hf_10_sc <- hf_10[hf_10$type=='l2_sc',]
hf_10_pc <- hf_10[hf_10$type=='l2_pc',] 

sc_1 <- filter(hf_10_sc, pulse=='pulse1'); sc_10 <- filter(hf_10_sc, pulse=='pulse10')
pc_1 <- filter(hf_10_pc, pulse=='pulse1'); pc_10 <- filter(hf_10_pc, pulse=='pulse10')

#compare variance between stellates/pyramidals for first/tenth pulse
ei_pulse1_ks <- ks.test(sc_1$ei, pc_1$ei)
ei_pulse10_ks <- ks.test(sc_10$ei, pc_10$ei)

#paired wilcoxon to compare epsp, ipsp & ei ratio across first & tenth pulse
wilcox_sc_epsp <- wilcox.test(sc_1$epsp, sc_10$epsp, paired = TRUE)
wilcox_sc_ipsp <- wilcox.test(sc_1$ipsp_abs, sc_10$ipsp_abs, paired = TRUE)
wilcox_sc_ei <- wilcox.test(sc_1$ei, sc_10$ei, paired = TRUE)
wilcox_pc_epsp <- wilcox.test(pc_1$epsp, pc_10$epsp, paired = TRUE)
wilcox_pc_ipsp <- wilcox.test(pc_1$ipsp_abs, pc_10$ipsp_abs, paired = TRUE)
wilcox_pc_ei <- wilcox.test(pc_1$ei, pc_10$ei, paired = TRUE)

#generate cum. freq plots for Figure 10, Panels B-C
cumulative_plot_overlay(hf_10_sc, ei, pulse, "E-I Ratio", "pulse_10_cum_freq_sc.eps")
cumulative_plot_overlay(hf_10_pc, ei, pulse, "E-I Ratio", "pulse_10_cum_freq_pc.eps")

# line plots of responses at pulse 1 v pulse 10 for Figure 10, Panel D-F
lp_hf(df_hf, ei, 'l2_sc', '10hz', y_label=expression("E-I Ratio"), "single", "ei") #Panel D
lp_hf(df_hf, ei, 'l2_pc', '10hz', y_label=expression("E-I Ratio"), "single", "ei")
lp_hf(df_hf, epsp, 'l2_sc', '10hz', y_label=expression(Delta~EPSP), "single", "epsp") #Panel E
lp_hf(df_hf, epsp, 'l2_pc', '10hz', y_label=expression(Delta~EPSP), "single", "epsp")
lp_hf(df_hf, ipsp_abs, 'l2_sc', '10hz', y_label=expression(Delta~IPSP), "single", "ipsp") #Panel F
lp_hf(df_hf, ipsp_abs, 'l2_pc', '10hz', y_label=expression(Delta~IPSP), "single", "ipsp")

# save results
results <- list(sum_stats_hf_epsp_l2_sc, sum_stats_hf_ipsp_l2_sc, sum_stats_hf_ei_l2_sc,
                sum_stats_hf_epsp_l2_pc, sum_stats_hf_ipsp_l2_pc, sum_stats_hf_ei_l2_pc,
                ei_pulse1_ks, ei_pulse10_ks,
                wilcox_sc_epsp, wilcox_sc_ipsp, wilcox_sc_ei,
                wilcox_pc_epsp, wilcox_pc_ipsp, wilcox_pc_ei)
capture.output(results, file="hf_stats.txt")


#-----------------------------------------------------------------------------------------------------------
#Figure 10, Supplement 1------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#20 Hz optogenetic stimulation of stellate/pyramidal cells

#create directory for supplemental plots/data
dir <- getwd()
create_dir(dir, 'Fig10_supplemental')
#subset data
hf_20_sc <- hf_20[hf_20$type=='l2_sc',]
hf_20_pc <- hf_20[hf_20$type=='l2_pc',] 

sc_20_1 <- filter(hf_20_sc, pulse=='pulse1'); sc_20_10 <- filter(hf_20_sc, pulse=='pulse10')
pc_20_1 <- filter(hf_20_pc, pulse=='pulse1'); pc_20_10 <- filter(hf_20_pc, pulse=='pulse10')

#summary statistics - 20 Hz, stellates
sum_stats_hf_epsp_l2_sc_20 <- df_hf %>%
  filter(frequency=='20hz' & type=='l2_sc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(epsp, na.rm=TRUE), 
            sd = sd(epsp, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp, na.rm=TRUE),
            iqr = IQR(epsp, na.rm=TRUE))

sum_stats_hf_ipsp_l2_sc_20 <- df_hf %>%
  filter(frequency=='20hz' & type=='l2_sc') %>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(ipsp_abs, na.rm=TRUE), 
            sd = sd(ipsp_abs, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ipsp_abs, na.rm=TRUE),
            iqr = IQR(ipsp_abs, na.rm=TRUE))

sum_stats_hf_ei_l2_sc_20 <- df_hf %>%
  filter(frequency=='20hz' & type=='l2_sc') %>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(ei, na.rm=TRUE), 
            sd = sd(ei, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ei, na.rm=TRUE),
            iqr = IQR(ei, na.rm=TRUE))

#summary statistics - 10 Hz, pyramidal cells
sum_stats_hf_epsp_l2_pc_20 <- df_hf %>%
  filter(frequency=='20hz' & type=='l2_pc') %>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(epsp, na.rm=TRUE), 
            sd = sd(epsp, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(epsp, na.rm=TRUE),
            iqr = IQR(epsp, na.rm=TRUE))

sum_stats_hf_ipsp_l2_pc_20 <- df_hf %>%
  filter(frequency=='20hz' & type=='l2_pc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(ipsp_abs, na.rm=TRUE), 
            sd = sd(ipsp_abs, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ipsp_abs, na.rm=TRUE),
            iqr = IQR(ipsp_abs, na.rm=TRUE))

sum_stats_hf_ei_l2_pc_20 <- df_hf %>%
  filter(frequency=='20hz' & type=='l2_pc')%>% 
  group_by(pulse) %>% 
  summarise(count = n(),
            mean = mean(ei, na.rm=TRUE), 
            sd = sd(ei, na.rm=TRUE),
            se = sd/(sqrt(count)),
            median = median(ei, na.rm=TRUE),
            iqr = IQR(ei, na.rm=TRUE))  

#compare distribution between stellates/pyramidals for first/tenth pulse
ei_pulse1_ks_20 <- ks.test(sc_20_1$ei, pc_20_1$ei)
ei_pulse10_ks_20 <- ks.test(pc_20_10$ei, sc_20_10$ei)

#20 Hz cumulative frequency plots - Panels B-C
cumulative_plot_overlay(hf_20_sc, ei, pulse, "E-I Ratio", "pulse_20_cum_freq_sc.eps")
cumulative_plot_overlay(hf_20_pc, ei, pulse, "E-I Ratio", "pulse_20_cum_freq_pc.eps")

#paired wilcoxon to compare epsp, ipsp & ei ratio across first & tenth pulse
wilcox_sc_epsp_20 <- wilcox.test(sc_20_1$epsp, sc_20_10$epsp, paired = TRUE)
wilcox_sc_ipsp_20 <- wilcox.test(sc_20_1$ipsp_abs, sc_20_10$ipsp_abs, paired = TRUE)
wilcox_sc_ei_20 <- wilcox.test(sc_20_1$ei, sc_20_10$ei, paired = TRUE)
wilcox_pc_epsp_20 <- wilcox.test(pc_20_1$epsp, pc_20_10$epsp, paired = TRUE)
wilcox_pc_ipsp_20 <- wilcox.test(pc_20_1$ipsp_abs, pc_20_10$ipsp_abs, paired = TRUE)
wilcox_pc_ei_20 <- wilcox.test(pc_20_1$ei, pc_20_10$ei, paired = TRUE)

# 20 Hz line plots - Panels D-F
lp_hf(df_hf, epsp, 'l2_sc', '20hz', y_label=expression(Delta~EPSP), "single", "epsp")
lp_hf(df_hf, epsp, 'l2_pc', '20hz', y_label=expression(Delta~EPSP), "single", "epsp")
lp_hf(df_hf, ipsp_abs, 'l2_sc', '20hz', y_label=expression(Delta~IPSP), "single", "ipsp")
lp_hf(df_hf, ipsp_abs, 'l2_pc', '20hz', y_label=expression(Delta~IPSP), "single", "ipsp")
lp_hf(df_hf, ei, 'l2_sc', '20hz', y_label=expression("E-I Ratio"), "single", "ei")
lp_hf(df_hf, ei, 'l2_pc', '20hz', y_label=expression("E-I Ratio"), "single", "ei")


#-------------------------------------------------------------------------------------------------------
#Figure 10, Supplement 2--------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
# comparison of e/i properties for single pulses versus first pulse of a stimulus train

#subset singe pulses
sc_single <- filter(hf_10_sc, pulse=="single") #10 Hz
pc_single <- filter(hf_10_pc, pulse=='single')
sc_20_single <- filter(hf_20_sc, pulse=='single') #20 Hz
pc_20_single <- filter(hf_20_pc, pulse=='single')

#paired wilcoxon to compare epsp, ipsp & ei ratio across single & first pulse - 10 Hz
wilcox_sc_ei_2 <- wilcox.test(sc_single$ei, sc_1$ei, paired = TRUE) #from peak
wilcox_pc_ei_2 <- wilcox.test(pc_single$ei, pc_1$ei, paired = TRUE)

# plots of values from single pulse v values from first pulse of train - 10 Hz, Panel A
lp_hf(df_hf, ei, 'l2_sc', '10hz',y_label=expression("E-I Ratio"), "pulse10", "og_ei")
lp_hf(df_hf, ei, 'l2_pc', '10hz', y_label=expression("E-I Ratio"), "pulse10", "og_ei")

#paired wilcoxon to compare epsp, ipsp & ei ratio across single & first pulse - 20 Hz
wilcox_sc_ei_20_2 <- wilcox.test(sc_20_single$ei, sc_20_1$ei, paired = TRUE) #from peak
wilcox_pc_ei_20_2 <- wilcox.test(pc_20_single$ei, pc_20_1$ei, paired = TRUE)

# plots of values from single pulse v values from first pulse of train - 20 Hz, Panel B
lp_hf(df_hf, ei, 'l2_sc', '20hz',y_label=expression("E-I Ratio"), "pulse10", "og_ei")
lp_hf(df_hf, ei, 'l2_pc', '20hz', y_label=expression("E-I Ratio"), "pulse10", "og_ei")

# save results
results <- list(sum_stats_hf_epsp_l2_sc_20, sum_stats_hf_ipsp_l2_sc_20, sum_stats_hf_ei_l2_sc_20,
                sum_stats_hf_epsp_l2_pc_20, sum_stats_hf_ipsp_l2_pc_20, sum_stats_hf_ei_l2_pc_20,
                ei_pulse1_ks_20, ei_pulse10_ks_20,
                wilcox_sc_epsp_20, wilcox_sc_ipsp_20, wilcox_sc_ei_20,
                wilcox_pc_epsp_20, wilcox_pc_ipsp_20, wilcox_pc_ei_20,
                wilcox_sc_ei_2, wilcox_pc_ei_2,
                wilcox_sc_ei_20_2, wilcox_pc_ei_20_2)
capture.output(results, file="hf_stats_supplemental.txt")
