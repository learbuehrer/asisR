###############################################################################
## load all packages needed:
###############################################################################
library(abn)
library(biostatUZH)
library(bnlearn)
library(brant)
library(broom.mixed)
library(car)
library(caret)
library(cobalt)
library(coda)
library(DHARMa)
library(forcats)
library(ggplot2)
library(ggExtra)
library(ggalluvial)
library(ggpubr)
library(glmmTMB)
library(nnet)
library(splines)
library(scales)
library(stringr)
library(sjPlot)
library(tableone)
library(xtable)
library(reshape2)
library(missForest)
library(mice)

# packages for mcmc bootstrap
library(foreach)
library(doParallel)
library(mcmcabn)

# load as last package to avoid masking with other packages:
library(dplyr)

#############################################################################
# define color palettes
#############################################################################

colorPalette <- c("#0064a6","#583119", 
                  "#f0b605", "#006633", "#83b819",
                  "#9a9a9c", "#6a205f", "#4b5d68",
                  "#d54e12", "darkred","black")

bluePalette <- c("#0028A5",  "#3353B7","#667EC9","#99a9db",
                 "#ccd4ed", "lightblue","#cfe8ec", "#a3adb7","#c8ced4", 
                 "#9ed0d9", "#6bb7c7" ,"#3c9fb6", "#0b82a0",
                 "#004573","#0064a6", "#0083d9","grey90")


colorblindPalette <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky Blue
  "#009E73",  # Bluish Green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",   # Reddish Purple
  "grey50"
)
darkColorPalette <- c(
  "#B07C00",  # dunkleres Orange
  "#3A8BCD",  # dunkleres Himmelblau
  "#007F5F",  # dunkleres Grün
  "#BBAF00",  # dunkleres Gelb
  "#005A8C",  # dunkleres Blau
  "#A04000",  # dunkleres Rot
  "#994C80",   # dunkleres Violett
  "grey50"
)
################################################################################
# source all needed scripts:
################################################################################
# read in the data:
dat_s <- readRDS("data/synthetic_dataset.rds")

source("scripts/00_functions.R")
source("scripts/00_mcmc_functions.R")
source("scripts/00_imputation.R")
