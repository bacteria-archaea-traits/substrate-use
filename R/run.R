
#################
# Load packages #
#################

library("ggpubr")
library("png")
library("readxl")
library("tidyverse")
library("Hmisc")
library("grid")
library("agricolae") #For New Tukeys method
library("hier.part") # Hieragical partitioning

options(stringsAsFactors = FALSE)

# Load custom functions for PCAs
source("R/PCA_functions.R")
scaleFUN1 <- function(x) sprintf("%.1f", x)
scaleFUN1 <- function(x) sprintf("%.2f", x)

###########
# Process #
###########

# Prepare data frames
source("R/prep.R")

# Create fiures
source("R/figures.R")

# Create tables
source("R/tables.R")
