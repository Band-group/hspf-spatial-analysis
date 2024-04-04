#code to reproduce the work
#install packages
list.of.packages <- c("tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#start timer to compute time to run session
tic()

source("code/priors.R",verbose=FALSE)
source("code/functions.R",verbose=FALSE)
source("code/HbS_model_fit.R",verbose=FALSE)
source("code/HbS_model_diagnostics.R",verbose=FALSE)

toc()#provide time used to run the code
