### Setup ###

# Make sure to clone all the files under your working directory and maintain the file structure as shown on the Github repository.
# Assign your working directory path string to the variable wd.
wd = "~/Documents/FDA/MKII/Code"

setwd(wd)

# load libraries
pckg = c("fda", "fdapace", "fdadensity", "R.matlab", "robCompositions", "colorspace", "plot.matrix", "Rcpp", "matlabr")
sapply(pckg, library, character.only = TRUE)

# source scripts that contain utility functions
util = paste("./R/", c("Wasserstein_Metric.R", "Finite_Difference.R", "Density_Estimation.R", "CLR.R", "Fisher_Rao.R", "RegularizeByAlpha_GPCA.R", "densFVE.R", "MdVar.R", "Functional_Regression.R"),
             sep = "")
sapply(util, source)

WRI = list.files("./Others/WRI-master/R")
WRI_path = paste(wd, "/Others/WRI-master/R/", WRI, sep = "")
sapply(WRI_path, source)

sourceCpp(file = paste(wd, "/Others/WRI-master/src/Rcpp_D_4cell.cpp", sep = ""))
sourceCpp(file = paste(wd, "/Others/WRI-master/src/Rcpp_D_cell.cpp", sep = ""))
sourceCpp(file = paste(wd, "/Others/WRI-master/src/kernel_partialF.cpp", sep = ""))
sourceCpp(file = paste(wd, "/Others/WRI-master/src/RcppExports.cpp", sep = ""))

