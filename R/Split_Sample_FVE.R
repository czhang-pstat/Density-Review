#####################
### Split Samples ###
#####################

#' Function for splitting samples.
#' 
#' @param  original_dens A matrix, the observed density functions.  Rows represent observations.
#' @param  times A numeric value, number of times to split the samples.
#' 
#' @return A list that contains the row numbers of training and testing sets.
#' @export
#' 
#######################################################

Split = function(original_dens, times = 100) {
  Sample_Size = nrow(original_dens)
  
  # sample seeds
  set.seed(930)
  seeds = sample(1:10000, size = times, replace = FALSE)
  
  training = lapply(seeds, function(xx){set.seed(xx); sample(c(1:Sample_Size), size = 20, replace = FALSE)})
  testing = lapply(training, function(xx)c(1:Sample_Size)[!(c(1:Sample_Size) %in% xx)])
  
  return(list(training, testing))
}
##########################################

########################
### Split Sample FVE ###
########################

#' Function for FVE calculation.
#' 
#' @param  original_dens A matrix, the observed density functions.  Rows represent observations.
#' @param  training_trsfmd_list A list of transformed densities for training.
#' @param  testing_trsfmd_list A list of transformed densities for testing.
#' @param  testing_ind A list of row numbers of the testing sets.
#' @param  trsfm A string, methods of transformation for analysis.  Available options are "clr", "lqd", and "fr".
#' @param  times A numeric value that indicates the number of times of sample splitting.
#' @param  npc A numeric value that indicates the number of principle components used in FPCA.
#' @param  dSup A numeric vector, the grid over which `original_dens` are evaluated.
#' @param  FVE_metric A string that indicates the metric based on which the FVE is calculated.  Available options are "L2", "Wasserstein","Aitchison", and "FisherRao".
#' @param  training_psi_list A list of the training set Frechet means under the Fisher-Rao metric.  It is provided to the function to imporve efficiency.
#' 
#' @return Boxplots of FVE analysis will be written into the current working directory.
#' @export
#' 
#######################################################

SplitFVE = function(original_dens, training_trsfmd_list, testing_trsfmd_list, testing_ind,
                    trsfm, times = 100, npc, dSup, FVE_metric, training_psi_list = NA) {
  Sample_Size = nrow(original_dens)
  
  testing_dens_list = lapply(testing_ind, function(xx) original_dens[xx, ])
  
  if(trsfm == "clr" | trsfm == "fr") {
    trsfmSup = dSup
  }
  
  if(trsfm == "lqd") {
    qSup = seq(0, 1, length.out = length(dSup))
    trsfmSup = qSup
  }
  
  # FVE out of sample
  result = rep(NA, time = times)
  
  for (i in 1:times) {
    # perform FPCA
    training_sample = training_trsfmd_list[[i]]
    
    testing_sample = testing_trsfmd_list[[i]]
    
    fpcaInput = MakeFPCAInputs(IDs = rep(1:20, each=90), tVec=rep(trsfmSup,20), t(training_sample))
    fpcaOutput = FPCA(fpcaInput$Ly, fpcaInput$Lt)
    
    training_sample_mu = fpcaOutput$mu
    
    if (trsfm == "fr") {
      training_sample_mu = rep(0, length(training_sample_mu))
    }
    
    testing_sample_ctrd = t(t(testing_sample) - training_sample_mu)
    eigenfcns = as.matrix(fpcaOutput$phi[,1:npc])
    
    fpca_scores = matrix(NA, 20, npc)
    for (j in 1:npc) {
      fpca_scores[,j] = apply(testing_sample_ctrd, 1, function(xx) trapzRcpp(trsfmSup, eigenfcns[,j]*xx))
    }
    fitted = t(t(fpca_scores %*% t(eigenfcns))+ training_sample_mu)
    
    # calculate FVE in density space
    testing_sample_dens = testing_dens_list[[i]]
    
    if (trsfm == "clr") {
      fitted_dens = t(apply(fitted, 1, function(xx) clr2dens(dSup, xx)))
    }
    
    else if (trsfm == "lqd") {
      fitted_dens = t(apply(fitted, 1, function(xx) lqd2dens(xx, dSup = dSup)))
    }
    
    else if (trsfm == "fr") {
      fitted_dens = t(apply(fitted, 1, function(xx) fr_Exp_Map(dSup, xx, training_psi_list[[i]])))^2
    }
    print(paste(trsfm, mthd,i))
    result[i] = densFVE(dens = testing_sample_dens, fitted_dens = fitted_dens, metric = FVE_metric, dSup = dSup, mean_psi = training_psi_list[[i]])
  }
  return(result)
}

#######################################################
#######################################################


#####################################################
### --------------------------------------------- ###

             ##########################
             ### Start FVE Analysis ###
             ##########################

### generate indices for training and testing samples

times = 100 # number of pairs of the split samples

original_dens = yr_2008_dens

Split_ind = Split(original_dens = original_dens)
training_ind = Split_ind[[1]]
testing_ind = Split_ind[[2]]

training_sample_list = lapply(training_ind, function(xx) original_dens[xx, ])
testing_sample_list = lapply(testing_ind, function(xx) original_dens[xx, ])

training_psi_list = lapply(training_sample_list, function(xx) Karcher_mean(xx, dSup, 1e-5, 1e-2)) # this is important to improve speed

##################
### FVE of clr ###
##################
M = ncol(original_dens)
num_knots = floor(M/2)
alpha_values = 10^seq(-3,3,by=1)

# transform the data set
# choosing smoothing penalty alpha
CV = smoothSplinesVal(k=3, l=2, alpha=alpha_values, data=original_dens, xcp=dSup, knots=num_knots, cores=2)

# clr transformation
dens_clr_coef = smoothSplines(k=3, l=2, alpha=alpha_values[which(CV$CVerror == min(CV$CVerror))], data=original_dens, xcp=dSup, knots=num_knots)

# evaluate the clr-transformed data
clr_basis = create.bspline.basis(c(20.5,109.5), nbasis = num_knots - 2 + 4, norder = 4)

clr = t(eval.fd(dSup, fd(t(dens_clr_coef$bspline), clr_basis)))

training_clr_list = lapply(training_ind, function(xx) clr[xx, ])
testing_clr_list = lapply(testing_ind, function(xx) clr[xx, ])

for(trsfm in c("clr")) {
  for(mthd in c("L2", "Wasserstein","Aitchison","FisherRao")) {
    for(npc in c(1:2)) {
      varname = paste("boxdf_", trsfm, "_PC", npc, "_", mthd, sep = "")
      print(varname)
      tmp = SplitFVE(original_dens = original_dens, training_trsfmd_list = training_clr_list, testing_trsfmd_list = testing_clr_list, testing_ind = testing_ind,
                     trsfm = trsfm, times = 100, npc = npc, dSup = dSup, FVE_metric = mthd, training_psi_list = training_psi_list)
      assign(varname, value = tmp)
    }
  }
}

##################
### FVE of lqd ###
##################

training_sample_list = lapply(training_ind, function(xx) original_dens[xx, ])
testing_sample_list = lapply(testing_ind, function(xx) original_dens[xx, ])

# transform the data set

lqd = t(apply(original_dens, 1, function(xx) dens2lqd(dens = xx, dSup = dSup)))

training_lqd_list = lapply(training_ind, function(xx) lqd[xx, ])
testing_lqd_list = lapply(testing_ind, function(xx) lqd[xx, ])

for(trsfm in c("lqd")) {
  for(mthd in c("L2", "Wasserstein","Aitchison","FisherRao")) {
    for(npc in c(1:2)) {
      varname = paste("boxdf_", trsfm, "_PC", npc, "_", mthd, sep = "")
      print(varname)
      tmp = SplitFVE(original_dens = original_dens, training_trsfmd_list = training_lqd_list, testing_trsfmd_list = testing_lqd_list, testing_ind = testing_ind,
                     trsfm = trsfm, times = 100, npc = npc, dSup = dSup, FVE_metric = mthd, training_psi_list = training_psi_list)
      assign(varname, value = tmp)
    }
  }
}

#################
### FVE of fr ###
#################

training_fr_list = list()
testing_fr_list = list()

for (i in 1:100) {
  training_fr_list[[i]] = t(apply(training_sample_list[[i]], 1, function(xx) fr_Log_Map(dSup, psi1 = training_psi_list[[i]], psi2 = sqrt(xx))))
  testing_fr_list[[i]] = t(apply(testing_sample_list[[i]], 1, function(xx) fr_Log_Map(dSup, psi1 = training_psi_list[[i]], psi2 = sqrt(xx))))
}

for(trsfm in c("fr")) {
  for(mthd in c("L2", "Wasserstein","Aitchison","FisherRao")) {
    for(npc in c(1:2)) {
      varname = paste("boxdf_", trsfm, "_PC", npc, "_", mthd, sep = "")
      print(varname)
      tmp = SplitFVE(original_dens = original_dens, training_trsfmd_list = training_fr_list, testing_trsfmd_list = testing_fr_list, testing_ind = testing_ind,
                     trsfm = trsfm, times = 100, npc = npc, dSup = dSup, FVE_metric = mthd, training_psi_list = training_psi_list)
      assign(varname, value = tmp)
    }
  }
}

### Generate Boxplots ###

trsfm_label = c(rep("clr", time = 100), rep("lqd", time = 100), rep("fr", time = 100))

boxcolors = c("#B3ADB8", "#CF9FAB","#E6E4E2")
for( mthd in c("L2","Wasserstein", "FisherRao", "Aitchison")) {
  for(i in 1:2) {
    filename = paste("box",mthd,i,"pc.pdf", sep ="_")
    dfnames = paste(c("boxdf_clr_PC","boxdf_lqd_PC", "boxdf_fr_PC"), i, "_", mthd,sep = "")
    tmp = data.frame(FVE =c(eval(parse(text = dfnames[1])),eval(parse(text = dfnames[2])), eval(parse(text = dfnames[3]))), method = trsfm_label)
    
    pdf(file = filename, 6.75,5)
    #png(file = filename, 640, 480)
    par(mar = c(5,5,2,2)+0.1)
    boxplot(data=tmp, FVE~method, cex.lab = 1.2, cex.axis = 1.2, col = boxcolors, ylim = c(0.8,1))
    dev.off()
  }
}
### --------------------------------------------- ###
#####################################################

