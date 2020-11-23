#' Function for modes of variation analysis.
#' 
#' @param  origin_dens Matrix of the original densities.
#' @param  dSup Density support, a numeric vector.
#' @param  method Method: "clr", "lqd", "GPCA".
#' @param  wd A string of the path of your current working directory, only needed for GPCA.
#' @param  numCores An integer, number of cores to be used when performing clr transformation; default value is 1.
#' 
#' @return The first 2 modes of variation plot will be written to the current working directory.
#' @export
#' 
#######################################################
MdVar = function(origin_dens, dSup, method, numCores = 1, wd = NA) {
  
  k = c(1:2) # modes of variation — one and two standard deviations
  npc = 2
  
  N = nrow(origin_dens)
  M = ncol(origin_dens)
  
  ###########################
  ### Transformation FPCA ###
  ###########################
  if (method != "GPCA") {
  
  ### Start: clr ###
  if (method == "clr") {
    num_knots = floor(M/2)
    alpha_values = 10^seq(-3,3,by=1)
    
    # choosing smoothing penalty alpha
    CV = smoothSplinesVal(k=3, l=2, alpha=alpha_values, data=origin_dens, xcp=dSup, knots=num_knots, cores=numCores)
    
    # clr transformation
    dens_clr_coef = smoothSplines(k=3, l=2, alpha=alpha_values[which(CV$CVerror == min(CV$CVerror))], data=origin_dens, xcp=dSup, knots=num_knots)
    
    # evaluate the clr-transformed data
    clr_basis = create.bspline.basis(range(dSup), nbasis = num_knots - 2 + 4, norder = 4)
    
    trsfmd_dens = t(eval.fd(dSup, fd(t(dens_clr_coef$bspline), clr_basis)))
  }
  ### End: clr ###

  ### Start: Fisher-Rao ###
  else if (method == "fr") {
    mean_psi = Karcher_mean(origin_dens, dSup, 1e-5,1e-2)
    trsfmd_dens = t(apply(sqrt(origin_dens), 1, function(xx) fr_Log_Map(dSup, mean_psi, xx)))
  }
  ### End: Fisher-Rao ###
 
  ### Start: lqd ###
  else if (method == "lqd") {
    qSup = seq(0, 1, length.out = length(dSup))
    trsfmd_dens = t(apply(origin_dens, 1, function(xx) dens2lqd(dens = xx, dSup = dSup)))
  }
  ### End: lqd ###
  
  ### Start: FPCA ###
  fpcaInput = MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(dSup,N), t(trsfmd_dens))
  fpcaOutput = FPCA(fpcaInput$Ly, fpcaInput$Lt)

  ### modes of variation ###
  trsfmd_mean = fpcaOutput$mu
    
  # standard deviation
  sd = sqrt(fpcaOutput$lambda)

  for (i in 1:npc) {
    trsfmd_plus_k_phi = trsfmd_mean + t(sd[i] * diag(k) %*% matrix(rep(as.vector(fpcaOutput$phi[, i]), time = length(k)), length(k), M, byrow = TRUE))
    trsfmd_minus_k_phi = trsfmd_mean - t(sd[i] * diag(k) %*% matrix(rep(as.vector(fpcaOutput$phi[, i]), time = length(k)), length(k), M, byrow = TRUE))

    # transform back to density
    if (method == "clr") {
      plus_k_phi_dens = apply(trsfmd_plus_k_phi, 2, function(xx) clr2dens(dSup, xx))
      minus_k_phi_dens = apply(trsfmd_minus_k_phi, 2, function(xx) clr2dens(dSup, xx))
      mean_dens = clr2dens(dSup, trsfmd_mean)

      trsfmd_fitted = fitted(fpcaOutput, K=i)
      fitted_dens = apply(trsfmd_fitted, 1, function(xx) clr2dens(dSup, xx))
    }
    
    else if (method == "fr") {
      plus_k_phi_dens = apply(trsfmd_plus_k_phi, 2, function(xx) fr_Exp_Map(dSup, xx, mean_psi))
      minus_k_phi_dens = apply(trsfmd_minus_k_phi, 2, function(xx) fr_Exp_Map(dSup, xx, mean_psi))
      mean_dens = fr_Exp_Map(dSup, trsfmd_mean, mean_psi)
      
      trsfmd_fitted = fitted(fpcaOutput, K=i)
      fitted_dens = apply(trsfmd_fitted, 1, function(xx) fr_Exp_Map(dSup, xx, mean_psi))
    }
    
    else if (method == "lqd") {
      plus_k_phi_dens = apply(trsfmd_plus_k_phi, 2, function(xx) lqd2dens(xx, dSup = dSup))
      minus_k_phi_dens = apply(trsfmd_minus_k_phi, 2, function(xx) lqd2dens(xx, dSup = dSup))
      mean_dens = lqd2dens(trsfmd_mean, dSup = dSup)
      
      trsfmd_fitted = fitted(fpcaOutput, K=i)
      fitted_dens = apply(trsfmd_fitted, 1, function(xx) lqd2dens(xx, dSup = dSup))
    }

    # upper bound of the plot
    ub = max(apply(cbind(plus_k_phi_dens, minus_k_phi_dens), 2, max))
    lb = min(apply(cbind(plus_k_phi_dens, minus_k_phi_dens), 2, min))
    args = list(type='n', x=dSup, y=seq(lb, ub, length.out = M),
                 main=bquote("Modes of Variation,"~alpha~"= 0, \U00B1 1, \U00B1 2"),
                 xlab='age', ylab='density', cex.lab = 1.2, cex.axis = 1.2)

    pdf(file = paste("MdVar_", method, "_dens_ModesofVariation_PC_", i, ".pdf", sep = ""), width = 6.75, height = 5)
    
    #png(file = paste("MdVar_", method, "_dens_ModesofVariation_PC_", i, ".png", sep = ""), width = 640, height = 480)

    par(mar = c(5, 5, 2, 2)+0.3)

    do.call(plot, args)

    grid()

    polygon(x=c(dSup, rev(dSup)), y = c(minus_k_phi_dens[,2], rev(plus_k_phi_dens[,2])),
            col='lightgrey', border=NA)

    polygon(x=c(dSup, rev(dSup)), y = c(minus_k_phi_dens[,1],rev(plus_k_phi_dens[,1])),
            col='darkgrey', border=NA)

    lines(x=dSup, y=mean_dens, col='blue', lwd = 2)

    dev.off()

    ### fitted densities ###
    pdf(file = paste("MdVar_", method, "_dens_fitted_PC_", i, ".pdf", sep = ""), width = 6.75, height = 5)
    #png(file = paste("MdVar_", method, "_dens_fitted_PC_", i, ".png", sep = ""), width = 640, height = 480)
    par(mar = c(5, 5, 2, 2)+0.3)
    matplot(dSup, fitted_dens, type='l', lty = 1, col = sequential_hcl(23, "ag_Sunset"),
            main=bquote("Fitted densities: LQD,"~.(i)~"PC(s)"),
            xlab='age', ylab='density', cex.lab = 1.2, cex.axis = 1.2)
    grid()
    dev.off()
  }
  }
  ### End: FPCA ###
  
  ####################
  ### Geodesic PCA ###
  ####################
  else {
    wdScript = paste(wd, "/GPCA/", sep = "")
    wdFile = paste(wd, "/GPCA/MdVar/", sep = "")
    
    origin_dens_drglrzd = apply(origin_dens, 1, function(xx) DeregulariseByAlpha(dSup, xx))
    origin_dens_GPCA = t(apply(origin_dens_drglrzd, 2, function(xx) RegulariseByAlpha_GPCA(dSup, xx)))
    write.csv(origin_dens_GPCA, paste(wdFile, "original_dens.csv", sep = ""))
    
    run_matlab_script(paste(wdScript,"MdVar.m", sep = ""))
    
    # write testing and training samples to the current working directory, and then feed the files to Matlab for iterative GPCA
    gpca_mean_dens = as.matrix(read.csv(paste(wdFile,"WassMean.csv", sep = ""), header = FALSE))[,1:M]
    method = "gpca"
    
    for (i in 1:npc) {
      # transform back to density
      gpca_plus_k_phi_dens = as.matrix(read.csv(paste(wdFile, "gpca_plus_k_phi_dens_PC", i, ".csv", sep = ""), header = FALSE))[2:91,]
      gpca_minus_k_phi_dens = as.matrix(read.csv(paste(wdFile, "gpca_minus_k_phi_dens_PC", i, ".csv", sep = ""), header = FALSE))[2:91,]
      
      # upper bound of the plot
      ub = max(apply(cbind(gpca_plus_k_phi_dens, gpca_minus_k_phi_dens), 2, max))
      lb = min(apply(cbind(gpca_plus_k_phi_dens,gpca_minus_k_phi_dens), 2, min))
      args <- list(type='n', x=dSup, y=seq(lb, ub, length.out = M),
                   main=bquote("Modes of Variation: GPCA;"~alpha~"= 0, \U00B1 1, \U00B1 2"),
                   xlab='age', ylab='density', cex.lab = 1.2, cex.axis = 1.2)
      
      pdf(file = paste("MdVar_", method, "_dens_ModesofVariation_PC_", i, ".pdf", sep = ""), width = 6.75, height = 5)
      
      #png(file = paste("MdVar_", method, "_dens_ModesofVariation_PC_", i, ".png", sep = ""), width = 640, height = 480)
      
      par(mar = c(5, 5, 2, 2)+0.3)
      
      do.call(plot, args)
      
      grid()
      
      polygon(x=c(dSup, rev(dSup)), y = c(gpca_minus_k_phi_dens[,2], 
                                          rev(gpca_plus_k_phi_dens[,2])), #density = 30, angle = 45,
              col='lightgrey', border=NA)
      
      polygon(x=c(dSup, rev(dSup)), y = c(gpca_minus_k_phi_dens[,1], 
                                          rev(gpca_plus_k_phi_dens[,1])),
              col='darkgrey', border=NA)
      
      lines(x=dSup, y=unlist(gpca_mean_dens), col='blue', lwd = 2)
      
      dev.off()
      
      ### fitted data ###
      
      fitted_gpca = as.matrix(read.csv(paste(wdFile, "gpca_fitted_tgnt_PC", i, ".csv", sep = ""), header = FALSE))[,2:91]
      
      fitted_gpca_dens = as.matrix(read.csv(paste(wdFile, "gpca_fitted_dens_PC", i, ".csv", sep = ""), header = FALSE))[,2:91]
      
      fitted_gpca_dens = apply(fitted_gpca_dens, 1, function(xx) xx/trapzRcpp(dSup,xx))
      
      pdf(file = paste("MdVar_", method, "_dens_fitted_PC_", i, ".pdf", sep = ""), width = 6.75, height = 5)
      #png(file = paste("MdVar_", method, "_dens_fitted_PC_", i, ".png", sep = ""), width = 640, height = 480)
      
      par(mar = c(5, 5, 2, 2)+0.3)
      matplot(dSup, fitted_gpca_dens, type='l', lty = 1, col = sequential_hcl(23, "ag_Sunset"),
              main=bquote("Fitted densities,"~.(i)~"PC(s)"), 
              xlab='age', ylab='density', cex.lab = 1.2, cex.axis = 1.2)
      grid()
      dev.off()
      
    }
  }
  return(0)
}

