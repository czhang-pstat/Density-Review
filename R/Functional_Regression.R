#' Function for functional regression.
#' 
#' @param  predictors A matrix, the design matrix (without intercept).
#' @param  EE A logical value, set to `TRUE` to produce the effect plot of the Eastern European indicator.
#' @param  dens A matrix, the density functions.  Rows represent observations.
#' @param  trsfm A string, method of regrassion analysis.  Available options are "clr_spline", "lqd_pw", and "Frechet".
#' @param  dSup A numeric vector, the density support.
#' @param  CI A string, method of uncertainty analysis.  Available options are "pw" and "bootstrap".
#' @param  new_df A matrix, new values of the predictors of the regression model.
#' 
#' @return The regression plots will be written under the current working directory and regression summaries will be printed to the screen.
#' @export
#' 
#######################################################
Regression_Analysis = function(predictors, EE = FALSE, dens, trsfm, dSup, CI, new_df) {
  
  ### Information shared by all cases ###
  N = nrow(dens) # number of replicates
  n = ncol(dens) # number of grid points over which densities are evaluated
  
  b0 = rep(1,N) # intercept
  X = as.matrix(cbind(b0,predictors)) # design matrix, predictors contains observations by row
  
  q = ncol(X) # number of covariates, including intercept
  
  ### bootstrap initial set-up
  boot = 200 # resample 200 times
  
  # sample seeds
  set.seed(930)
  seeds = sample(1:boot, size = boot, replace = FALSE)
  
  boot_sample = array(NA, dim = c(N,n,boot), dimnames = list(paste("obs.", 1:N), paste("grid", 1:n), paste("bootstrap rep.", 1:boot)))
  boot_Beta = array(NA, dim = c(n,q,boot), dimnames = list(paste("grid", 1:n), paste("covariate", 1:q), paste("bootstrap rep.", 1:boot)))
  
  effect_color = c("#009B9F", "#00ACAF", "#35e8ee", "#4ef9ff", "#fb9bb0","#ea8098", "#e1728b", "#d65b77","#ba3654")
  
  if (trsfm == "clr_spline") {
    
    ################################
    ### Compositional Regression ###
    ################################
    
    # Following the method in Talska 2018, we only perform the regression with coefficients of smoothing splines.
    # choosing roughness penalty alpha
    cv = smoothSplinesVal(k=3, l=2, alpha=10^seq(-3,3,by=1), data=dens, xcp=dSup, knots=n/2, cores=2)
    
    nbasis = n/2-2+4 # number of internal knots: n/2-2, B-spline order: 3+1
    
    alpha_cv = 10^seq(-3,3,by=1)[which(cv$CVerror == min(cv$CVerror))]
    
    # clr transformation
    dens_coef = smoothSplines(k=3, l=2, alpha=alpha_cv, data=dens, xcp=dSup, knots=n/2)
    
    # Estimate regression coefficient matrix B
    # dim(B) = q x (g+k+1), where g is the number of internal knots, k is the spline degree, and recall that q is the number of covariates, including the intercept.
    B = solve(t(X)%*%X)%*%t(X)%*%dens_coef$bspline
    
    # Evaluate functional data objects.
    # create basis functions
    clr_basis = create.bspline.basis(range(dSup), nbasis = nbasis, norder = 4)
    
    # evaluate clr transformed data
    Y = t(eval.fd(dSup, fd(t(dens_coef$bspline), clr_basis)))
    
    Beta.fd = fd(t(B), clr_basis)
    
    # evaluate betas
    Beta = eval.fd(dSup, Beta.fd)
    
    ### Calculate point-wise 95% confidence intervals for the regression functions ###
    
    # nbasis-by-n matrix
    Theta = t(eval.basis(dSup, clr_basis))
    
    # fitted Y_hat
    Y_hat = X %*% B %*% Theta
    
    dens_fitted = t(apply(Y_hat, 1, function(xx) clr2dens(dSup, xx)))
    
    # effect
    Y_effect = as.matrix(cbind(rep(1, time = nrow(new_df)),new_df)) %*% B %*% Theta
    
    Y_effect_dens = t(apply(Y_effect, 1, function(xx) clr2dens(dSup, xx)))
    
    # residuals
    res = Y - Y_hat
    
    # estimate Cov(Y)
    cov_Y = cov(res)*(N+1)/(N-q)
    
    # roughtness penalty term, nbasis-by-n matrix
    L_Theta = t(eval.basis(evalarg = dSup, basisobj = clr_basis, Lfdobj = 2))
    
    R = matrix(NA, nbasis, nbasis)
    for(i in 1:nbasis) {
      for(j in 1:nbasis) {
        R[i,j] = trapzRcpp(dSup, L_Theta[i,]*t(L_Theta)[,j])
      }
    }
    
    J = matrix(NA, nbasis, nbasis)
    for(i in 1:nbasis) {
      for(j in 1:nbasis) {
        J[i,j] = trapzRcpp(dSup, Theta[i,]*t(Theta)[,j])
      }
    }
    
    # calculate S_phi (Y2CMap)
    S_phi = solve(Theta%*%t(Theta) + (1/alpha_cv)*R)%*%Theta # the equation on p240 in Ramsay's 2005 book has an extra \Phi in the front
    
    # calculate S_beta (C2BMap)
    S_beta = solve(kronecker(J, t(X)%*%X)+kronecker(R, (1/alpha_cv)*diag(q)))%*%kronecker(t(J), t(X))
    
    Var_vec_Beta = kronecker(t(Theta), diag(q)) %*% S_beta %*% kronecker(S_phi, diag(N))%*%kronecker(cov_Y, diag(N)) %*% t(kronecker(S_phi, diag(N))) %*% t(S_beta) %*% t(kronecker(t(Theta), diag(q)))
    
    var_Beta = matrix(NA, n, q)
    
    for (j in 1:q) {
      for (i in 1:n) {
        var_Beta[i,j] = Var_vec_Beta[q*i-(q-j),q*i-(q-j)]
      }
    }
    
    # standard deviation of Beta
    dev_Beta = sqrt(var_Beta)
    
    # bootstrap
    for(i in 1:boot) {
      set.seed(seeds[i])
      indx = sample(c(1:N), size = N, replace = TRUE)
      boot_sample[,,i] = X %*% B %*% Theta + res[indx,]
    }
    
    for(i in 1:boot) {
      boot_B = solve(t(X) %*% X) %*% t(X) %*% t(smooth.basis(dSup, t(boot_sample[,,i]),fdPar(fdobj=clr_basis, Lfdobj=2, lambda=1/alpha_cv))$fd$coefs)
      boot_Beta[,,i] = eval.fd(dSup, fd(t(boot_B), clr_basis))
    }
    
    trsfmSup = dSup
    trsfmBasis = clr_basis
    trsfmLab = "age"
    args <- list(xlab="age", ylab="", cex.lab = 1.2, cex.axis = 1.2)
    
    # R^2
    numerator_clr = sum(apply((Y - Y_hat)^2, 1, function(xxx) trapzRcpp(dSup, xxx)))
    demoninator_clr = sum(apply((t(Y) - colMeans(Y))^2, 2, function(xxx) trapzRcpp(dSup, xxx)))
    R2_clr = 1 - numerator_clr/demoninator_clr
    
    qSup = seq(0,1, length.out = length(dSup))
    quantile = t(apply(dens, 1, function(xx) dens2quantile(xx, dSup = dSup)))
    mu_quantile = colMeans(quantile)
    
    demoninator_wass = sum(apply(quantile, 1, function(xx) trapzRcpp(X = qSup, Y = (xx - mu_quantile)^2)))
    
    fitted_quantile = t(apply(dens_fitted, 1, function(xx) dens2quantile(xx, dSup = dSup)))
    numerator_wass = sum(apply((quantile -  fitted_quantile)^2, 1, function(xx) trapzRcpp(X = qSup, Y = xx)))
    R2_wass = 1 - numerator_wass/demoninator_wass
    
  }
  else if (trsfm == "lqd_pw") {
    ###########################################
    ### Regression in lqd Space, Point-Wise ###
    ###########################################
    
    # get quantile support
    qSup = seq(0, 1, length.out = length(dSup))
    # get lqd transformed data
    Y = t(apply(dens, 1, function(xx) dens2lqd(xx, dSup = dSup)))
    
    # pointwise regression
    Beta = matrix(NA, n, q)
    
    for(i in 1:n) {
      tmp_data = cbind(data.frame(pw = Y[,i]), X)
      Beta[i, ] = (lm(pw ~ .-1, data = tmp_data))$coefficients
    }
    
    ### Calculate point-wise 95% confidence intervals for the regression functions.
    # fitted Y_hat
    Y_hat = X %*% t(Beta)
    
    dens_fitted = t(apply(Y_hat, 1, function(xx) lqd2dens(xx, dSup = dSup)))
    
    # effect
    Y_effect = as.matrix(cbind(rep(1, time = nrow(new_df)),new_df)) %*% t(Beta)
    
    Y_effect_dens = t(apply(Y_effect, 1, function(xx) lqd2dens(xx, dSup = dSup)))
    
    res = Y - Y_hat
    
    # estimate Cov(Y)
    cov_Y = cov(res)*(N+1)/(N-q)
    
    # variance of regression functions
    var_Beta = matrix(NA, n, q)
    
    for (i in 1:n) {
      var_Beta[i,] = diag(solve(t(X)%*%X)*cov_Y[i,i])
    }
    
    dev_Beta = sqrt(var_Beta)
    
    # bootstrap
    for(i in 1:boot) {
      set.seed(seeds[i])
      indx = sample(c(1:N), size = N, replace = TRUE)
      boot_sample[,,i] = X %*% t(Beta) + res[indx,]
    }
    
    for(i in 1:boot) {
      for(j in 1:n) {
        tmp_data = cbind(data.frame(pw = boot_sample[,j,i]), X)
        boot_Beta[j,,i] = (lm(pw ~ .-1, data = tmp_data))$coefficients
      }
    }
    
    trsfmSup = qSup
    trsfmLab = "quantile"
    args = list(xlab="age", ylab="", cex.lab = 1.2, cex.axis = 1.2)
    
    # R^2
    numerator_clr = sum(apply((Y - Y_hat)^2, 1, function(xxx) trapzRcpp(qSup, xxx)))
    demoninator_clr = sum(apply((t(Y) - colMeans(Y))^2, 2, function(xxx) trapzRcpp(qSup, xxx)))
    R2_clr = 1 - numerator_clr/demoninator_clr
    
    qSup = seq(0,1, length.out = length(dSup))
    quantile = t(apply(dens, 1, function(xx) dens2quantile(xx, dSup = dSup)))
    mu_quantile = colMeans(quantile)
    
    demoninator_wass = sum(apply(quantile, 1, function(xx) trapzRcpp(X = qSup, Y = (xx - mu_quantile)^2)))
    
    fitted_quantile = t(apply(dens_fitted, 1, function(xx) dens2quantile(xx, dSup = dSup)))
    numerator_wass = sum(apply((quantile -  fitted_quantile)^2, 1, function(xx) trapzRcpp(X = qSup, Y = xx)))
    R2_wass = 1 - numerator_wass/demoninator_wass
  }
  
  else if (trsfm == "Frechet") {
    
    #################################
    ### Frechet Global Regression ###
    #################################
    Y_fit = wass_regress(rightside_formula = ~., Xfit_df = predictors, Ymat = dens, bandwidth = 0.01, Ytype = 'density', Sup = dSup)
    
    print.summary.WRI(summary(Y_fit))
    
    dens_fitted = confidenceBands(Y_fit, Xpred_df = predictors, type = 'density', figure = FALSE)$den_list$fpred
    
    confidence_Band = confidenceBands(Y_fit, Xpred_df = new_df, type = 'density', figure = FALSE)
    efct = confidence_Band$den_list$fpred
    LB = confidence_Band$den_list$f_lx
    UB = confidence_Band$den_list$f_ux
    
    args <- list(xlab="age", ylab="", cex.lab = 1.2, cex.axis = 1.2)
    
    if (EE) {
    # EE effect plot
    pdf(file = paste("Reg_efct_", trsfm, "E.pdf", sep = ""), width = 6.75, height = 5)
    #png(file = paste("Reg_efct_", trsfm, "E.png", sep = ""), width = 640, height = 480)
    
    par(mar = c(5, 5, 2, 10)+0.1)
    
    lb = range(as.vector(LB))[1]
    ub = range(as.vector(UB))[2]
    
    do.call(matplot, c(list(type='n'), list(x=dSup[-c(1,n)]), list(main = "Effect Plot: EE"), list(y=seq(lb, ub, length.out = length(dSup[-c(1,n)]))), args))
    
    grid()
    
    polygon(x=c(dSup[-c(1,n)], rev(dSup[-c(1,n)])), y = c(LB[1,], rev(UB[1,])), col='darkgrey', border=NA)
    
    polygon(x=c(dSup[-c(1,n)], rev(dSup[-c(1,n)])), y = c(LB[2,], rev(UB[2,])), col='lightgrey', border=NA)
    
    matlines(dSup[-c(1,n)], t(efct), lty = 1, col = c("#009B9F", "#EE3A8C"), lwd = 1.5)
    
    par(xpd = TRUE)
    
    labels = c("East. Europe", "Others")
    
    legend("right",inset=c(-0.5,0), legend=sapply(labels, as.expression), col=c("#009B9F", "#EE3A8C"), lty=1, lwd=1.5, cex = 1)
    
    dev.off()
    }
    
    else {
      # effect plot
      pdf(file = paste("Reg_efct_", trsfm, ".pdf", sep = ""), width = 6.75, height = 5)
      
      #png(file = paste("Reg_efct_", trsfm, ".png", sep = ""), width = 640, height = 480)
      
      par(mar = c(5, 5, 2, 9)+0.1)
      
      lb = range(as.vector(efct))[1]
      ub = range(as.vector(efct))[2]
      
      do.call(matplot, c(list(type='n'), list(x=dSup[-c(1,n)]), list(main = "Effect Plot"), list(y=seq(lb, ub, length.out = length(dSup[-c(1,n)]))), args))
      
      grid()
      
      matlines(dSup[-c(1,n)], t(efct), lty = 1, col = effect_color, lwd = 1.5)
      
      par(xpd = TRUE)
      
      labels = c("10%", "20%","30%","40%","50%","60%","70%","80%","90%")
      
      legend("right",inset=c(-0.4,0), legend=sapply(labels, as.expression), col=effect_color, lty=1, lwd=1.5, cex = 1.4)
      
      dev.off()
      
      # confidence band
      
      for( i in 1:9) {
        pdf(file = paste("Reg_efct", trsfm, "_CI_", i,"0th.pdf", sep = ""), width = 6.75, height = 5)
        #png(file = paste("Reg_efct", trsfm, "_CI_", i,"0th.png", sep = ""), width = 640, height = 480)
        
        par(mar = c(5, 5, 2, 2)+0.1)
        
        lb = min(LB[i,])
        ub = max(UB[i,])
        
        do.call(plot, c(list(type='n'), list(x=dSup[-c(1,n)]), list(main = bquote("Effect Plot with Confidence Band: inf at"~.(i*10)*"% quantile")),
                        list(cex.main = 1.5),
                        list(y=seq(lb, ub, length.out = length(dSup[-c(1,n)]))), args))
        
        grid()
        
        polygon(x=c(dSup[-c(1,n)], rev(dSup[-c(1,n)])), y = c(LB[i,], rev(UB[i,])), col='lightgrey', border=NA)
        
        lines(x=dSup[-c(1,n)], y=efct[i,], col=effect_color[i], lwd = 1.5)
        
        dev.off()
      }
    }
    # R^2
    new_dSup = dSup[-c(1,n)]
    new_fitted_Y = t(apply(dens_fitted, 1, function(xx) dens2clr(new_dSup,xx)))
    #new_fitted_Y_E = t(apply(dens_fitted_E, 1, function(xx) dens2clr(new_dSup,xx)))
    Y = t(apply(dens[,2:89], 1, function(xx) dens2clr(new_dSup,xx)))
    
    numerator_clr = sum(apply((Y - new_fitted_Y)^2, 1, function(xxx) trapzRcpp(new_dSup, xxx)))
    demoninator_clr = sum(apply((t(Y) - colMeans(Y))^2, 2, function(xxx) trapzRcpp(new_dSup, xxx)))
    R2_clr = 1 - numerator_clr/demoninator_clr
    
    R2_wass = summary(Y_fit)$r.square
    
    
    # fitted density plot
    pdf(file = paste("Reg_fitted_", trsfm, ".pdf", sep = ""), width = 6.75, height = 5)
    #png(file = paste("Reg_fitted_", trsfm, ".png", sep = ""), width = 640, height = 480)
    
    par(mar = c(5, 5, 2, 9)+0.1)
    
    lb = range(as.vector(dens_fitted))[1]
    ub = range(as.vector(dens_fitted))[2]
    
    do.call(matplot, c(list(type='n'), list(x=dSup[-c(1,n)]), list(main = "Fitted Densities"), list(y=seq(lb, ub, length.out = length(dSup[-c(1,n)]))), args))
    
    grid()
    
    matlines(dSup[-c(1,n)], t(dens_fitted), col = ifelse(ctry_abr[-38] %in% c("BLR", "BGR", "CZE", "HUN", "POL", "RUS", "SVK", "UKR", "EST", "LTU", "LVA"), "#009B9F", "#EE3A8C"), lwd = 1.2,
             lty = 1)
    
    par(xpd = TRUE)
    
    labels = c("East. Europe", "Others")
    
    legend("right",inset=c(-0.43,0), legend=sapply(labels, as.expression), col=c("#009B9F", "#EE3A8C"), lty=1, lwd=1.5, cex = 1)
    
    dev.off()
    
    
    trsfmSup = dSup[-c(1,n)]
    trsfmLab = "age"
  }
  
  
  ### --------------- ###
  ### Plotting Module ###
  ### --------------- ###
  
  if(trsfm != "Frechet") {
    
    ### Beta ###
    pdf(file = paste("Reg_", trsfm, "_beta0.pdf", sep = ""), width = 6.75, height = 5)
    
    par(mar = c(5, 5, 2, 2)+0.1)
    
    matplot(trsfmSup, Beta[,1], type = "l", ylab = "", xlab = trsfmLab, col=c("black"), lty = 1, lwd = 3, cex.lab = 1.2, cex.axis=1.2)
    
    abline(h = 0, col = "magenta", lty = 4, lwd = 1.5)
    
    grid()
    
    dev.off()
    
    pdf(file = paste("Reg_", trsfm, "_Beta.pdf", sep = ""), width = 6.75, height = 5)
    
    par(mar = c(5, 5, 2, 7)+0.1)
    
    matplot(trsfmSup, Beta[,-1], type = "l", ylab = "", xlab = trsfmLab, col = c(1:(q-1)), lty = c(1:(q-1)), lwd = rep(3, time = (q-1)), cex.lab = 1.2, cex.axis=1.2)
    
    abline(h = 0, col = "magenta", lty = 4, lwd = 1.5)
    
    grid()
    
    par(xpd = TRUE)
    
    labels = c()
    
    for (i in 1:(q-1)) {
      labels = c(labels, bquote(hat(beta)[.(i)]))
    }
    
    legend("right",inset=c(-0.25,0), legend=sapply(labels, as.expression), col=1:(q-1), lty = c(1:(q-1)),lwd = rep(3, time = (q-1)), cex = 1.2)
    
    dev.off()
    
    ### Uncertainty Estimation ###
    
    if(CI == "pw") {
      ### 95% Point-Wise CI
      for(i in 1:q) {
        lb = min(Beta[,i] - 2*dev_Beta[,i])
        ub = max(Beta[,i] + 2*dev_Beta[,i])
        
        par(mar = c(5, 5, 2, 2)+0.1, xpd = FALSE)
        
        pdf(file = paste("Reg_CI_",CI,"_", trsfm, "_beta", i-1, ".pdf", sep = ""), width = 6.75, height = 5)
        
        do.call(plot, c(list(type='n'), list(x=trsfmSup), list(y=seq(lb, ub, length.out = length(trsfmSup))), args))
        
        grid()
        
        polygon(x=c(trsfmSup, rev(trsfmSup)), y = c(Beta[,i] - 2*dev_Beta[,i],
                                                    rev(Beta[,i] + 2*dev_Beta[,i])), col='lightgrey', border=NA)
        
        lines(x=trsfmSup, y=Beta[,i], col='red', lwd = 1.5)
        
        abline(h = 0, col = "magenta", lty = 4, lwd = 1.5)
        
        dev.off()
      }
    }
    
    if(CI == "bootstrap") {
      ### Bootstrap
      for(j in 1:q) {
        pdf(file = paste("Reg_CI_",CI,"_", trsfm, "_beta", j-1, ".pdf", sep = ""),width = 6.75, height = 5)
        #png(file = paste("Reg_CI_",CI,"_", trsfm, "_beta", j-1, ".png", sep = ""),width = 640, height = 480)
        
        par(mar = c(5, 5, 2, 2)+0.1)
        
        lb = range(as.vector(boot_Beta[,j,]))[1]
        ub = range(as.vector(boot_Beta[,j,]))[2]
        
        do.call(matplot, c(list(type='n'), list(x=trsfmSup),list(main = bquote("Bootstrap Confidence Band:"~hat(beta)[.(j-1)])), list(y=seq(lb, ub, length.out = length(trsfmSup))), args))

        grid()
        matlines(trsfmSup, boot_Beta[,j,], col = "lightgrey")
        matlines(trsfmSup, as.matrix(Beta[,j]), col=c("red"), lty = 1, lwd = 3)
        abline(h = 0, col = "magenta", lty = 4, lwd = 1.5)
        dev.off()
      }
    }
    
    ### Effect Plot ###
    
    args <- list(xlab="age", ylab="", cex.lab = 1.2, cex.axis = 1.2)
    
    # effect plot
    pdf(file = paste("Reg_efct_", trsfm, ".pdf", sep = ""), width = 6.75, height = 5)
    #png(file = paste("Reg_efct_", trsfm, ".png", sep = ""), width = 640, height = 480)
    
    par(mar = c(5, 5, 2, 7)+0.1)
    
    lb = range(as.vector(Y_effect_dens))[1]
    ub = range(as.vector(Y_effect_dens))[2]
    
    do.call(matplot, c(list(type='n'), list(main = "Effect Plot"), list(x=dSup), list(y=seq(lb, ub, length.out = length(dSup))), args))
    
    grid()
    
    if (EE) {
      matlines(dSup, t(Y_effect_dens), lty = 1, col = c("#009B9F", "#EE3A8C"), lwd = 1.2)
      
      par(xpd = TRUE)
      
      labels = c("East. Europe", "Others")
      legend("right",inset=c(-0.3,0), legend=sapply(labels, as.expression), col=c("#009B9F", "#EE3A8C"), lty=1, lwd=1.5, cex = 0.75)
    }
    
    else {
      matlines(dSup, t(Y_effect_dens), lty = 1, col = effect_color, lwd = 1.5)
      
      par(xpd = TRUE)
      
      labels = c("10%", "20%","30%","40%","50%","60%","70%","80%","90%")
      legend("right",inset=c(-0.25,0), legend=sapply(labels, as.expression), col=effect_color, lty=1, lwd=1.5, cex = 1)
    }
    dev.off()
    
    #############################################
    # effect plot legend
    #pdf(file = paste("Reg_efct_legend.pdf", sep = ""), width = 6.75, height = 5)
    
    #par(mar = c(5, 5, 2, 3)+0.1)
    
    #plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    
    #labels = c("10%", "20%","30%","40%","50%","60%","70%","80%","90%")

    #legend("topleft", inset = c(0.3,0.18), legend=sapply(labels, as.expression), col=effect_color, lty=1, lwd=1.5, cex = 1.2, title = "Quantile")
    
    #dev.off()
    
    #############################################
    
    ### Fitted Densities ###
    pdf(file = paste("Reg_Fitted_dens_", trsfm, ".pdf", sep = ""), width = 6.75, height = 5)
    #png(file = paste("Reg_Fitted_dens_", trsfm, ".png", sep = ""), width = 640, height = 480)
    
    par(mar = c(5, 5, 2, 8.8)+0.1)
    
    lb = range(as.vector(dens_fitted))[1]
    ub = range(as.vector(dens_fitted))[2]
    
    do.call(matplot, c(list(type='n'), list(main = "Fitted Densities"), list(x=dSup), list(y=seq(lb, ub, length.out = length(dSup))), args))
    
    grid()
    
    matlines(dSup, t(dens_fitted), col = ifelse(ctry_abr[] %in% c("BLR", "BGR", "CZE", "HUN", "POL", "RUS", "SVK", "UKR", "EST", "LTU", "LVA"), "#009B9F", "#EE3A8C"), lwd = 1.2,
             lty = 1)
    
    par(xpd = TRUE)
    
    labels = c("East. Europe", "Others")
    
    legend("right", inset=c(-0.425,0), legend=sapply(labels, as.expression), col=c("#009B9F", "#EE3A8C"), lty=1, lwd=1.5, cex = 1)
    dev.off()
  }
  print(paste("clr R square: ", R2_clr, sep = ""))
  print(paste("Wasserstein R square: ", R2_wass, sep = ""))
}

