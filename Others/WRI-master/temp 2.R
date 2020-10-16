data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve
xpred = colMeans(predictor)

res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)
confidence_Band = confidenceBands(res, Xpred_df = predictor[200:230, ], type = 'both')

i = i +1
plot(dSup, densityCurves[i, ], type = 'l')
lines(confidence_Band$den_list$Qpred[i-199, ], confidence_Band$den_list$fpred[i-199, ], col = 'red')

confidence_Band$den_list$fpred

FMat <- t(sapply(1:31, function(i) fdapace::cumtrapzRcpp(Y = confidence_Band$den_list$fpred[i,], X = confidence_Band$den_list$Qpred[i, ])))
i = 1

i = i + 1
plot(seq(0, 1, length.out = 101), confidence_Band$quan_list$Qpred[i, ], type = 'l')
lines(FMat[i, ] + 0.05,  confidence_Band$den_list$Qpred[i, ], col = 'red')



if (any(diff_t - mean(diff_t) > 1e-7)){
        t_equal = seq(0, 1, length.out = length(t))
        QMat_raw = t(sapply(1:n, function(i) {approx(t, QMat_raw[i, ], xout = t_equal, rule = 2)$y}))
        t = t_equal
}

################################ 0215   ########################################

data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve
t = seq(0, 1, length.out = 101)
QMat <- t(sapply(1:nrow(predictor), function(i) fdadensity::dens2quantile(densityCurves[i,], dSup, t)))
xpred = colMeans(predictor)

res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)
confidence_Band = confidenceBands(res, Xpred_df = predictor[c(1,4,2,45,67,255), ], type = 'both')


res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = QMat, Ytype = 'quantile', Sup = t)
confidence_Band = confidenceBands(res, Xpred_df = predictor[c(1,4,2,45,67,255), ], type = 'both')


par(mfrow = c(1,2))
i = i + 1
plot(dSup, densityCurves[i, ], type = 'l')
plot(smoothY$Qobs[i, ], smoothY$fobs[i, ], type = 'l')


par(mfrow = c(1,2))
i = i + 1
plot(t, QMat[i, ], type = 'l')
plot(t, res$Qobs[i, ], type = 'l')

plot(res$Qobs[i, ], t, type = 'l')



data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve
xpred = colMeans(predictor)

res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, bandwidth = 0.01, Ytype = 'density', Sup = dSup)
confidence_Band = confidenceBands(res, Xpred_df = xpred, type = 'density')


data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve

res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)
globalF_table = globalFtest(res, alpha = 0.05, permutation = TRUE, numPermu = 200)

################################ 0405   ########################################
smooth_fit = ksmooth(dSup_equal, den_approx, kernel = "normal", bandwidth = diff(range(dSup))*bandwidth, x.points = dSup_equal)

x = seq(-10, 10, length.out = 100)
y = x^2+rnorm(100, sd = 10)
plot(x, y)

bandwidth = 0.1
smooth_fit = ksmooth(x, y, kernel = "normal", bandwidth = 1e-10, x.points = x)

plot(x, y)
plot(x, smooth_fit$y)
smooth_fit$y - y


data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve

t = c(seq(0, 0.5, length.out = 50), seq(0, 1, length.out = 20))

t = seq(0, 1, length.out = length(dSup))
Ysmooth = den2Q_qd(densityCurves = densityCurves, dSup = dSup, bandwidth = 0, t = t)#c(seq(0, 0.1, length.out = 10),seq(0.11, 0.89, length.out = 80),seq(0.9, 1, length.out = 10) ))

plot(smoothY$Qobs_t[i, ], smoothY$fobs_t[i, ])
lines(smoothY$Qobs[i, ], smoothY$fobs[i, ], col = 3)
lines(dSup, densityCurves[i, ], col= 2)
i = i+1


######
data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve
xpred = colMeans(predictor)

res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, bandwidth = 0.05, Ytype = 'density', Sup = dSup)
confidence_Band = confidenceBands(res, Xpred_df = xpred, type = 'both')

###

xx = matrix(rnorm(120), 20)

xx%*% solve(t(xx)%*% xx) %*% t(xx)%*% runif(20)


data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve

res <- wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)
summary(res)

write.csv(res$Ysmooth$Qobs_t, file = 'QMat.csv')
write.csv(res$Ysmooth$qobs_t, file = 'qdMat.csv')
write.csv(res$Ysmooth$qobs_prime_t, file = 'qprimeMat.csv')
write.csv(predictor, file = 'predictor.csv')
###
data(strokeCTdensity)
predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve
xpred = colMeans(predictor)

res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, bandwidth = 0.05, Ytype = 'density', Sup = dSup)
confidence_Band = confidenceBands(res, Xpred_df = predictor[1:10, ], type = 'both', delta = 0.01)

write.csv(confidence_Band$quan_list$Q_lx, 'Q_lx.csv')
write.csv(confidence_Band$quan_list$Q_ux, 'Q_ux.csv')
write.csv(confidence_Band$den_list$f_lx, 'f_lx.csv')

write.csv(confidence_Band$den_list$f_ux, 'f_ux.csv')

####

plot(seq(0, 1, length.out = 120), colMeans(equal_t$Qobs_t))
