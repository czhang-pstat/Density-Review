# density time series

chk_USA_df = data.frame(dens = as.vector(t(dens[[40]])), age = rep(dSup, time = dim(dens[[40]])[1]), year = rep(yrs[[40]], each = length(dSup)))

USA_dens = dens[[40]]

USA_dens = USA_dens[1:(85-9),]

var_prop = 0.9

USA_lqd = apply(USA_dens, 1, function(xx) dens2lqd(xx, dSup = dSup))
qSup = seq(0,1, length.out = length(dSup))


# expanding window

# lqd
ts_dens_lqd = matrix(NA, 90, 5)

for( i in 71:75){
ftsm_fitting = ftsm(y = fts(x = qSup, y = USA_lqd[, 1:i]), order = 10)
ftsm_ncomp = head(which(cumsum((ftsm_fitting$varprop^2)/sum((ftsm_fitting$varprop^2))) >= var_prop),1)
den_fore = forecast(object = ftsm(y = fts(x = qSup, y = USA_lqd[, 1:i]), order = ftsm_ncomp), h = 1, method = "arima")$mean$y
res_fore = lqd2dens(lqd = den_fore, lqdSup = qSup, dSup = dSup)
ts_dens_lqd[, (i-70)] = res_fore
}

# WARp
ts_dens_WARp = matrix(NA, 90, 5)
ts_qntl_WARp = matrix(NA, 90, 5)

USA_quantile = apply(USA_dens, 1, function(xx) dens2quantile(xx, dSup = dSup))
qSup = seq(0,1,length.out = length(dSup))

for( i in 71:75){
USA_acvfs = WARp_acvfs(i, i, USA_quantile, quantile.grid = qSup, 10)
den_fore_WARp = WARp_Forecast(10, USA_acvfs)

# Check monotonicity of the forecast, if so, we can use it as quantile function directly.
# print(sum(diff(den_fore_WARp)<=0))

ts_qntl_WARp[, (i-70)] = den_fore_WARp

res_fore_WARp = finite.differences(den_fore_WARp, qSup)
ts_dens_WARp[, (i-70)] = res_fore_WARp
}

for(i in 1:5){
pdf(file = paste("ts_year_", 2003+i, ".pdf", sep = ""),
    width = 6.75, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mar = c(5, 4, 2, 6)+0.1)

matplot(cbind(dSup, dSup, ts_qntl_WARp[,i]), cbind(USA_dens[(71+i),], ts_dens_lqd[,i], ts_dens_WARp[,i]), type = "l", col = c("black", "#c913ed","#009B9F"), lty = c(4,1,1), lwd = c(3,2,2), xlab = "age", ylab = "density", cex.lab = 1.2,
        cex.axis = 1.2)
grid()

par(xpd = TRUE)

labels = c("Obs.","LQD","WARp")

legend("right",inset=c(-0.24,0), legend=sapply(labels, as.expression), col=c("black", "#c913ed","#009B9F"), lty=c(4,1,1), lwd=1.5, cex = 0.9)

dev.off()
}

pdf(file = paste("ts_USA.pdf", sep = ""),
    width = 6.75, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mar = c(5, 4, 2, 5)+0.1)

do.call(matplot, list(type='n', x=dSup, y=t(USA_dens), xlab="age", ylab="density", cex.lab=1.2, cex.axis=1.2))

grid()

matlines(dSup, t(USA_dens), type = "l", col = sequential_hcl(90, "ag_Sunset"), lty = 1, lwd = 1)

par(xpd = TRUE)
colkey (side = 4, add = TRUE, col = sequential_hcl(90, "ag_Sunset"), clim = c(1933,2008), length = 0.3, clab = "year", at = c(1933, seq(1950,1990,20), 2008), adj.clab = 0)

dev.off() 
