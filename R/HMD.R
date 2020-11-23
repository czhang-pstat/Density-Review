########################
### Import densities ###
########################

fls = list.files(paste(wd, "/Data/HMD/dens_est", sep = ""))[-3] # get file names and discard Belgium, 40 countries in total
ctry_abr = substr(fls, 1, 3) # get country names

dens = dens_nonGPCA = dens_GPCA = yrs = age = list()

# alpha value, by which the densities will be regularized to avoid numerical issues
a = 0.01

for (i in 1:40) {
  tmp = readMat(paste("./Data/HMD/dens_est/", fls[i], sep = "")) # Date folder is stored under Code
  yrs[[i]] = as.character(tmp$yr)
  age[[i]] = tmp$x
  dens[[i]] = t(tmp$dns)
  dens[[i]] = t(apply(dens[[i]], 1, function(xx) RegulariseByAlpha(x=age[[i]], y=xx, alpha = a))) # all densities are regularized by alpha = 0.01
}

# density support — all densities are on this common support
dSup = as.vector(age[[1]])
qSup = seq(0, 1, length.out = length(dSup))

#########################################################
### Prepare Data Sets for Modes of Variation Analysis ###
#########################################################

# use year 2008 for modes of variation analysis
yr = 2008
loc = sapply(yrs, function(xx) which(xx == as.character(yr)))

yr_2008_dens = matrix(NA, 40, 90)
for(i in 1:40){
  yr_2008_dens[i,] = dens[[i]][loc[i],]
}


pdf(file = paste("original_dens.pdf", sep = ""), width = 6.75, height = 5)
par(mar = c(5, 4, 2, 2)+0.3)
matplot(dSup, t(yr_2008_dens), type='l', lty = 1, col = sequential_hcl(23, "ag_Sunset"),
        main = NULL, xlab='age', ylab='density', cex.lab = 1.2, cex.axis = 1.2)
grid()
dev.off()

pdf(file = paste("original_dens_drg.pdf", sep = ""), width = 6.75, height = 5)
par(mar = c(5, 4, 2, 2)+0.3)
matplot(dSup, apply(yr_2008_dens, 1, function(xx) DeregulariseByAlpha(dSup, xx)), type='l', lty = 1, col = sequential_hcl(23, "ag_Sunset"),
        main = NULL, xlab='age', ylab='density', cex.lab = 1.2, cex.axis = 1.2)
grid()
dev.off()

#################################################
### Prepare Data Sets for Regression Analysis ###
#################################################

rownames(yr_2008_dens) = ctry_abr

yr_2008_dens_reg = yr_2008_dens[-which(ctry_abr == "TWN"), ]


GDP = read.csv("./Data/Econ/GDP_WB.csv", header = TRUE)
INF = read.csv("./Data/Econ/INF_WB.csv", header = TRUE)
UEM = read.csv("./Data/Econ/UEM_WB.csv", header = TRUE)

rownames(GDP) = GDP[,2]
GDP = GDP[,-c(1, 2)]
colnames(GDP) = substr(colnames(GDP), 2, 5)

rownames(INF) = INF[,2]
INF = INF[,-c(1, 2)]
colnames(INF) = substr(colnames(INF), 2, 5)

rownames(UEM) = UEM[,2]
UEM= UEM[,-c(1, 2)]
colnames(UEM) = substr(colnames(UEM), 2, 5)

GDP_aval = !is.na((GDP)[ctry_abr,])
INF_aval = !is.na((INF)[ctry_abr,])
UEM_aval = !is.na((UEM)[ctry_abr,])

# GDP
yr_1998_GDP = (GDP[ctry_abr,])[-which(ctry_abr == "TWN"), "1998"]

# inflation rate
yr_1998_INF = (INF[ctry_abr,])[-which(ctry_abr == "TWN"), "1998"]

# unemployment rate
yr_1998_UEM = (UEM[ctry_abr,])[-which(ctry_abr == "TWN"), "1998"]

EE = as.numeric(ctry_abr[-38] %in% c("BLR", "BGR", "CZE", "HUN", "POL", "RUS", "SVK", "UKR", "EST", "LTU", "LVA"))

pred_df = data.frame(gdp = yr_1998_GDP, inf = yr_1998_INF, uem = yr_1998_UEM)
pred_df_E = data.frame(gdp = yr_1998_GDP, inf = yr_1998_INF, uem = yr_1998_UEM, EE = EE)

# effect of inf
effect_inf = data.frame(gdp = rep(mean(pred_df$gdp), 9),
                     inf = quantile(pred_df$inf, prob = seq(0.1, 0.9, by=0.1)),
                     uem = rep(mean(pred_df$uem), 9))

# effect of uem
effect_uem = data.frame(gdp = rep(mean(pred_df$gdp), 9),
                     inf = rep(mean(pred_df$inf), 9),
                     uem = quantile(pred_df$uem, prob = seq(0.1, 0.9, by=0.1)))
# effect of Est
effect_EE = data.frame(gdp = rep(mean(pred_df$gdp), 2),
                     inf = rep(mean(pred_df$inf), 2),
                     uem = rep(mean(pred_df$uem), 2),
                     EE = c(1,0))

# original Eastern European countries plot
pdf(file = paste("original_dens_E.pdf", sep = ""), width = 6.75, height = 5)
#png(file = paste("original_dens_E.png", sep = ""), width = 640, height = 480)

par(mar = c(5, 5, 2, 9)+0.1)

lb = range(as.vector(yr_2008_dens_reg))[1]
ub = range(as.vector(yr_2008_dens_reg))[2]

do.call(matplot, list(type='n', x=dSup, y=seq(lb, ub, length.out = length(dSup)), main = "Observed Sample Densities", xlab="age", ylab="", cex.lab = 1.2, cex.axis = 1.2))

grid()

matlines(dSup, t(yr_2008_dens_reg), col = ifelse(ctry_abr[] %in% c("BLR", "BGR", "CZE", "HUN", "POL", "RUS", "SVK", "UKR", "EST", "LTU", "LVA"), "#009B9F", "#EE3A8C"), 
         lwd = 1.2,
         lty = 1)

par(xpd = TRUE)

labels = c("East. Europe", "Others")

legend("right",inset=c(-0.43,0), legend=sapply(labels, as.expression), col=c("#009B9F", "#EE3A8C"), lty=1, lwd=1.5, cex = 1)

dev.off()
