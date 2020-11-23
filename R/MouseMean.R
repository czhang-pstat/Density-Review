### This script contains the steps to 
### estimate the densities and means 
### for the mouse data sets.

### Import Data Sets ###
mouse1 = read.delim("./Data/Mouse/Mouse1.txt") # normal
mouse2 = read.delim("./Data/Mouse/Mouse2.txt") # Ts1Cje

mouse1_list = list()
mouse2_list = list()

for (i in 1:6) {
  mouse1_list[[i]] = mouse1[which(is.finite(mouse1[,i])), i]
  mouse2_list[[i]] = mouse2[which(is.finite(mouse2[,i])),i]
}

# Estimate and regularize densities.
mouse1_dens_list = density_estimation(mouse1_list, band_choice = "Silverman", kernel = "gaussian")

mouse2_dens_list = density_estimation(mouse2_list, band_choice = "Silverman", kernel = "gaussian")

mouse1_dens = t(mouse1_dens_list[[1]])
mouse1_dSup = mouse1_dens_list[[2]]

mouse1_dens = t(apply(mouse1_dens, 1, function(xx) RegulariseByAlpha(mouse1_dSup, xx, alpha = 0.01)))

mouse2_dens = t(mouse2_dens_list[[1]])
mouse2_dSup = mouse2_dens_list[[2]]

mouse2_dens = t(apply(mouse2_dens, 1, function(xx) RegulariseByAlpha(mouse2_dSup, xx, alpha = 0.01)))

### Evaluate Mean Densities ###

###########
### CLR ###
###########

# Normal
# choosing smoothing penalty alpha
CV1 = smoothSplinesVal(k=3, l=2, alpha=10^seq(-3,3,by=1), data=mouse1_dens, xcp=mouse1_dSup, knots=100, cores=2)

# clr transformation
dens_clr_coef1 = smoothSplines(k=3, l=2, alpha=10^seq(-3,3,by=1)[which(CV1$CVerror == min(CV1$CVerror))], data=mouse1_dens, xcp=mouse1_dSup, knots=100)

# evaluate the clr-transformed data
clr_basis1 = create.bspline.basis(range(mouse1_dSup), nbasis = 100 - 2 + 4, norder = 4)

clr1 = t(eval.fd(mouse1_dSup, fd(t(dens_clr_coef1$bspline), clr_basis1)))

clr_mean1_dens = clr2dens(mouse1_dSup, colMeans(clr1))


# Ts1Cje
# choosing smoothing penalty alpha
CV2 = smoothSplinesVal(k=3, l=2, alpha=10^seq(-3,3,by=1), data=mouse2_dens, xcp=mouse2_dSup, knots=100, cores=2)

# clr transformation
dens_clr_coef2 = smoothSplines(k=3, l=2, alpha=10^seq(-3,3,by=1)[which(CV2$CVerror == min(CV2$CVerror))], data=mouse2_dens, xcp=mouse2_dSup, knots=100)

# evaluate the clr-transformed data
clr_basis2 = create.bspline.basis(range(mouse2_dSup), nbasis = 100 - 2 + 4, norder = 4)

clr2 = t(eval.fd(mouse2_dSup, fd(t(dens_clr_coef2$bspline), clr_basis2)))

clr_mean2_dens = clr2dens(mouse2_dSup, colMeans(clr2))

##################
### Fisher-Rao ###
##################

# Normal
fr_mean1_dens = (Karcher_mean(mouse1_dens, mouse1_dSup, 1e-5,1e-2))^2

# Ts1Cje
fr_mean2_dens = (Karcher_mean(mouse2_dens, mouse2_dSup, 1e-5,1e-2))^2

###################
### Wasserstein ###
###################

qSup = seq(0,1,length.out = length(mouse1_dSup))

### Normal ###
wass_mean1 = rowMeans(apply(mouse1_dens, 1, function(xx) dens2quantile(xx, dSup = mouse1_dSup)))

wass_mean1_dens = Finite_Diff(wass_mean1, qSup)

wass_mean1_dens = wass_mean1_dens/trapzRcpp(wass_mean1, wass_mean1_dens)

### Ts1Cje ###
wass_mean2 = rowMeans(apply(mouse2_dens, 1, function(xx) dens2quantile(xx, dSup = mouse2_dSup)))

wass_mean2_dens = Finite_Diff(wass_mean2, qSup)

wass_mean2_dens = wass_mean2_dens/trapzRcpp(wass_mean2, wass_mean2_dens)

#######################
### Cross-Sectional ###
#######################

### Normal ###
cs_mean1 = colMeans(mouse1_dens)

### Ts1Cje ###
cs_mean2 = colMeans(mouse2_dens)


pdf(file = paste("EDA_mouse1_means.pdf", sep = ""), width = 6.75, height = 5)

par(mar = c(4.5, 5, 2, 9.2)+0.1)

matplot(mouse1_dSup, t(mouse1_dens), type = "l", lty = 1, lwd = 1, col = "grey", xlab = expression(paste(log,"(PM)", sep = "")), ylab = "density", cex.lab = 1.2,
        cex.axis = 1.2)
grid()

matlines(cbind(mouse1_dSup,mouse1_dSup,wass_mean1,mouse1_dSup), cbind(clr_mean1_dens, fr_mean1_dens, wass_mean1_dens, cs_mean1), type = "l", 
         col = c("#009B9F", "blue","#c913ed", "#fa9b1e"), lty = 1, lwd = 2)



par(xpd = TRUE)

labels = c("CLR","Fisher-Rao","Wasserstein","Cross-Sectional")

legend("right",inset=c(-0.46,0), legend=sapply(labels, as.expression), col=c("#009B9F", "blue","#c913ed","#fa9b1e"), lty=1, lwd=1.5, cex = 0.9)


dev.off()

pdf(file = paste("EDA_mouse2_means.pdf", sep = ""), width = 6.75, height = 5)

par(mar = c(4.5, 5, 2, 9.2)+0.1)


matplot(mouse2_dSup, t(mouse2_dens), type = "l", lty = 1, lwd = 1.5, col = "grey", xlab = expression(paste(log,"(PM)", sep = "")), ylab = "density", cex.lab = 1.2,
        cex.axis = 1.2)
grid()

matlines(cbind(mouse2_dSup,mouse2_dSup,wass_mean2,mouse2_dSup), cbind(clr_mean2_dens, fr_mean2_dens, wass_mean2_dens, cs_mean2), type = "l",
         col = c("#009B9F", "blue","#c913ed","#fa9b1e"), lty = 1, lwd = 2)


par(xpd = TRUE)

labels = c("CLR","Fisher-Rao","Wasserstein","Cross-Sectional")

legend("right",inset=c(-0.46,0), legend=sapply(labels, as.expression), col=c("#009B9F", "blue","#c913ed","#fa9b1e"), lty=1, lwd=1.5, cex = 0.9)

dev.off()

pdf(file = paste("EDA_mouse1_dens.pdf", sep = ""), width = 6.75, height = 5)

par(mar = c(5, 4, 2, 2)+0.1)

matplot(mouse1_dSup, t(mouse1_dens), type = "l", col = sequential_hcl(23, "ag_Sunset"), lty = 1, lwd = 2, xlab = expression(paste(log,"(PM)", sep = "")), ylab = "density", cex.lab = 1.2,
        cex.axis = 1.2)
grid()

dev.off()

pdf(file = paste("EDA_mouse2_dens.pdf", sep = ""), width = 6.75, height = 5)

par(mar = c(5, 4, 2, 2)+0.1)

matplot(mouse2_dSup, t(mouse2_dens), type = "l", col = sequential_hcl(23, "ag_Sunset"), lty = 1, lwd = 2, xlab = expression(paste(log,"(PM)", sep = "")), ylab = "density", cex.lab = 1.2,
        cex.axis = 1.2)
grid()

dev.off()

#############################################################
pdf(file = paste("EDA_mouse1n2_means.pdf", sep = ""), width = 6.75, height = 5)

par(mar = c(4, 5, 1, 5)+0.1, mfrow = c(2,1), oma = c(0, 0, 0, 4))

matplot(mouse1_dSup, t(mouse1_dens), type = "l", lty = 1, lwd = 1, col = "grey", xlab = NA, 
        ylab = "density", cex.lab = 1,
        cex.axis = 1.2)
grid()

matlines(cbind(mouse1_dSup,mouse1_dSup,wass_mean1,mouse1_dSup), cbind(clr_mean1_dens, fr_mean1_dens, wass_mean1_dens, cs_mean1), type = "l", 
         col = c("#009B9F", "blue","#c913ed", "#fa9b1e"), lty = 1, lwd = 2)

matplot(mouse2_dSup, t(mouse2_dens), type = "l", lty = 1, lwd = 1.5, col = "grey", xlab = expression(paste(log,"(PM)", sep = "")), ylab = "density", cex.lab = 1,
        cex.axis = 1.2)
grid()

matlines(cbind(mouse2_dSup,mouse2_dSup,wass_mean2,mouse2_dSup), cbind(clr_mean2_dens, fr_mean2_dens, wass_mean2_dens, cs_mean2), type = "l",
         col = c("#009B9F", "blue","#c913ed","#fa9b1e"), lty = 1, lwd = 2)

par(xpd = NA)

labels = c("CLR","Fisher-Rao","Wasserstein","Cross-Sectional")

legend(x=10.6, y=3, legend=sapply(labels, as.expression), col=c("#009B9F", "blue","#c913ed","#fa9b1e"), lty=1, lwd=1.5, cex = 0.9)


dev.off()
