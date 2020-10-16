# Geodesic PCA

# export dens
write.csv(original_dens, file = "original_dens.csv")

# read itr
GPCA_itr_dens_PC1 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_PC1.csv", header = FALSE)
GPCA_itr_dens_PC2 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_PC2.csv", header = FALSE)
GPCA_itr_dens_PC3 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_PC3.csv", header = FALSE)
GPCA_itr_dens_PC4 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_PC4.csv", header = FALSE)

GPCA_itr_dens_PC1 = t(apply(GPCA_itr_dens_PC1[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))
GPCA_itr_dens_PC2 = t(apply(GPCA_itr_dens_PC2[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))
GPCA_itr_dens_PC3 = t(apply(GPCA_itr_dens_PC3[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))
GPCA_itr_dens_PC4 = t(apply(GPCA_itr_dens_PC4[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))


GPCA_itr_dens_Surf2 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_Surf2.csv", header = FALSE)
GPCA_itr_dens_Surf3 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_Surf3.csv", header = FALSE)
GPCA_itr_dens_Surf4 = read.csv(file = "../Data/GPCA_itr/GPCA_itr_dens_Surf4.csv", header = FALSE)


GPCA_itr_dens_Surf2 = t(apply(GPCA_itr_dens_Surf2[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))
GPCA_itr_dens_Surf3 = t(apply(GPCA_itr_dens_Surf3[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))
GPCA_itr_dens_Surf4 = t(apply(GPCA_itr_dens_Surf4[,-91], 1, function(xx) xx/trapzRcpp(dSup, xx)))

GPCA_itr_mean_dens = read.csv(file = "../Data/GPCA_itr/GPCA_itr_Barycenter.csv", header = FALSE)
plot(dSup, GPCA_itr_mean_dens)

GPCA_itr_PCs = read.csv(file = "../Data/GPCA_itr/GPCA_itr_PCs.csv", header = FALSE)

plot(dSup, GPCA_itr_mean_dens)

matplot(dSup, t(GPCA_itr_dens_PC1[,-91]), type = "l")
matplot(dSup, t(GPCA_itr_dens_PC2[,-91]), type = "l")
matplot(dSup, t(GPCA_itr_dens_PC3[,-91]), type = "l")
matplot(dSup, t(GPCA_itr_dens_PC4[,-91]), type = "l")

matplot(qSup, t(GPCA_itr_PCs[,-91]), type = "l")

matplot(dSup, t(GPCA_itr_dens_Surf2[,-91]), type = "l")
matplot(dSup, t(GPCA_itr_dens_Surf3[,-91]), type = "l")
matplot(dSup, t(GPCA_itr_dens_Surf4[,-91]), type = "l")
