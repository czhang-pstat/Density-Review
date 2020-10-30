addpath('toolbox')

clear all
close all
%%
% initialize
% Support of the data
Omega = 20.5:109.5;
OmegaExt = [Omega(1)-1 Omega];
% sample size
n = 40;

method = 'pchip';
n_inv = 100000;
%% 

% work under the directory GPCA.  Split, Split_Results are sub-directories under GPCA
disp('Importing data sets...')
input = 'MdVar/original_dens.csv';
dens = readmatrix(input);
dens(:,1) = [];

disp('Normalizing data sets...')
for i = 1:size(dens,1)
    dens(i,:) = dens(i,:)/sum(dens(i,:));
end
disp('Done.')
%%
[Bs, FBs] = wasserstein_barycenter_1D_smooth(dens,Omega,method,n_inv);
f = [Bs(1) Bs]; % smooth histogram of the barycenter

V_tngt = zeros(size(dens,1),length(Omega)+1);
for j = 1:size(dens,1)
    V_tngt(j,:) = logMap(dens(j,:),FBs,Omega,method); % log maps of the training data at the barycenter
end

mean_tngt = logMap(Bs,FBs,Omega,method);
writematrix(mean_tngt, 'MdVar/TangentMean.csv')
writematrix(f, 'MdVar/WassMean.csv')
writematrix(V_tngt, "MdVar/Tngt_Vec.csv")
%%
L = 2;
% Choose initialization
rng('default');
rng(1)
V0 = rand(L,length(Omega)+1); % random initialization
range_t0 = -0.1:0.01:0.1;
[PCs_gpca_iter, Scores_gpca_iter,t0_iter,residual_iter,W_residual_iter] = algo_GPCA_1D_iter(V_tngt,OmegaExt,L,V0,f,range_t0);
outpt1 = 'MdVar/PCs.csv';
writematrix(PCs_gpca_iter, outpt1);
writematrix(Scores_gpca_iter, 'MdVar/Scores.csv')

%%
sd1 = std(Scores_gpca_iter(:,1));
sd2 = std(Scores_gpca_iter(:,2));

K1 = diag([sd1, sd1*2]);
K2 = diag([sd2, sd2*2]);

GPC1 = PCs_gpca_iter(1,:);
GPC2 = PCs_gpca_iter(2,:);

gpca_plus_k_phi_PC1 = mean_tngt + K1 * [GPC1;GPC1];
gpca_minus_k_phi_PC1 = mean_tngt - K1 * [GPC1;GPC1];

gpca_plus_k_phi_PC2 = mean_tngt + K2 * [GPC2;GPC2];
gpca_minus_k_phi_PC2 = mean_tngt - K2 * [GPC2;GPC2];
%%
gpca_plus_k_phi_dens_PC1 = zeros(91, 2);
gpca_minus_k_phi_dens_PC1 = zeros(91, 2);

for i = 1:2
    gpca_plus_k_phi_dens_PC1(:,i) = pushforward_density_new(OmegaExt+gpca_plus_k_phi_PC1(i,:),f,OmegaExt, 11);
    gpca_minus_k_phi_dens_PC1(:,i) = pushforward_density_new(OmegaExt+gpca_minus_k_phi_PC1(i,:),f,OmegaExt, 11);
end

writematrix(gpca_plus_k_phi_dens_PC1, "MdVar/gpca_plus_k_phi_dens_PC1.csv")
writematrix(gpca_minus_k_phi_dens_PC1, "MdVar/gpca_minus_k_phi_dens_PC1.csv")

gpca_plus_k_phi_dens_PC2 = zeros(91, 2);
gpca_minus_k_phi_dens_PC2 = zeros(91, 2);

for i = 1:2
    gpca_plus_k_phi_dens_PC2(:,i) = pushforward_density_new(OmegaExt+gpca_plus_k_phi_PC2(i,:),f,OmegaExt, 11);
    gpca_minus_k_phi_dens_PC2(:,i) = pushforward_density_new(OmegaExt+gpca_minus_k_phi_PC2(i,:),f,OmegaExt, 11);
end

writematrix(gpca_plus_k_phi_dens_PC2, "MdVar/gpca_plus_k_phi_dens_PC2.csv")
writematrix(gpca_minus_k_phi_dens_PC2, "MdVar/gpca_minus_k_phi_dens_PC2.csv")
%%
dens_fitted_PC1 = zeros(size(dens,1),length(OmegaExt));
for i=1:size(dens,1)
    T_iter = OmegaExt+Scores_gpca_iter(i,1)*PCs_gpca_iter(1,:);
    dens_fitted_PC1(i,:) =pushforward_density_new(T_iter,f,OmegaExt,11);
end

tgnt_fitted = zeros(size(dens,1),length(OmegaExt));
for i=1:size(dens,1)
    tgnt_fitted(i,:) = logMap(dens_fitted_PC1(i,:),FBs,Omega,method);
end

writematrix(tgnt_fitted,"MdVar/gpca_fitted_tgnt_PC1.csv")
writematrix(dens_fitted_PC1,"MdVar/gpca_fitted_dens_PC1.csv")
%%
dens_fitted_PC2 = zeros(size(dens,1),length(OmegaExt));
for i=1:size(dens,1)
    T_iter = OmegaExt+Scores_gpca_iter(i,1:2)*PCs_gpca_iter(1:2,:);
    dens_fitted_PC2(i,:) =pushforward_density_new(T_iter,f,OmegaExt,11);
end


tgnt_fitted = zeros(size(dens,1),length(OmegaExt));
for i=1:size(dens,1)
    tgnt_fitted(i,:) = logMap(dens_fitted_PC2(i,:),FBs,Omega,method);
end

writematrix(tgnt_fitted,"MdVar/gpca_fitted_tgnt_PC2.csv")
writematrix(dens_fitted_PC2,"MdVar/gpca_fitted_dens_PC2.csv")

%%


