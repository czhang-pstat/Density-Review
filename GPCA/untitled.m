train1 = readmatrix("MdVar/original_dens.csv");
train1(:,1) = [];
train2 = readmatrix("MdVar/original_dens(old).csv");
train2(:,1) = [];

disp('Normalizing data sets...')
for i = 1:40
    train1(i,:) = train1(i,:)/sum(train1(i,:));
    train2(i,:) = train2(i,:)/sum(train2(i,:));
end
disp('Done.')

Sinv = 0:1/n_inv:1;

q1 = histogram_pseudo_inverse_smooth(train1(1,:),Omega,Sinv, "linear");
q2 = histogram_pseudo_inverse_smooth(train2(1,:),Omega,Sinv, "linear");


figure;

subplot(1,2,1);
plot(Sinv, q1)
title("Quantile Function of AUS (Regularized)");

subplot(1,2,2);
plot(Sinv, q2)
title("Quantile Function of AUS (Deregularized)");

%%
method = 'pchip';
n_inv = 100000; 
[Bs1, FBs1] = wasserstein_barycenter_1D_smooth(train1,Omega,method,n_inv);
[Bs2, FBs2] = wasserstein_barycenter_1D_smooth(train2,Omega,method,n_inv);
f1 = [Bs1(1) Bs1];
f2 = [Bs2(1) Bs2];

writematrix(Bs1, 'MdVar/Bs1.csv')

subplot(1,2,1)
plot(OmegaExt, f1)
title("Barycenter (Regularized)");

subplot(1,2,2)
plot(OmegaExt, f2)
title("Barycenter (Deregularized)");
%% Compute log maps

display('Log-maps of the data at the barycenter')

V1 = zeros(size(train1,1),length(Omega)+1);
for i = 1:size(train1,1)
    V1(i,:) = logMap(train1(i,:),FBs1,Omega,method); % log maps of the data at the barycenter
end

V2 = zeros(size(train2,1),length(Omega)+1);
for i = 1:size(train2,1)
    V2(i,:) = logMap(train2(i,:),FBs2,Omega,method); % log maps of the data at the barycenter
end

figure;

subplot(1,2,1)
plot(OmegaExt, V1)
title("Tangent Vector (Regularized)");
subplot(1,2,2)
plot(OmegaExt, V2)
title("Tangent Vector (Deregularized)");

%%
L = 1;
% Choose initialization
rng('default');
rng(1)
V0 = rand(L,length(Omega)+1); % random initialization
range_t0 = -0.1:0.01:0.1;
[v_gpca_iter1, t_gpca_iter1,t0_iter1,residual_iter1,W_residual_iter1] = algo_GPCA_1D_iter(V1,OmegaExt,L,V0,f1,range_t0);
[v_gpca_iter2, t_gpca_iter2,t0_iter2,residual_iter2,W_residual_iter2] = algo_GPCA_1D_iter(V2,OmegaExt,L,V0,f2,range_t0);

display("done")

%%
T1 = zeros(size(train1,1),length(Omega)+1);
h_iter1 = zeros(size(train1,1),length(OmegaExt));
for i=1:size(train1,1)
    T_iter1 = OmegaExt+t_gpca_iter1(i,1)*v_gpca_iter1(1,:);
    T1(i,:) = T_iter1 + OmegaExt;
    h_iter1(i,:) =pushforward_density_new(T_iter1,f1,OmegaExt, 11);
end

T2 = zeros(size(train2,1),length(Omega)+1);
h_iter2 = zeros(size(train2,1),length(OmegaExt));
for i=1:size(train2,1)
    T_iter2 = OmegaExt+t_gpca_iter2(i,1)*v_gpca_iter2(1,:);
    T2(i,:) = T_iter2 + OmegaExt;
    h_iter2(i,:) =pushforward_density_new(T_iter2,f2,OmegaExt, 11);
end

figure;

subplot(1,2,1)
plot(OmegaExt, T1-OmegaExt)

subplot(1,2,2)
plot(OmegaExt, T2-OmegaExt)
writematrix(h_iter1, 'MdVar/h1.csv')
figure;

subplot(1,2,1)
plot(OmegaExt, h_iter1)
title("Fitted Density (Regularized)");
subplot(1,2,2)
plot(OmegaExt, h_iter2)
title("Fitted Density (Deregularized)");