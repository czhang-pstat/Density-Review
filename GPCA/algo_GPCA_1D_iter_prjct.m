function [fitted] = algo_GPCA_1D_iter_prjct(Omega,w,L,v_opt,f,t0)
% Function:
%      Compute principal geodesics via Iterative Geodesic approach in [Geodesic PCA versus Log-PCA of histograms in the Wasserstein
%      space. E. Cazelles, V. Seguy and al.]
%
% Usage:
%       [v_gpca_opt, t_gpca_opt,t0_opt,residual_opt,W_residual_g] = algo_GPCA_1D_iter(w,Omega,L,w0,f,range_t0)
%   input:
%       w    = log maps of the data at the barycenter f
%       Omega     = support of the data
%       L   = number of component to estimate
%       v_opt   = principle component
%       f = Wasserstein barycenter of the data
%       range_t0 = values taken by the centering variable t0
%       iter_max = Maximum number of iterations
%       iter_sub_max = Maximum number of iterations for the proximity operator
%   output:
%       v_gpca  = principal component
%       t_gpca  = position of the projection of the data onto principal
%       components
%       residual = residual error
% Authors:
%       Elsa Cazelles, Institut de Mathématiques de Bordeaux, Université
%       Bordeaux.
%       Vivien Seguy, Graduate School of Informatics, Kyoto University.
%       Jérémie Bigot, Institut de Mathématiques de Bordeaux, Université
%       Bordeaux.
%       Marco Cuturi, CREST, ENSAE, Université de Paris Saclay.
%       Nicolas Papadakis, Institut de Mathématiques de Bordeaux, CNRS.
%
% Copyright 2017 Elsa Cazelles, Vivien Seguy

[m, n]=size(w);  %w contain m data of dimension n

delta = diff(Omega);

delta = delta(1);

t_cnd = -1:0.001:1;

t_eval = zeros(m,size(t_cnd,2),L);

t_opt = zeros(m,L);

for l=1:L  %estimation of the component l=1...L
    display(['Projecting on component ' num2str(l) '/' num2str(L)]);
    v_l = v_opt(l,:);
    t0_l = t0(l);
    
    for i = 1:m
        w_i = w(i,:);

        for k = 1:size(t_cnd,2)
            t_eval(m, k,l) = sum((w_i - (t0_l + t_cnd(k))*v_l).^2 .* f .* delta);
        end
        
        t_rslt = t_eval(m, :,l);
        [~, min_loc] = min(t_rslt);
        t_opt(i,l) = t_cnd(min_loc);
    end
end

coef_opt = repmat(t0, [m 1]) + t_opt;

fitted = zeros(m,n,L);

for p=1:L
    fitted(:,:,p) = coef_opt(:,1:p) * v_opt(1:p,:);
end
