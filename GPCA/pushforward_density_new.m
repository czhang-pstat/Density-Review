function g = pushforward_density_new(T,f,OmegaExt,F)
% Function:
%      Compute the pushforward density of a histogram by a map non-
%      necessary increasing
%
% Usage:
%       h = pushforward_density(T,f,OmegaExt)
%   input:
%       T    = the map used to push the histogram
%       f     = the histogram to transport
%       OmegaExt   = the support of the histogram f to transport
%
%   output:
%       hfinal   = the histogram transported by the map T from the measure
%       f
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


epsilon = 1e-5;
int = OmegaExt(2)-OmegaExt(1);
m = 2;

% Divide the support OmegaExt into intervals for which T is monotonic
fine_OmegaExt = OmegaExt(1):1/m:OmegaExt(end);

fine_len = size(fine_OmegaExt, 2);


fine_f = interp1(OmegaExt, f, fine_OmegaExt, 'linear', 'extrap');
fine_T = interp1(OmegaExt, T, fine_OmegaExt, 'linear', 'extrap');

fine_dR = -savitzkyGolayFilt(fine_T,1,1,F)/(int/m); % int/m is the sub interval length of density support.

N = (F+1)/2; 
N_r = fine_len-N+1;
L_ind = 1:N;
R_ind = N_r:fine_len;

weights = (0:(N-1))/(N-1);

if fine_dR(N) == 1
    fine_dR(L_ind) = 1;
else
    xStart = fine_dR(N);
    dx = fine_dR(N) - 1;
    x = xStart - weights.*dx;
    fine_dR(L_ind) = flip(x);
end

if fine_dR(N_r) == 1
    fine_dR(R_ind) = 1;
else
    xStart = fine_dR(N_r);
    dx = fine_dR(N_r) - 1;
    x = xStart - weights.*dx;
    fine_dR(R_ind) = x;

end
fine_g = fine_f./abs(fine_dR);
g = interp1(fine_OmegaExt, fine_g, OmegaExt,'linear',0.01);
return
