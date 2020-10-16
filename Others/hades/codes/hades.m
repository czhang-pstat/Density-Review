% function [haz,dens]=hades(x,censor,Delta,kern,method,hbw,hout,trunc,dbw,dout)
% 
% Purpose: Main routine for hazard and density smoothing.
%          Nonparametric estimates of the hazard function haz and the density
%          function dens from
%          one-dimensional lifetimes x(n) and censoring indicators
%          censor(n). Hazard function haz is estimated on
%          output grid hout. Density function dens is estimated on output grid
%          dout.
% 
% Method:  Locally weighted least squares with linear local line
%          fitting applied to raw mortality ratios. Use centralized
%          mortality function estimate q_c (method=1) or the
%          transformation psi(x)=log(4/(2-x*Delta) - 1) (method=2,
%          default).
%          Linear local line fitting applied to histogram to get density
%          estimation.
% 
% References: 
%     o  From lifetables to hazard rates: The transformation approach.
%        M\"uller, H.G., Wang, J.L. and Capra, W. B. (1997) Biometrika 84, 881-892.
%     o  Analysis of oldest-old mortality: Lifetables revisited.
%        Wang, J.L., M\"uller, H.G. and Capra, W. B. (1998) Ann. Statist. 26, 126-163.
% 
% Bandwidth Choices: User defined or cross-validation. 
%                    For hazard smoothing, the bandwidth for the
%                    transformation method (method=2) is variance-stabilized
%                    method as described in the references.
%                    Case weights used are the number of subjects at risk.
% 
% ======= 
% Input:
% =======
% 
%  x(n):           lifetimes.
%                  require x(i) > 0 i=1,...,n.
%  censor(n):      integer censoring indicator.
%                      censor(i) = 0  ==> x(i) is right censored liftime
%                      censor(i) = 1  ==> x(i) is uncensored liftime
%  Delta:          binning parameter to create lifetable bins (hazard
%                  estimation) or histogram bins (density estimation).
%                  If the lifetimes are on a natural grid, for example,
%                  Delta should be the natural grid, size (one day or one year) is in days or in years. 
%                      Delta > 0:  choose this user defined Delta
%                      Delta = 0:  calculate Delta based on the input data as
%                                  Delta = max[min(|x(i)-x(j)|,for |x(i)-x(j)|>0),range(x)/1000]
%                                  (default)
%  kern:           a character string for the kernel to be used
%                      'epan' : epanechikov kernel 
%                      'gauss': gaussian kernel  (default)
%                      'rect' : rectangular kernel
%                      'quar' : quartic kernel
%  method:         indicator for mortality estimate
%                      method=1 returns q_c 
%                      method=2 returns psi(q_c) (default)
%  hbw:            bandwidth for hazard function
%                      hbw > 0:  chooses this user defined bandwidth
%                      hbw = 0:  chooses cross-validation bandwidth (default)
%  hout:           output grid for smoothed hazard function.
%                  default value is a grid of length 101 as
%                  0:max(x)/100:max(x).                    
%  trunc:          positive integer; truncate the smoothed hazard function
%                  at the time when the number of
%                  subjects at risk is smaller than the threshold trunc.
%                  default is no trunction applied (trunc=0).
%  dbw:            bandwidth for density function
%                      dbw > 0:  chooses this user defined bandwidth
%                      dbw = 0:  chooses cross-validation bandwidth (default)
%  dout:           output grid for smoothed density function.
%                  default value is a grid of length 101 as
%                  min(x):range(x)/100:max(x). 
% 
% 
% Details:   i) There are no default values for the first 2 arguments, that
%               is, they are always part of the input arguments. 
%           ii) Any non-used or optional arguments can be set to "[]" for
%               default values.
%          iii) If a lifetable is provided as input, rather than the vectors of
%               lifetimes and censoring indicators, the function lifetabsm.m can
%               be used to estimate the hazard function.
% 
% ========
% Output:
% ========
% 
%  haz:            a structure for the hazard estimation, includes the 
%                  following values:
%                      hazfun:    hazard function estimated at hout
%                      hout:      output grid where hazard is estimated
%                      hbw:       bandwidth used in smoothing hazard function
%                      mratio:    empirical motality ratios
%                      mr_grid:   grid of the empirical morality ratios
%  dens:           a structure for the density estimation, includes the
%                  following values:
%                      dens:      density function estimated at dout
%                      dout:      output grid where density is estimated
%                      dbw:       bandwidth used in smoothing density function
% 
% Includes the following steps:
%   1)  Calculation of lifetable
%   2)  Calculation of empirical mortality ratios
%   3)  Smoothing step for hazard function 
%   4)  Calculate transformation psi
%   5)  Calculation of histogram
%   6)  Smoothing step for density function

function [haz,dens]=hades(x,censor,Delta,kern,method,hbw,hout,trunc,dbw,dout)

    % sorting the input lifetimes
    [x,ind] = sort(x);
    censor = censor(ind);
    
    % set defaults
    if isempty(Delta) || Delta == 0
        diffx = diff(x);
        Delta = max(min(diffx(diffx>0)),range(x)/1000);
    end
    if isempty(kern)
        kern = 'gauss';
    end
    if isempty(method)
        method = 2;
    end
    if isempty(hbw)
        hbw = 0;
    end
    if isempty(hout)
        hout = 0:max(x)/100:max(x);
    end
    if isempty(trunc)
        trunc = 0;
    end
    if isempty(dbw)
        dbw = 0;
    end
    if isempty(dout)
        dout = min(x):range(x)/100:max(x);
    end
    
    %  Calculation of Life Table
    %  -----------------------
    %  aux1 is lifetable grid
    %  aux2 is number of deaths
    %  aux3 is risk set

    n = length(x);
    aux1 = 0:Delta:max(x);
    if max(aux1) < max(x)
        aux1 = [aux1,max(x)];
    end
    m = length(aux1);
    aux2 = zeros(1,m);
    aux3 = zeros(1,m);
    aux3(1) = n;
    for i = 2:m
        ind1 = find(x>aux1(i-1));
        ind2 = find(x<=aux1(i));
        aux2(i) = sum(censor(intersect(ind1,ind2)));
        aux3(i) = n-ind1(1)+1;
    end
    lifetab = [aux1; aux2; aux3];

    [haz] = lifetabsm(lifetab,method,kern,hbw,hout,trunc);
    
    [dens] = hist2fx(x, Delta, dbw, kern, dout);

    
