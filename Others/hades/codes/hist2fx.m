% function [dens] = hist2fx(x, Delta, dbw, kern, dout)
% 
% hist2fx is a function that first generates the histogram of the data and then
% smooths the histogram to get a density estimate fx
% 
% Input:
%   x:          1 * n vector of which the density function is estimated.
%   Delta:      histogram bins (density estimation).
%               If the lifetimes are on a  natural grid, for example daily or annual,
%               Delta should be the natural grid span. 
%                   Delta > 0:  choose this user defined Delta
%                   Delta = 0:  calculate Delta based on the input data as
%                               Delta = max[min(|x(i)-x(j)|,for |x(i)-x(j)|>0),range(x)/1000]
%                               (default)
%   dbw:        bandwidth for density function
%                   dbw > 0:  chooses this user defined bandwidth
%                   dbw = 0:  chooses cross-validation bandwidth (default)
%   kern:       a character string for the kernel to be used
%                   'epan' : epanechikov kernel 
%                   'gauss' : gaussian kernel   (default)
%                   'rect' : rectangular kernel
%                   'quar' : quartic kernel
%   dout:       output grid for smoothed density function.
%               default value is a grid of length 101 as
%               min(x):range(x)/100:max(x). 
% 
% Output:
%   dens:       a structure for the density estimation, includes the
%               following values:
%                   dens:      density function estimated at dout
%                   dout:      output grid where density is estimated
%                   dbw:       bandwidth used in smoothing density function

function [dens] = hist2fx(x, Delta, dbw, kern, dout)

    if isempty(kern)
       kern = 'gauss';
    end

    nbins = range(x)/Delta;
    if nbins ~= ceil(nbins)
        Delta = range(x)/ceil(nbins);
    end
    
    histgrid = min(x)-Delta/2:Delta:max(x)+Delta/2;
    xhist_n = histc(x,histgrid);
    yin = xhist_n(1:end-1)/length(x)/Delta;
    xin = min(x):Delta:max(x);
        
    if dbw == 0
        dbw = cv_lwls(xin,yin,ones(1,length(xin)),kern,1,1,0);
    end
    
    [invalid, fxtmp] = lwls(dbw,kern,1,1,0,xin,yin',ones(1,length(xin)),dout);
    fx = fxtmp/trapz(dout,fxtmp);
    
    dens = struct('dens',fx,'dout',dout,'dbw',dbw);
    
    
