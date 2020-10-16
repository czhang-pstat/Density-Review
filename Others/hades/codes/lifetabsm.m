% function [haz] = lifetabsm(lifetab,method,kern,hbw,hout,trunc)
% 
% Input:
%   lifetab:  lifetable, a matrix as [aux1;aux2;aux3] which are
%                 aux1:  grid on which the lifetable is calculated
%                 aux2:  number of deaths
%                 aux3:  risk set
%   method:   indicator for mortality estimate
%                 method=1 returns q_c 
%                 method=2 returns psi(q_c) (default)
%   kern:     a character string for the kernel to be used
%                 'epan' : epanechikov kernel 
%                 'gauss' : gaussian kernel  (default)
%                 'rect' : rectangular kernel
%                 'quar' : quartic kernel
%   hbw:      bandwidth for hazard function
%                 hbw > 0:  user defined bandwidth
%                 hbw = 0:  cross-validation bandwidth (default)
%   hout:     output grid for smoothed hazard function.
%             default value is a grid of length 101 as
%             0:max(x)/100:max(x).                    
%   trunc:    positive integer; truncate the smoothed hazard function
%             at the time when the number of
%             subjects at risk is smaller than the threshold trunc.
%             default is no truncation applied (trunc=0).
% 
% Output:
%   haz:      a structure for the hazard estimation, includes the 
%             following values:
%                 hazfun:    hazard function estimated at hout
%                 hout:      output grid where hazard is estimated
%                 hbw:       bandwidth used in smoothing hazard function
%                 mratio:    empirical motality ratios
%                 mr_grid:   grid of the empirical morality ratios

function [haz] = lifetabsm(lifetab,method,kern,hbw,hout,trunc)

    if isempty(kern)
      kern = 'gauss';
    end
    aux1 = lifetab(1,:);
    aux2 = lifetab(2,:);
    aux3 = lifetab(3,:);
    nbreak = length(aux1);  % input timepoints

    % Calculation of observed mortality ratio
    % aux2 is the present observed mortality ratio
    % aux3 is the risk set

    d = [aux1(2)-aux1(1),diff(aux1)];
    for i = 1:nbreak
        if ((aux3(i) > 0) && (d(i) > 0))
            if (i < nbreak)
                aux2(i) = 2*aux2(i)/(aux3(i) + aux3(i+1));
            else
                aux2(i) = 2*aux2(i)/(2*aux3(i) - aux2(i));
            end
        end
    end
    aux2 = aux2./d;

    % obtain the bandwidth by cross validation

    if hbw == 0
        hbw = cv_lwls(aux1,aux2,aux3,kern,1,1,0);
    end
    [invalid, haz1] = lwls(hbw,kern,1,1,0,aux1,aux2',aux3,hout);
    
%     haz_1 = haz1;
%     for i = 1:length(haz1)
%         if hout(i) > max(aux1)
%             haz_1(i) = NaN;
%         end
%         if (haz1(i) < 0) && (isnan(haz1(i)) == 0)
%             haz_1(i) = NaN;
%         end
%     end
    
    if method == 2    % calculate bandwidth factor for transformation psi
        
        var = varest1d(nbreak, aux1, aux2);
        bpsi = 0;
        sumvar = 0;
        for i = 1:length(haz1)
            if ((haz1(i) ~= -99) && (haz1(i) < 0))
                haz1(i) = 0;
            end
            if ((var(i) > 0) && (haz1(i) < (2/d(i))) && (haz1(i) ~= -99))
                sumvar = sumvar + var(i);
                temp = (4 - (haz1(i)*d(i))^2)^2;
                bpsi = bpsi + var(i)/temp;
            end
        end
        if sumvar > 0
            bpsi = ((16*bpsi/sumvar)^0.2);
            if (bpsi > 0)
                hbw = hbw*bpsi;
            end
        end
        [invalid, haz1] = lwls(hbw,kern,1,1,0,aux1,aux2',aux3,hout);
        hazfun = haz1*0;
        for i = 1:length(haz1)
            if (i == 1)
                d = hout(2) - hout(1);
            else
                d = hout(i) - hout(i-1);
            end
            if (haz1(i) >= 0) && (haz1(i) < (2/d))
                hazfun(i) = log((2+d*haz1(i))/(2-d*haz1(i)))/d;
            else
                hazfun(i) = NaN;
            end
        end
        
    else
        
        hazfun = haz1;
        for i = 1:length(haz1)
            if hout(i) > max(aux1)
                hazfun(i) = NaN;
            end
            if (haz1(i) < 0) && (isnan(haz1(i)) == 0)
                hazfun(i) = NaN;
            end
        end
        
    end

    if trunc > 0
        ind = find(aux3<trunc,1)-1;
        ind1 = find(hout>=aux1(ind),1);
        if ~isempty(ind1)
            hazfun(ind1:end) = NaN;
        end
    end

    haz = struct('hazfun',hazfun,'hout',hout,'hbw',hbw,'mratio',aux2,'mr_grid',aux1);
    
    
