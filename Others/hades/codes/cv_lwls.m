% function [bopt] = cv_lwls(xin,yin,win,kernel,nwe,npoly,nder)
% Purpose: Cross-validation routine used to select smoothing
%          parameter bopt of locally weighted least squares subroutine.
%          Good only for nder=0 (no derivative) and npoly=1 (local linear fitting).
%          Subroutine for dens1d and regr1d.
% Method:  Weighted least squares fitted locally in windows
%          [xou-bw, xou+bw]. Bandwidth bopt estimated by
%          cross validation.
% Details: Locally weighted least squares routine estimates
%          a function you(m)=f[xou(m)] from one-dimensional data
%          pairs [xin(n),yin(n)], with case weights win(n).
%          Weight Function: See below.
% 
% input xin(1,n)    vector of of x-coordinate (predictor)
% input yin(1,n)    vector of y-coordinates of input data
% input win(1,n)    vector of case weights for input data
% input kern:       a character string for the kernel to be used
%                       'epan' : epanechikov kernel
%                       'gauss' : gaussian kernel
%                       'rect' : rectangular kernel
%                       'quar' : quartic kernel
% input nwe:        degree of opt polynomial
% input npoly       order of polynomial to be fitted locally
%                   Requirement: npoly>nder and
%                                npoly-nder=odd
%                                npoly <= 5
% input nder        order of derivative to be estimated
% 
% output bopt       optimum cross-validation bandwidth
%
% Includes the following steps:
% 1) Basic setup
% 2) Cross-Validation step 
% 3) Finding the final bandwidth(if preliminary chosen bw isn't first or last one)

function bopt = cv_lwls(xin,yin,win,kernel,nwe,npoly,nder)

    n = length(xin);
    r = range(xin);
    dstar = minb(xin,3);
    if dstar > r/4
        dstar = dstar*0.75;
        display(['Warning: the min bandwidth choice is too big, reduce to ' num2str(h0) '!']);
    end
    
    % creat 10 bandwidth candidates
    nbw=11;
    bw=zeros(1,nbw-1);
    for i=1:nbw-1
        bw(i)=2.5*r/n*(n/2.5)^((i-1)/(nbw-1));
    end   
    bw = bw-min(bw)+dstar;
    

    % Cross-Validation step 
    ilow=1;
    low=10000000;
    w = win;
    cv = [];
    for ibw = 1:nbw-1
        cv(ibw)=0;
        count=0;
        f=0;
        icv=1;
        while (icv <= n) && (f==0)
            if win(icv) ~= 0
               win=w;
               win(icv)=0;
               xcv=xin(icv);
               [temp,ypred]=lwls(bw(ibw),kernel,nwe,npoly,nder,xin,yin',win,xcv);
               if (temp ~= 1)
                  count=count+1;
                  ypred=w(icv)*(yin(icv)-ypred)*(yin(icv)-ypred);
                  cv(ibw)=cv(ibw) + ypred;
               else
                  f=1;
               end
            end 
            icv=icv+1;
        end 
        if (count > 0) && (f == 0) 
           cv(ibw)=cv(ibw)/count;
        else 
           cv(ibw)=10000000000;
        end   
        if cv(ibw) < low  
           low=cv(ibw);
           ilow=ibw;
        end
    end
    
    % Finding the final bandwidth
    ymat = [];
    if (ilow ~=1) && (ilow ~=nbw-1)
        for i=1:3
            j = ilow - 2 + i;
            xmat(i,1)=1;
            xmat(i,2)=bw(j);
            xmat(i,3)=bw(j)*bw(j);
            ymat(i)=cv(j);
        end
        coeff=pinv(xmat'*xmat)*xmat'*ymat';
        bopt=-0.5*coeff(2)/coeff(3);

        if bopt < bw(ilow-1)
            bopt=bw(ilow-1);
        end
        if bopt > bw(ilow+1)
            bopt=bw(ilow+1);
        end
    else
        bopt=bw(ilow);
    end

    fprintf(1,['CV choice for bandwidth: ' num2str(bopt) '\n']);

