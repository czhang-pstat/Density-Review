% example for hades

p = path;
isExist = regexp(p, 'codes');
if isempty(isExist) == 1
  addpath(genpath('../codes/'));
end

nsub = 1000;
b = 0.01;
c = exp(0.2);
x = randraw('gompertz',[b,c],nsub)';
censor = random('binomial',1,0.8,1,nsub);

kern = 'epan';
[haz,dens]=hades(x,censor,[],kern,2,0,[],0,0,[]);


subplot(1,2,1)
plot(haz.hout,b*c.^haz.hout,'k-',haz.hout,haz.hazfun,'r--')
legend('true','estimated')
title('Hazard Function')
subplot(1,2,2)
plot(dens.dout,b*c.^dens.dout.*exp(-b*(c.^dens.dout-1)/log(c)),'k-',dens.dout,dens.dens,'r--')
title('Density Function')
legend('true','estimated')

