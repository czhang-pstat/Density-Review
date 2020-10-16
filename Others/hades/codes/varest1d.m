function var = varest1d(n, t, q)
% Method to estimate variance of q_c(t_j) based on paper:
% Residual variance and residual pattern in nonlinear regression
% by Gasser, Sroka & Jennen-Steinmetz (1986)

a(1) = 1;
b(1) = 0;

for i = 2:(n-1)
    a(i) = (t(i+1) - t(i))/(t(i+1) - t(i-1));
    b(i) = (t(i) - t(i-1))/(t(i+1) - t(i-1));
end

a(n) = 0;
b(n) = 1;

resid(1) = 0;
for i = 2:(n-1)
    resid(i) = a(i)*q(i-1) + b(i)*q(i+1)-q(i);
end
resid(n) = 0;
var = 1./(a.^2 + b.^2 + 1).*resid.^2;