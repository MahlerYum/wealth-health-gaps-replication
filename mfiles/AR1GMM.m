function Crit = AR1GMM(alpha,eta,rho,eps,Mhat,na)

na = na;
alpha = alpha; eta=eta ; rho = rho; eps = eps;

% Approach following STY 2004 (and their NBER 2000 paper)
% 1st: Calculate Theoretical Variance Vector for all ages 
rhos = zeros(na,1);
rhos(1) = 0;
rhos(2) = 1;
for j = 2:(na-1)
    rhos(j+1) = rhos(j) + rho^(2*(j-1));
end

Vbar = zeros(na,1);
for l = 1:na
    Vbar(l) = alpha + eps + eta*(rhos(l));
end

% 2nd: Calculate Theoratical Variance of Sum of three consecutive residuals
Sbar = zeros(na-2,1);
for l = 1:(na-2)
    Sbar(l) = 9*alpha + 3*eps + ((1+rho)^2+1)*eta + eta*rhos(l);
end
    


Mbar = [Vbar ; Sbar(1)];
N = length(Mbar);
Crit = (Mhat - Mbar)' * eye(N) * (Mhat - Mbar);

end
