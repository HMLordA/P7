%-----------------------------------------------
%- FORMULE de Merton.
%-----------------------------------------------

%- Formule de Merton pour le put europeen vanille (base sur BS)
function P=Merton(t,s,K,r,sigma,N)
% Here t is Tau = T-t, so t = T-tau
global kappa lambda T mu gamma
P=zeros(size(s));
for n=0:N;
    r_n = r - lambda*kappa+n*mu/t;
    sigma_n = sqrt(sigma^2+n*gamma^2/t);
    sn = s*exp(n*gamma^2/2);
    P = P + exp(-lambda*t)*(lambda*t)^n/factorial(n)*exp(r_n*t)*BS(t,sn,K,r_n,sigma_n);
end;
P = P * exp(-r*t);
end