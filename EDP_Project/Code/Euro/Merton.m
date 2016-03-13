%-----------------------------------------------
%- FORMULE de Merton.
%-----------------------------------------------

%- Formule de Black et scholes pour le put europeen vanille.
function P=Merton(t,s,K,r,sigma,N)
% Here t is Tau = T-t, so t = T-tau
global kappa lambda T mu gamma
P=zeros(size(s));
for n=0:N;
    r_n = r - lambda*kappa+n*mu/t;
    sigma_n = sqrt(sigma^2+n*gamma^2/t);
    P = P + exp(-lambda*t)*(lambda*t)^n/factorial(n)*exp(r_n*t)*BS(t,s,K,r_n,sigma_n);
end;
P = P * exp(-r*t);
end