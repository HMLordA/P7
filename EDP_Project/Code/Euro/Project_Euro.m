%- DF - Black et Scholes 
%- European options
%- Olivier Bokanowski
clear

%------------------------
%- DONNEES FINANCIERES / FINANCIAL DATA
%------------------------
global  K r sigma T Smin Smax lambda mu gamma kappa
K=100; sigma=0.2; r=0.1; T=1;  Smin=20; Smax=200; lambda = 0.0; mu = 0.0; gamma = 1.0; 

kappa = exp(mu+gamma^2/2)-1; % JCD : expectancy of eta, which is log-normal

%- IC : Initial Condition   = function  u0
%- BD : Boundary Conditions = functions ul, ur
%- ==> COMPLETE definition of functions u0, ul, ur (inline definitions: see below)
global ul ur g
u0= @(s) max(K-s,0);		%- Initial values (payoff function)
ul= @(t) K*exp(-r*t)-Smin;	%- ul= left  value, at Smin
ur= @(t) 0;			%- ur= right value, at Smax
g= @(eta) exp(-(log(eta)-mu)^2/(2*gamma^2))/(sqrt(2*pi)*gamma*eta) ;


%------------------------
%- DONNEES NUMERIQUES / NUMERICAL DATA
%------------------------
global I N p nMerton
I=20; N=40; p = 10; nMerton = 10;
%I=2*10; N=I*I/10; 

SCHEMA='EE'; 		%- 'EE' or 'EI' or 'CN' 
CENTRAGE='CENTRE'; 	%- 'CENTRE', 'DROIT', 'GAUCHE' 

%- Parameters for the graphics:
global Xmin Xmax Ymin Ymax
Xmin=log(Smin/K); Xmax=log(Smax/K); Ymin=-20; Ymax=K;
err_scale=0; %- Echelle pour le graphe d'erreur.
deltan=N/10; %- Eventuellement, Affichage uniquement tous les deltan pas.

%- Printing some data:
fprintf('sigma=%5.2f, r=%5.2f, Smax=%5.2f\n',sigma,r,Smax);
fprintf('Maillage I= %5i, N=%5i\n',I,N);
fprintf('CENTRAGE : %s\n',CENTRAGE)
fprintf('SCHEMA: %s\n',SCHEMA)

%-----------------------------------------------
%- Black - Scholes formula
%-----------------------------------------------
%- Formule de Black et scholes pour le put europeen vanille : fichier BS.m
%- function y=BS(t,s,K,r,sigma)


%--------------------
%- MAILLAGE / MESH
%--------------------
%- FILL : dt, h, s (time step, mesh step, mesh) 
global h
dt=T/N; 		%- pas de temps  / time step
h=(Xmax-Xmin)/(I+1); 	%- pas d'espace  / mesh step
%s=Smin+(1:I)'*h; 	%- maillage      / mesh : column vector of size I, containing the mesh values s_i = Smin + i*h
s=K*exp(Xmin+(1:I)'*h);

global x
x= @(i) Xmin + i*h;

%- CFL COEFFICIENT 
%COMPLETE
cfl=dt/h^2 * (sigma*Smax)^2;
fprintf('CFL : %5.3f\n',cfl); 


%--------------------------
%- Initializations:  matrix A, function q(t)
%--------------------------

switch CENTRAGE

case 'CENTRE';  %- CENTERED APPROXIMATION

  %- FILL IN / COMPLETER matrice A, vecteurs alpha et bet de taille I, et fonction q
  A=zeros(I,I);
  G=zeros(I,I);
  %FILL the values of A(i,i), A(i,i-1), A(i,i+1)
  alpha=sigma^2/2/h^2;
  bet=(r-lambda*kappa+sigma^2/2)/h;
  for i=1:I;   A(i,i) = 2*alpha + r + lambda; end;
  for i=2:I;   A(i,i-1) = -alpha + bet/2; end;
  for i=1:I-1; A(i,i+1) = -alpha - bet/2; end;
  
  for i=1:I-1;
      %G(i,i)=0;
      for j=1:I-1;
          G(i,j) = g(exp(x(j-i)))*exp(x(j-i));
      end;
  end;
          

  % FILL IN
  q = @(t) [(-alpha + bet/2)*ul(t);  zeros(I-2,1);  (-alpha - bet/2)* ur(t)];

  
% JCD : not implemented
% case 'DROIT';	%- FORWARD DIFFERENCES
% 
%   A=zeros(I,I);
%   alpha=sigma^2/2 * s.^2 /h^2;
%   bet=r*s/h;
%   for i=1:I;   A(i,i) = 2*alpha(i) + bet(i) + r; end;
%   for i=2:I;   A(i,i-1) = -alpha(i) ; end;
%   for i=1:I-1; A(i,i+1) = -alpha(i) - bet(i); end;
% 
%   q = @(t) [(-alpha(1))* ul(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* ur(t)];
% 
% case 'GAUCHE';	%- BACKWARD DIFFERENCES
% 
%   A=zeros(I,I);
%   alpha=sigma^2/2 * s.^2 /h^2;
%   bet=r*s/h;
%   for i=1:I;   A(i,i) = 2*alpha(i) - bet(i) + r; end;
%   for i=2:I;   A(i,i-1) = -alpha(i) + bet(i) ; end;
%   for i=1:I-1; A(i,i+1) = -alpha(i) ; end;
% 
%   q = @(t) [(-alpha(1) +bet(1))* ul(t);  zeros(I-2,1);  (-alpha(end))* ur(t)];

otherwise 

  fprintf('this CENTRAGE not programmed !'); abort

end
 
%--------------------
%- Initialiser P et graphique
%--------------------
P=u0(s);
ploot(0,s,P);
fprintf('waiting for ''Enter'''); input('');


%--------------------
%- BOUCLE PRINCIPALE / MAIN LOOP
%--------------------
%- starting cputime counter
tic(); 

Id=eye(size(A));

for n=0:N-1

  t=n*dt;

  %- Schema
  switch SCHEMA
      case 'EE';
          % COMPLETER
          Tug_=Tug(t);
          Tud_=Tud(t);
          q_=q(t);
          m_=(Id - dt*(A-h*lambda*G))*P;
          m__=Id - dt*(A-h*lambda*G);
          P =  (Id - dt*(A-h*lambda*G))*P - dt*(q(t)+Tug(t)+Tud(t));
          
          % JCD : not yet
          %   case 'EI';
          %     % COMPLETER
          %     t1=t+dt;
          %     P = (Id + dt*A)\(P-dt*q(t1));
          %
      case 'CN';
          % COMPLETER
          q0=q(t);
          q1=q(t+dt);
          Tug0 = Tug(t);
          Tug1 = Tug(t+dt);
          Tud0 = Tud(t);
          Tud1 = Tud(t+dt);
          P = (Id+dt/2*(A-h*lambda*G)) \ ( (Id - dt/2*(A-h*lambda*G)) * P - dt/2*((q0+q1+Tug0+Tug1+Tud0+Tud1)));


  otherwise
    fprintf('SCHEMA not programmed'); abort;

  end

  if mod(n+1,deltan)==0; 	%- Printings at each deltan steps.

   %- Graphs at time t_{n+1}
   t1=(n+1)*dt; 
   ploot(t1,s,P); pause(1e-3);
 
   %- Error computations:
   %COMPLETER errLI
   Pex=Merton(t1,s,K,r,sigma,nMerton);		%- Merton
   errLI=norm(P-Pex,'inf');	%- Linfty error
   fprintf('t=%5.2f; iteration n=%4i; Err.Linf=%8.5f',t1,n+1,errLI);  
   fprintf('\n');
   %input('');
  end

end
t_total=toc();
fprintf('total time = %5.2f\n',t_total);
fprintf('program ended normaly');

