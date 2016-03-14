%- DF - Black et Scholes 
clear

%------------------------
%- DONNEES / FINANCIAL DATA
%------------------------
global  K r sigma T Smin Smax
K=100; sigma=0.2; r=0.1; T=1;  Smin=0; Smax=200;  

%- IC : Initial Condition   = function  u0
%- BD : Boundary Conditions = functions ul, ur
%- ==> COMPLETE definition of functions u0, ul, ur (inline definitions: see below)
global ul ur u0
u0= @(s) 0.0;		%- Initial values (payoff function) - typically u0=@(s) max(K-s,0)
ul= @(t) 0.0;		%- ul= left  value, at Smin
ur= @(t) 0.0;		%- ur= right value, at Smax


%------------------------
%- DONNEES NUMERIQUES / NUMERICAL DATA
%------------------------
I=10;
N=10;

SCHEMA='EE'; 		%- 'EE' or 'EI' or 'CN' 
CENTRAGE='CENTRE'; 	%- 'CENTRE', 'DROIT', 'GAUCHE' 

%- Parameters for the graphics:
global Xmin Xmax Ymin Ymax
Xmin=Smin; Xmax=Smax; Ymin=-20; Ymax=K;
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
dt=0; 		%- pas de temps  / time step
h=0; 		%- pas d'espace  / mesh step
s=0;  	 	%- maillage      / mesh : column vector of size I, containing the mesh values s_i = Smin + i*h, i=1,...,I

%- CFL COEFFICIENT 
%COMPLETE
cfl=0;
fprintf('CFL : %5.3f\n',cfl); 


%--------------------------
%- Initializations:  matrix A, function q(t)
%--------------------------

switch CENTRAGE

case 'CENTRE';  %- CENTERED APPROXIMATION

  %- FILL IN / COMPLETER matrice A, vecteurs alpha et bet de taille I, et fonction q
  A=zeros(I,I);
  %FILL the values of A(i,i), A(i,i-1), A(i,i+1)

  % FILL IN
  q = @(t) [0;  zeros(I-2,1);  0]; %- COMPLETE CORRECTLY FIRST AND LAST VALUES

case 'DROIT';	%- FORWARD DIFFERENCES

  fprintf('not programmed !'); 

case 'DROIT';	%- BACKWARD DIFFERENCES

  fprintf('not programmed !'); 

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

for n=0:N-1

  t=n*dt;

  %- Schema
  switch SCHEMA 
  case 'EE'; 
    % COMPLETER
    P = P;

  case 'EI'; 
    % COMPLETER
    P = P;

  case 'CN';
    % COMPLETER
    P = P;

  otherwise
    fprintf('SCHEMA not programmed'); abort;

  end

  if mod(n+1,deltan)==0; 	%- Printings at each deltan steps.

   %- Graphs at time t_{n+1}
   t1=(n+1)*dt; 
   ploot(t1,s,P); pause(1e-3);
 
   %- Error computations:
   %COMPLETER errLI
   Pex=BS(t1,s);		%- Black and Scholes
   errLI=0.0; 			%- Linfty error
   fprintf('t=%5.2f; iteration n=%4i; Err.Linf=%8.5f',t1,n+1,errLI);  
   fprintf('\n');
   %input('');
  end

end
t_total=toc();
fprintf('total time = %5.2f\n',t_total);
fprintf('program ended normaly');

