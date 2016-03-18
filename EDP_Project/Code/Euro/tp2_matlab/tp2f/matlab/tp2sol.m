%- Finite Difference scheme for Black and Scholes :
%- American options.
%- Nov 2012, Olivier Bokanowski
clear

%------------------------
%- FINANCIAL DATA
%------------------------
global  K r sigma T Smin Smax
K=100; sigma=0.2; r=0.1; T=1;  Smin=0; Smax=200;  

%- CDI: Initial Condition
%- ==> COMPLETE the definition of P0.m, ul.m, ur.m



%------------------------
%- NUMERICAL DATA
%------------------------
N=1000;
I=50-1;

% choose SCHEME value in: {'EE' , 'EI' , 'CN' , 'EE-AMER' , 'EI-AMER-UL' , 'EI-AMER-NEWTON' , 'EI-AMER-SPLITTING'}
SCHEME='EE-AMER';
%SCHEME='EI-AMER-UL'; 
%SCHEME='EI-AMER-NEWTON'; 
%SCHEME='EI-AMER-SPLITTING';
CENTRAGE='CENTRE'; % 'CENTRE', 'DROIT', or 'GAUCHE' 


%- Paramètres pour la fenetre graphique:
global Xmin Xmax Ymin Ymax
Xmin=Smin; Xmax=Smax; Ymin=-20; Ymax=K; 
err_scale=0; %- Echelle pour le graphe d'erreur.
deltan=N/10; %- Eventuellement, Affichage uniquement tous les deltan pas.


%- Affichage des données:
fprintf('sigma=%5.2f, r=%5.2f, Smax=%5.2f\n',sigma,r,Smax);
fprintf('Mesh I= %5i, N=%5i\n',I,N);
fprintf('CENTRAGE : %s\n',CENTRAGE)
fprintf('SCHEME: %s\n',SCHEME)

%-----------------------------------------------
%- FORMULE de Black et Scholes.
%-----------------------------------------------
%- Formule de Black et scholes pour le put europeen vanille : fichier BS.m
%- function y=BS(t,s,K,r,sigma)


%--------------------
%- Maillage
%--------------------
%- COMPLETER: dt, h (pas d'espace), s (maillage: vecteur colonne), 
dt=T/N; 		%- pas de temps 
h=(Smax-Smin)/(I+1); 	%- pas d'espace = s_{i+1}-s_i
s=Smin+(1:I)'*h; 	%- vecteur colonne de taille i contenant les s_i

%- COEFFICIENT CFL
%COMPLETER
cfl=dt/h^2 * (sigma*Smax)^2;
fprintf('CFL : %5.3f\n',cfl); 


%//--------------------------
%//- Initialiser la matrice A
%//--------------------------
%//COMPLETER:  matrice A , fonction q(t)

switch CENTRAGE

case 'CENTRE';

  %COMPLETER 
  %alpha(i)=
  %bet(i)=
  %A(i,i), A(i,i-1), A(i,i+1)
  A=zeros(I,I);
  alpha=sigma^2/2 * s.^2 /h^2;
  bet=r*s/(2*h);
  for i=1:I;   A(i,i) = 2*alpha(i) + r; end;
  for i=2:I;   A(i,i-1) = -alpha(i) + bet(i); end;
  for i=1:I-1; A(i,i+1) = -alpha(i) - bet(i); end;

  %q=inline('[(-alpha(1) + bet(1))* ul(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* ur(t)]','t','alpha','bet','I');
  q= @(t) [(-alpha(1) + bet(1))* ul(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* ur(t)];


case 'DROIT' ;

  A=zeros(I,I);
  alpha=sigma^2/2 * s.^2 /h^2;
  bet=r*s/(h);
  for i=1:I;   A(i,i) = 2*alpha(i) + bet(i) + r; end;
  for i=2:I;   A(i,i-1) = -alpha(i) ; end;
  for i=1:I-1; A(i,i+1) = -alpha(i) - bet(i); end;

  q=inline('[(-alpha(1))* ul(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* ur(t)]','t','alpha','bet','I');

case 'GAUCHE';

  A=zeros(I,I);
  alpha=sigma^2/2 * s.^2 /h^2;
  bet=r*s/(h);
  for i=1:I;   A(i,i)   =2 *alpha(i) - bet(i) + r; end;
  for i=2:I;   A(i,i-1) =  -alpha(i) + bet(i) ; end;
  for i=1:I-1; A(i,i+1) =  -alpha(i) ; end;
  q=inline('[(-alpha(1)+bet(1))* ul(t);  zeros(I-2,1);  (-alpha(end))* ur(t)]','t','alpha','bet','I');

otherwise 

  fprintf('CENTRAGE not programmed !'); break

end
 
%--------------------
%- Initialiser P et graphique
%--------------------
P=P0(s);
ploot(0,s,P);
fprintf('appuyer sur la touche Return'); input('');


%--------------------
%- BOUCLE PRINCIPALE
%--------------------
%- demarrage compteur temps
tic(); 

Id=eye(size(A));

for n=0:N-1

  t=n*dt;

  %- Schema
  switch SCHEME 
  case 'EE'; 
    P =  (Id - dt*A)*P - dt*q(t);

  case 'EI'; 
    t1=t+dt; 
    P = (Id + dt*A)\(P-dt*q(t1));

  case 'CN';
    q0=q(t   );
    q1=q(t+dt);
    P = (Id + dt/2*A) \ ( (Id - dt/2*A) * P - dt*(q0+q1)/2 );
 
  case 'EE-AMER'
    % COMPLETE
    P = max(P0(s), (Id - dt*A)*P - dt*q(t));

  case 'EI-AMER-UL';
    if n==0
      B=Id+dt*A; [U,L]=uldecomp_sol(B);
      fprintf('Verification: norm(B-UL)=%10.5f\n', norm(B-U*L));
    end
    %- pb: min(Bx-b,x-g)=0, b=Pold-dt*q(t1), g=P0(s);
    % COMPLETE:
    t1=t+dt;
    Pold=P-dt*q(t1);
    c=montee(U,P-dt*q(t1));
    P=descente_p(L,c,P0(s));

    %- Verification:
    err=norm(min(B*P-Pold,P-P0(s)));
    fprintf('Verification: err=%10.5f\n',err);

  case 'EI-AMER-NEWTON';
    if (n==0); B=Id+dt*A; end;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    % COMPLETE 
    Pold=P;
    b=P; x0=P; g=P0(s); eps=1e-10; kmax=50;
    [P,k]=newton_sol(B,b,g,x0,eps,kmax);
    %- Verification
    err=norm(min(B*P-Pold,P-P0(s)));
    fprintf('Verif: err=%10.5f\n',err);

  case 'EI-AMER-SPLITTING'
    % COMPLETE
    t1=t+dt; 
    P = max(P0(s), (Id + dt*A)\(P-dt*q(t1)));

  otherwise
    fprintf('SCHEME not programmed'); abort;

  end
 
  if mod(n+1,deltan)==0; %- Affichage tous les deltan pas.

   %- Graphiques:
   t1=(n+1)*dt; 
   ploot(t1,s,P); pause(1e-3);
 
   %- Calculs d'erreurs:
   %COMPLETER
   Pex=BS(t1,s);   %- Appel de la solution Black et Scholes
   errLI=norm(P-Pex,'inf');        %- Calcul erreur Linfty
   fprintf('t=%5.2f; Err.Linf=%8.5f',t1,errLI);  
   fprintf('\n');
   %input('');
  end

end
t_total=toc();
fprintf('total time = %5.2f\n',t_total);
fprintf('program ended normaly');

