%- DF - Black et Scholes 
%- European options
%- Jean-Cristophe DIETRICH & Nazar KOSTYUCHYK
clear

%------------------------
%- DONNEES FINANCIERES / FINANCIAL DATA
%------------------------
global  K r sigma T Smin Smax lambda mu gamma kappa
K=100; sigma=0.15; r=0.05; T=1; lambda = 0.1; mu = 0.0; gamma = 0.75; %Smin=30; Smax=500;

kappa = exp(mu+gamma^2/2)-1; % JCD : expectancy of eta, which is log-normal

%------------------------
%- DONNEES NUMERIQUES / NUMERICAL DATA
%------------------------
global I N p nMerton
I=160; N=160; p = 160; nMerton = 100;
%I=2*10; N=I*I/10; 

SCHEMA='CN-FFT'; 		%- 'EE' or 'EI' or 'CN' or 'EI-AMER-UL' or 'EI-AMER-NEWTON' or 'CN-AMER-UL' or 'CN-AMER-NEWTON' 'CN-FFT'
CENTRAGE='CENTRE'; 	%- 'CENTRE', 'DROIT', 'GAUCHE' 

%- Parameters for the graphics:
global Xmin Xmax Ymin Ymax
%Xmin=log(Smin/K); Xmax=log(Smax/K); Ymin=-20; Ymax=K;
Xmin=-2.0; Xmax=2.0; Ymin=-20; Ymax=K;
Smin=K*exp(Xmin);
Smax=K*exp(Xmax);
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
h=(Xmax-Xmin)/(I+1); %- pas d'espace  / mesh step
s=K*exp(Xmin+(1:I)'*h); 	%- maillage      / mesh : column vector of size I, containing the mesh values s_i = Smin + i*h

global x
x= @(i) Xmin + i*h;
%- IC : Initial Condition   = function  u0
%- BD : Boundary Conditions = functions ul, ur
%- ==> COMPLETE definition of functions u0, ul, ur (inline definitions: see below)
global ul ur g
u0= @(s) max(K-s,0);		%- Initial values (payoff function)
ul= @(t,i) K*exp(-r*t)-K*exp(x(i));	%- ul= left  value, at Smin        EUROP
%ul= @(t,i) K*exp(0*t)-K*exp(x(i)); %- ul= left  value, at Smin            AMERICAIN
ur= @(t) 0;			%- ur= right value, at Smax
g= @(eta) exp(-(log(eta)-mu)^2/(2*gamma^2))/(sqrt(2*pi)*gamma*eta) ;

%%%%%%%%%%%%%%%%
% A remplacer comme avant
%%%%%%%%%%%%%%%%

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
  bet=(r-lambda*kappa-sigma^2/2)/h;
  for i=1:I;   A(i,i) = 2*alpha + r + lambda; end;
  for i=2:I;   A(i,i-1) = -alpha + bet/2; end;
  for i=1:I-1; A(i,i+1) = -alpha - bet/2; end;
  
  for i=1:I;
      %G(i,i)=0;
      for j=1:I;
          G(i,j) = g(exp(x(j-i)-Xmin))*exp(x(j-i)-Xmin);
      end;
  end;
          

  % FILL IN
  q = @(t) [(-alpha + bet/2)*ul(t,0);  zeros(I-2,1);  (-alpha - bet/2)* ur(t)];
 % q = @(t) [-0.65*ul(t,0);  zeros(I-2,1);  (-alpha - bet/2)* ur(t)];
  
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
          
      case 'CN-FFT';
          % COMPLETER
          epss=0.01;
          q0=q(t);
          q1=q(t+dt);
          
          Vn1=P;
          Vn=P;
          
          %Xi=(((-p):p)'-p-1)*h;
          
          %Xi=((-p):p)'*h;
          %VG=exp(Xi);
          %for i=1:(2*p+1);
          %  VG(i) =VG(i)* g(exp(Xi(i)));
          %end;
          
          Xi=((-p):(p+I-1))'*h;
          VG=exp(Xi);
          for i=1:(2*p+I);
            VG(i) =VG(i)* g(exp(Xi(i)));
          end;
          %sum(VG)
          %FFTVG = transpose(fft(VG));
          FFTVG = fft(VG);
          %myFFTVG = FFTVG;
          %for l=2:I;
              %myFFTVG = [myFFTVG; FFTVG];
          %    myFFTVG = [myFFTVG FFTVG];
          %end;
          myMVn = zeros(2*p+I,1);
          for m=1:p;
              myMVn(m)=ul(t,m-p);
          end;
          for m=1:I;
              myMVn(m+p)=Vn(m);
          end;
          for m=1:p;
              myMVn(m+p+I)=ur(t);
          end;
              
          %VnConst=ifft(fft(MVn(Vn,t)).*conj(FFTVG));
          VnConst1=ifft(fft(myMVn).*conj(FFTVG));
          %VnConst = VnConst1((p+1):(p+I));
          VnConst = VnConst1((1):(I));
          %conjug=conj(myFFTVG);
          first=true;
          while(first|(max(abs(Vn1-Vn)./max(1,abs(Vn1))))>epss)
              first=false; 
              
              Vn=Vn1;
              %myMVn = MVn(Vn,t);
              myMVn = zeros(2*p+I,1);
              for m=1:p;
                  myMVn(m)=ul(t,m-p);
              end;
              for m=1:I;
                  myMVn(m+p)=Vn(m);
              end;
              for m=1:p;
                  myMVn(m+p+I)=ur(t);
              end;    
              %VnVar = ifft(fft(myMVn).*conj(myFFTVG)); 
              %VnVar = VnVar1((p+1):(p+I));
              
              %JCD : Convolution computation by FFT (comment the next 2 lines
              %to test the simple convolution)
              VnVar1 = ifft(fft(myMVn).*conj(FFTVG));  
              VnVar = VnVar1((1):(I));
              
              %JCD : Simple convolution computation (uncomment the next 2 lines
              %to test the simple convolution)
              %VnVar1 = conv(myMVn,VG);
              %VnVar = VnVar1((2*I+1):(3*I));
              
              %Vn1 = (Id+dt/2*A) \ ( (Id - dt/2*A) * P - dt/2*((q0+q1) - h*lambda*sum(transpose(VnVar),2) - h*lambda*sum(transpose(VnConst),2)) );
              Vn1 = (Id+dt/2*A) \ ( (Id - dt/2*A) * P - dt/2*((q0+q1) - h*lambda*VnVar - h*lambda*VnConst) );         
              
          end
          
          P=Vn1;
          
       case 'EI-AMER-UL';
          if n==0
              B=Id+dt*(A-h*lambda*G); [U,L]=uldecomp_sol(B);
              fprintf('Verification: norm(B-UL)=%10.5f\n', norm(B-U*L));
          end
          %- pb: min(Bx-b,x-g)=0, b=Pold-dt*q(t1), g=P0(s);
          % COMPLETE:
          t1=t+dt;
          Pold=P+dt*(q(t1)+Tug(t1)+Tud(t1));
          c=montee(U,P+dt*(q(t1)+Tug(t1)+Tud(t1)));
          P=descente_p(L,c,P0(s));
          
          %- Verification:
          err=norm(min(B*P-Pold,P-P0(s)));
          fprintf('Verification: err=%10.5f\n',err);
          
       case 'EI-AMER-NEWTON';
          if (n==0); B=Id+dt*(A-h*lambda*G); end;
          %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
          % COMPLETE 
          t1=t+dt;
          Pold=P;
          b=P+dt*(q(t1)+Tug(t1)+Tud(t1)); x0=P; gfinal=P0(s); eps=1e-10; kmax=50;
          [P,k]=newton_sol(B,b,gfinal,x0,eps,kmax);
          %- Verification
          err=norm(min(B*P-Pold,P-P0(s)));
          fprintf('Verif: err=%10.5f\n',err);
     
      case 'CN-AMER-UL';
          if n==0
              B=(Id+dt/2*(A-h*lambda*G)); [U,L]=uldecomp_sol(B);
              fprintf('Verification: norm(B-UL)=%10.5f\n', norm(B-U*L));
          end
          %- pb: min(Bx-b,x-g)=0, b=Pold-dt*q(t1), g=P0(s);
          % COMPLETE:
          t1=t+dt;
          q0=q(t1);
          q1=q(t1+dt);
          Tug0 = Tug(t1);
          Tug1 = Tug(t1+dt);
          Tud0 = Tud(t1);
          Tud1 = Tud(t1+dt);
          %(Id+dt/2*(A-h*lambda*G)) \ ( (Id - dt/2*(A-h*lambda*G)) * P - dt/2*((q0+q1+Tug0+Tug1+Tud0+Tud1)));
          Pold= (Id - dt/2*(A-h*lambda*G)) * P + dt/2*((q0+q1+Tug0+Tug1+Tud0+Tud1));
          b= (Id - dt/2*(A-h*lambda*G)) * P + dt/2*((q0+q1+Tug0+Tug1+Tud0+Tud1));
          c=montee(U,b);
          P=descente_p(L,c,P0(s));
          
          %- Verification:
          err=norm(min(B*P-Pold,P-P0(s)));
          fprintf('Verification: err=%10.5f\n',err);
          
      case 'CN-AMER-NEWTON';
          if (n==0); B=(Id+dt/2*(A-h*lambda*G)); end;
          %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
          % COMPLETE
          q0=q(t);
          q1=q(t+dt);
          Tug0 = Tug(t);
          Tug1 = Tug(t+dt);
          Tud0 = Tud(t);
          Tud1 = Tud(t+dt);
          
          Pold=P; %TODO: check why not in tp2
          b=(Id - dt/2*(A-h*lambda*G))*P + dt/2*((q0+q1+Tug0+Tug1+Tud0+Tud1)); 
          x0=P; gfinal=P0(s); eps=1e-10; kmax=50; %TODO : check value for x0 (coherent with b ?)
          % g is already a function name
          [P,k]=newton_sol(B,b,gfinal,x0,eps,kmax);
          %- Verification
          err=norm(min(B*P-Pold,P-P0(s)));
          fprintf('Verif: err=%10.5f\n',err);
    
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
   z=P-Pex;
   zz=norm(P-Pex,'inf');
   zzz=norm(z(4:length(z)),'inf');
   %input('');
  end

end
t_total=toc();
fprintf('total time = %5.2f\n',t_total);
fprintf('program ended normaly');

