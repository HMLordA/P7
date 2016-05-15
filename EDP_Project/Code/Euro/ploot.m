function ploot(t,s,P)
global K r sigma nMerton
global Xmin Xmax Ymin Ymax
global Smin Smax

PexBS = BS(t,s,K,r,sigma);
Pex=Merton(t,s,K,r,sigma,nMerton);

figure(1);
clf;
axis =[Xmin,Xmax,Ymin,Ymax];

%- copying here the ul, ur, u0 functions
global ul ur u0

sgraph=[Smin;s;Smax];
Pgraph=  [ul(t,0);P;  ur(t)];
Pexgraph=[ul(t,0);Pex;ur(t)];
PexgraphBS=[ul(t,0);PexBS;ur(t)];
PexgraphINTRINSIC=[ul(t,0); max(K-s,0) ;ur(t)]; 
PexgraphError=[ul(t,0);Pex-PexBS;ur(t)];

plot(sgraph,Pexgraph,'magenta.--'); hold on; %exact
plot(sgraph,PexgraphBS,'green.-.'); hold on; %exact
plot(sgraph,Pgraph,'blue.-'); hold on; %edp
legend ('\fontname{Times} \fontsize{12} Exact Merton','\fontname{Times} \fontsize{12} Theorical BS','\fontname{Times} \fontsize{12} Explicit Euler');
%legend ('\fontname{Times} \fontsize{12} Exact Merton','\fontname{Times} \fontsize{12} Theorical BS','\fontname{Times} \fontsize{12} Crank-Nicolson');
%legend ('\fontname{Times} \fontsize{12} Exact Merton','\fontname{Times} \fontsize{12} Theorical BS','\fontname{Times} \fontsize{12} Crank-Nicolson FFT');
%-------------------
%plot(sgraph,PexgraphINTRINSIC,'red.-','Linewidth',2); %intrinsic value
%legend ('\fontname{Times} \fontsize{12} Intrinsic Value','\fontname{Times} \fontsize{12} EI Brennan Schwartz');
%legend ('\fontname{Times} \fontsize{12} Intrinsic Value','\fontname{Times} \fontsize{12} EI Newton');
%legend ('\fontname{Times} \fontsize{12} Intrinsic Value','\fontname{Times} \fontsize{12} CN Brennan Schwartz');
%legend ('\fontname{Times} \fontsize{12} Intrinsic Value','\fontname{Times} \fontsize{12} CN Newton');
%-------------------
%plot(sgraph,100*(Pgraph-Pexgraph),'red.-');
titre=strcat('t=',num2str(t)); title(titre);
xlabel('S');
ylabel('Price');

grid;

