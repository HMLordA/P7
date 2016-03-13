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
Pgraph=  [ul(t);P;  ur(t)];
Pexgraph=[ul(t);Pex;ur(t)];
PexgraphBS=[ul(t);PexBS;ur(t)];

plot(sgraph,Pexgraph,'black.-'); hold on; %exact
plot(sgraph,PexgraphBS,'green.-'); hold on; %exact
plot(sgraph,Pgraph,'blue.-'); %edp
%plot(sgraph,100*(Pgraph-Pexgraph),'red.-');
titre=strcat('t=',num2str(t)); title(titre);
xlabel('s');
ylabel('prix');
grid;

