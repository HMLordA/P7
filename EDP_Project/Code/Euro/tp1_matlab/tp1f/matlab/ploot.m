function ploot(t,s,P)
global K r sigma
global Xmin Xmax Ymin Ymax
global Smin Smax

Pex=BS(t,s);

figure(1);
clf;
axis =[Xmin,Xmax,Ymin,Ymax];

%- copying here the ul, ur, u0 functions
global ul ur u0

sgraph=[Smin;s;Smax];
Pgraph=  [ul(t);P;  ur(t)];
Pexgraph=[ul(t);Pex;ur(t)];


plot(sgraph,Pexgraph,'black.-'); hold on;
plot(sgraph,Pgraph,'blue.-');
%plot(sgraph,100*(Pgraph-Pexgraph),'red.-');
titre=strcat('t=',num2str(t)); title(titre);
xlabel('s');
ylabel('prix');
grid;

