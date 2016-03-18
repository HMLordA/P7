function ploot(t,s,P)
global K r sigma
global Xmin Xmax Ymin Ymax
global Smin Smax

Pex=BS(t,s);

figure(1);
clf;
axis =[Xmin,Xmax,Ymin,Ymax];

sgraph  =[Smin;s;Smax];
Pgraph  =[ul(t);P;  ur(t)];
Pexgraph=[ul(t);Pex;ur(t)];

%- European case:
%plot(sgraph,P0(sgraph),'red--'); % Payoff function
%hold on;
%plot(sgraph,Pexgraph,'black.-'); 
%plot(sgraph,Pgraph,'blue.-');
%- American case:
PAYOFF=plot(sgraph,P0(sgraph),'blue.-','Linewidth',2); % Payoff function
hold on;
NUM=plot(sgraph,Pgraph,'black.-','Linewidth',2); 
legend([PAYOFF,NUM],'Payoff','Scheme','Location','Best');

titre=strcat('t=',num2str(t)); title(titre);
xlabel('s');
ylabel('price');
grid;

