function y=ul(t)
  global K r Smin
  %y=K*exp(-r*t)-Smin;	% European
  y=K-Smin; 		% American

