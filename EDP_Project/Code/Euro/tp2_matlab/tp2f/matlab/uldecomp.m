function [U,L]=uldecomp(A)
%---------------------------------------------------
%- DECOMPOSITION  A = U * L
%- case A= tridiag(.,.,.); 
%- with L=tridiag(l,d,0), U=tridiag(0,1,u)
%---------------------------------------------------
n=size(A,1);
d=zeros(n,1); 
l=zeros(n-1,1);
u=zeros(n-1,1);
%- TO COMPLETE: diagonal terms, lower diagonal, upper diagonal
%-   d(1) u(1)   0  0    ...
%-   l(1) d(2) u(2) 0    ...
%-     0  l(2) d(3) u(3) ...
%-    ...
%d(n) =
%for i=n-1:-1:1
%  l(i)=
%  u(i)=
%  d(i)=
%end
L=diag(d)+diag(l,-1);
U=eye(size(A))+diag(u,1);

