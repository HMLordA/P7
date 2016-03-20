function y=descente(L,b)
%-- solving Ly=b: "descente"
n=length(b);
y=zeros(b);
y(1)=b(1)/L(1,1);
for k=2:n; y(k)=(b(k)-L(k,k-1)*y(k-1))/L(k,k); end


