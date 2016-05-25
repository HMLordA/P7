function M = MVn(Vn,t)
%function creating the Un vector we consider (with correct boundary terms)
global p x lambda h ul ur I g Xmin
M = zeros(2*p+1,I);
for i=1:(2*p+1);
    for j=(1:I);
        if ((i+j-p-1)<=0);
            M(i,j) = ul(t,i-j-p+1);
        elseif ((i+j-1)>(I+p));
            M(i,j) = ur(t);
        else
            M(i,j)=Vn(i+j-p-1);
        end;
    end;
end;
end