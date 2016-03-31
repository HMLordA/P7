function res = Tug(t)
global p x lambda h ul I g Xmin
res = [zeros(I,1)];
for i=1:I;
    for j=(-p-i+1:-i);
        res(i) = res(i) + g(exp(x(j)-Xmin))*exp(x(j)-Xmin)*ul(t,i+j);
    end;
    res(i) = (-h*lambda)*res(i);
end;
end