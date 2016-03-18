function res = Tug(t)
global p x lambda h ul I g
res = [zeros(I,1)];
for i=1:I;
    for j=(-p-i+1:-i);
        res(i) = res(i) + g(exp(x(j)))*exp(x(j));
    end;
    res(i) = (-h*lambda*ul(t))*res(i);
end;
end