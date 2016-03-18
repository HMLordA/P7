function res = Tud(t)
global p x lambda h ur I g
res = zeros(I,1);
for i=1:I;
    for j=(I-i+1:p+I-i);
        res(i) = res(i) + g(exp(x(j)))*exp(x(j));
    end;
    res(i) = (-h*lambda*ur(t))*res(i);
end;
end