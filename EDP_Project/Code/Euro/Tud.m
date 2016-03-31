function res = Tud(t)
global p x lambda h ur I g Xmin
res = zeros(I,1);
for i=1:I;
    for j=(I-i+1:p+I-i);
        res(i) = res(i) + g(exp(x(j)-Xmin))*exp(x(j)-Xmin)*ur(t);
    end;
    res(i) = (-h*lambda)*res(i);
end;
end