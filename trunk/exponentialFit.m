function y = exponentialFit(x, t)

numCoefs = length(x);
y = x(1);

for i = 2:(numCoefs-1)
    y = y + ( x(i).*exp((-x(i+1).*t)) );        
end