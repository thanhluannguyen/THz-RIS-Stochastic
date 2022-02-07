function out = pochhammerA(a,n)
    out = a^n*(1+(n-1)*n/(2*a)+(n-1)*n*(3*n^2-7*n+2)/24/a^2);
end