function F = alphamuestm(x,C)
    F(1) = pochhammerA(x(1)+x(2),x(2)) / pochhammerA(x(1),x(2))-(1+1/C(1));
    F(2) = pochhammerA(x(1)+2*x(2),2*x(2)) / pochhammerA(x(1),2*x(2))-(1+1/C(2));
%     F(3) = x(1)^x(2) / pochhammerA(x(1),x(2)) * C(3) - x(3);
    F(3) = x(1)^x(2) / pochhammerA(x(1),x(2)) - x(3) / C(3);
%     F(3) = x(1)^x(2) - x(3) / C(3) * pochhammerA(x(1),x(2));
    
%     F(1) = gamma(x(1))/gamma(x(1)+  x(2)) * gamma(x(1)+2*x(2))/gamma(x(1)+  x(2))-1-1/C(1);
%     F(2) = gamma(x(1))/gamma(x(1)+2*x(2)) * gamma(x(1)+4*x(2))/gamma(x(1)+2*x(2))-1-1/C(2);
%     F(3) = x(1)^x(2) * gamma(x(1)) * C(3) / gamma(x(1)+x(2)) - x(3);
end