function [innerProd] = innerproduct_optmzd(X,Y,PSD)
%INNERPRODUCT 
N = length(X);
innerProd = (1/N)*sum((X./PSD).*conj(Y));

end

