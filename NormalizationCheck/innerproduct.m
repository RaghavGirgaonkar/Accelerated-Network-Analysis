function [innerProd] = innerproduct(X,Y,PSD)
%INNERPRODUCT 
N = length(X);
innerProd = (1/N)*sum((fft(X)./PSD).*conj(fft(Y)));

end

