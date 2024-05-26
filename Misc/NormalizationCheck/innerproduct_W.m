function [innerProd] = innerproduct_W(X,Y,PSD)
%INNERPRODUCT 
 N = length(X);
% N = 4096*length(X);
% Fs = 4096;
% N = length(X)/Fs;
innerProd = (1/N)*sum((fft(X)./PSD).*conj(fft(Y)));
innerProd = real(innerProd);
end

