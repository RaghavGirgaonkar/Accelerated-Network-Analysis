function [innerProd] = innerproduct(X,Y,PSD, sampFreq)
%INNERPRODUCT 
% N = length(X);
N = sampFreq*length(X);
% Fs = 4096;
% N = length(X)/Fs;
innerProd = (1/N)*sum((fft(X)./PSD).*conj(fft(Y)));
innerProd = real(innerProd);
end

