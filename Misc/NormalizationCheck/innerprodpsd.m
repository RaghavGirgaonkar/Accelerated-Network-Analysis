function innProd = innerprodpsd(fftY,fftXbyPSD)
%Calculates the inner product of vectors X and Y for the case of Gaussian
% stationary noise having a specified power spectral density (PSD) Sn(f). 
% Input: fftY = Normlized FFT vector of the phase term of quadrature template 
%        fftXbyPSD = FFT vector of data divided by PSD Sn(f) and
%                    multiplied by the frequency amplitude term A. 
% Output: innProd: inner product value 
%The scalar normalization factor (1/N) is derived from Parseval's theorem.

%Raghav Girgaonkar, Apr 2023

nSamples = length(fftXbyPSD);
dataLen = nSamples;
innProd = (1/dataLen)*(fftXbyPSD)*fftY';
innProd = real(innProd);
