function innProd = innerprodpsd(fftY,sampFreq,fftXbyPSD)
%P = INNERPRODPSD(X,Y,Fs,Sn)
%Calculates the inner product of vectors X and Y for the case of Gaussian
%stationary noise having a specified power spectral density. Sn is a vector
%containing PSD values at the positive frequencies in the DFT of X
%and Y. The sampling frequency of X and Y is Fs.

%Soumya D. Mohanty, Mar 2019

nSamples = length(fftXbyPSD);
dataLen = sampFreq*nSamples;
innProd = (1/dataLen)*(fftXbyPSD)*fftY';
innProd = real(innProd);
