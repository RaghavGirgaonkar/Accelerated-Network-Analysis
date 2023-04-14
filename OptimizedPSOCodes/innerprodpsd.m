function innProd = innerprodpsd(fftY,fftXbyPSD)
%P = INNERPRODPSD(fftY,fftXbyPSD)
%Calculates the inner product of vectors X and Y for the case of Gaussian
%stationary noise having a specified power spectral density. fftY is the
%FFT of Y and fftXbyPSD is the FFT of vector X divided by Sn(f) and
%multiplied by the frequency amplitude term A. In this case, fftY is the fourier domain normalized
%frequency-phase term of the quadrature templates.
%The scalar factor (1/N) is derived from Parseval's theorem.

%Raghav Girgaonkar, Apr 2023

nSamples = length(fftXbyPSD);
dataLen = nSamples;
innProd = (1/dataLen)*(fftXbyPSD)*fftY';
innProd = real(innProd);
