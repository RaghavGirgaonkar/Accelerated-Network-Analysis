function timesVec = matchedfiltering(fftXbyPSD, fftY)
%Function to perform FFT correlation based matched filtering 
% Input: fftXbyPSD = FFT vector of data which has been multiplied by
%                    frequency magnitude vector A and divided by total PSD 
%        fftY = Normalized FFT phase vector of quadrature template
% Output: timesVec = MatchedFiltering timeseries vector

%Raghav Girgaonkar, April 2023

innProd = fftXbyPSD.*conj(fftY);
timesVec = ifft(innProd);

