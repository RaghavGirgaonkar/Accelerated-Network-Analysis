function timesVec = matchedfiltering(fftXbyPSD, fftY)
%function to perform FFT correlation based matched filtering using a signal
%y = signal + noise and q a unit normalised template starting at t = 0,
%assumes two-sided psd is provided
%returns the matchedfiltering timeseries

%Raghav Girgaonkar, April 2023

innProd = fftXbyPSD.*conj(fftY);
timesVec = ifft(innProd);

