function out=innerproduct(x,y,sampling_freq,psd)
%Function to calculate the inner product of x and y given by 
%<x,y> = ((X./S)*Y')/N
%N is the length in time, X and Y are the fourier transforms of x and y and S is the two sided
%Power spectral density (assuming the two sided PSD is provided as input)

nsamples = length(x);
%Length of signal in time
siglen = nsamples*sampling_freq;

X = fft(x);
Y = fft(y);

prod = (1/siglen)*(X./psd)*Y';
out = real(prod);

