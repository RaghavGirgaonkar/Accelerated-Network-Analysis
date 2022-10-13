function normsig = normsignal(x, sampling_freq, psd, A)
%Function to normalise a signal to a specified snr = A
%Assuming two sided psd is provided

%nsamples = length(x);

normsigsqr = innerproduct(x,x,sampling_freq,psd);

N = A/sqrt(normsigsqr);

normsig = N*x;
