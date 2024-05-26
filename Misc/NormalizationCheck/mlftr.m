function [mftimeseries] = mlftr(X,Y,PSD)
%MLFTR 
mftimeseries = (1/4096)*ifft((fft(X)./PSD).*conj(fft(Y)));
end

