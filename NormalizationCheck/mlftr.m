function [mftimeseries] = mlftr(X,Y,PSD)
%MLFTR 
mftimeseries = ifft((fft(X)./PSD).*conj(fft(Y)));
end

