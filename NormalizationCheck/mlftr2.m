function [mftimeseries] = mlftr2(X,Y,TFtotal)
%MLFTR on whitened X
mftimeseries = ifft((fft(X).*TFtotal).*conj(fft(Y)));
end

