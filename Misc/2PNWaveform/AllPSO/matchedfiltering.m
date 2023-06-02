function timesVec = matchedfiltering(xVec,yVec,~,psdVals)
%function to perform FFT correlation based matched filtering using a signal
%xVec = signal + noise and yVec a unit normalised template starting at t = 0,
%assumes two-sided psd is provided
%returns the matched filtering timeseries of the template yVec in the signal xVec 

%Raghav Girgaonkar, Apr 2023

nSamples = length(xVec);
if length(yVec) ~= nSamples
    error('Vectors must be of the same length');
end
kNyq = floor(nSamples/2)+1;
if length(psdVals) ~= kNyq
    error('PSD values must be specified at positive DFT frequencies');
end

fftX = fft(xVec);
fftY = fft(yVec);
%We take care of even or odd number of samples when replicating PSD values
%for negative frequencies
negFStrt = 1-mod(nSamples,2);
psdVec4Norm = [psdVals,psdVals((kNyq-negFStrt):-1:2)];

% dataLen = sampFreq*nSamples;
% innProd = (1/dataLen)*(fftX./psdVec4Norm)*fftY';
innProd = (fftX./psdVec4Norm).*conj(fftY);
timesVec = ifft(innProd);
