function [outNoise, interPSD] = LIGOnoise(T_sig, num, Fs)
%Load PSD 
y = load('iLIGOSensitivity.txt','-ascii');
freqs = y(:,1);
sqrtPSD = y(:,2);

P = sqrtPSD.^2;

%Freq
% Fs = 2048;
% num=3;
% T_sig=54;

N = num*T_sig*Fs; %Total number of Time Samples

T = N/Fs;
timeVec = (0:(N-1))/Fs;

% fvec = freqs(1):(1/T):freqs(end);
fvec = 0:(1/T):(Fs/2);

%% Interpolation
interPSD = interp1(freqs,P, fvec);

noisePSD = interPSD;

%% Modifications for matched filtering PSD
minidx = find(fvec==15);
maxidx = find(fvec==700);

Sn30 = interPSD(minidx);
Sn700 = interPSD(maxidx);

interPSD(1:minidx) = Sn30;
interPSD(maxidx:end) = Sn700;


%% Modifications on PSD for making colored noise

minidx = find(fvec==15);
maxidx = find(fvec==700);

Sn15 = noisePSD(minidx);
Sn700 = noisePSD(maxidx);

noisePSD(1:minidx) = Sn15;
noisePSD(maxidx:end) = Sn700;

%% Make colored Noise
fltrOrdr = 500;
% 
outNoise_t = statgaussnoisegen(N+fltrOrdr,[fvec(:), noisePSD(:)],fltrOrdr,Fs);

outNoise = outNoise_t(fltrOrdr+1:end);




