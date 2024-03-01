function [outNoise, PSD] = LIGOnoise(N, Fs, noise_num, noisefile)
%Function to create colored noise using LIGO Design Sensitivities 
% Design PSD is modified between 15 Hz and 700Hz.
% Input: N = Total number of samples,
%        Fs = Sampling Frequency,
%        (Optional) noise_num = noise realization number from a pre-created noise realizations file
%        (Optional) noisefile = pre-created noise realizations filename
% Output: outNoise = colored noise vector,
%         PSD = two-sided PSD vector for positive DFT frequencies

% Raghav Girgaonkar, April 2023

%Load PSD 
y = load('iLIGOSensitivity.txt','-ascii');
% freqs = y(:,1);
% sqrtPSD = y(:,2);

%Turn one-sided sensitivity to two-sided
y(:,2) = (1/sqrt(2))*y(:,2);

% Interpolate sensitivity curve to positive DFT frequencies
minF = min(y(:,1));
maxF = max(y(:,1));
if minF ~= 0
% f=0 does not exist, put it in
y = [0, y(1,2);...
                  y];
end
if maxF < Fs/2
    error('High frequency limit requested is higher than supplied');
end


%Positive DFT frequencies
kNyq = floor(N/2)+1;
fvec = (0:(kNyq-1))*Fs/N;

%% Interpolation
interPSD = interp1(y(:,1),y(:,2), fvec);

%% Modifications, change cutoff frequencies as needed 
minidx = find(fvec<=30, 1, 'last' );
maxidx = find(fvec<=700, 1, 'last' );

Sn50 = interPSD(minidx);
Sn700 = interPSD(maxidx);
 
interPSD(1:minidx) = Sn50;
interPSD(maxidx:end) = Sn700;

%Convert Amplitude spectral density to Power Spectral Density
PSD = interPSD.^2;

%% Make colored Noise
fltrOrdr = 10000;

outNoise_t = statgaussnoisegen(N,[fvec(:),PSD(:)],fltrOrdr,Fs, noise_num, noisefile);

outNoise = outNoise_t(fltrOrdr+1:end - fltrOrdr);




