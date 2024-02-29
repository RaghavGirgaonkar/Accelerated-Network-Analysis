function [PSD, fvec]=createPSD(sampFreq, Tsig, welchPSD, freqs)
%% Script to create interpolated PSD from SHAPES and PWELCH Estimates

%% Data Parameters
Fs = sampFreq;
T = Tsig;
N = Fs*T;

kNyq = floor(N/2);
fvec = (0:(kNyq))*Fs/N;

% nyqfreq = winlen*Fs;

% freqs = (0:nyqfreq)*(Fs/(2*nyqfreq));

%% 1-D Interpolation
logwelchPSD = log10(welchPSD);
loginterPSD = interp1(freqs, logwelchPSD, fvec);

% %% Antilog

PSD = (10.^loginterPSD);
% PSD = (loginterPSD)./2;
