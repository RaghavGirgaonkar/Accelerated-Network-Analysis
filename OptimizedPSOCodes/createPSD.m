%Script to create interpolated PSD from SHAPES and PWELCH Estimates

%% Provide Filename containing log-estimates for SHAPES and PWELCH
S = load('~/Documents/linetest_Job_col_PSD_dataSz4Up2N1024SNRR5e-03_C.mat');

%% Get One-sided estimates

welchPSD = S.dataY; 
shapesPSD = S.estS;

%% Data Parameters
Fs = 4096;
T = 4096;
N = Fs*T;

kNyq = floor(N/2);
fvec = (0:(kNyq))*Fs/N;

minidx = find(fvec<=15, 1, 'last' );
maxidx = find(fvec<=700, 1, 'last' );

nsamp = kNyq - minidx + 1;

%% Complete the PSD Vectors

winlen = 4; %Window length of Welch Estimate
n = winlen*Fs;

shapespsd_temp = [ones(1, n+1-length(shapesPSD))*shapesPSD(1), shapesPSD'];
welchpsd_temp = [ones(1, n+1-length(shapesPSD))*welchPSD(1), welchPSD'];

nyqfreq = n;
freqs = (0:nyqfreq)*(Fs/(2*n));

%% 1-D Interpolation

WPSD = interp1(freqs, welchpsd_temp, fvec);
SPSD = interp1(freqs, shapespsd_temp, fvec);

%% Antilog

SHAPESPSD = (10.^SPSD)./2;
WELCHPSD = (10.^WPSD)./2;