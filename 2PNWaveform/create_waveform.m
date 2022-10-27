%% Time of Signal based on Chirp Parameter Tau0 for both masses = 1.4 Solar Masses
% T_sig = Tau0 (in this case 27 seconds)
T_sig = 27;
%% Sampling Frequency
Fs = 2048;
%% Number of samples = 3*Fs*T_sig
N = 3*T_sig*Fs; 
%% Min and Max Frequencies
fmin = 30;
fmax = 700;
%% Masses in Solar Mass 
m1 = 1.4;
m2 = 1.4;
%% Initial Time of Arrival and Phase 
ta = 0;
phase = 0;
%% Positive Frequency Vector
fpos = (0:floor(N/2))*(1/T_sig);
%% Negative Frequency Vector
fneg = (N -1 - floor(N/2):-1:1)*(-1/T_sig);
%% Total Frequency vector
fvec = [fpos, fneg];
%% Create 2PN Waveform in Fourier Domain
fwave = waveform(fvec,ta,phase,fmin,m1,m2);
plot(fwave);
%% Set all indices to zero except for the ones corresponding to [fmin, fmax]
min_index  = floor(T_sig*fmin) + 1;
max_index = floor(T_sig*fmax) + 1;

fwave(1:min_index-1) = 0;
fwave(max_index+1: floor(N/2) + floor(N/2) - max_index + 1) = 0;
fwave(floor(N/2) + floor(N/2) - min_index + 2:end) = 0;


wave = ifft(fwave);
%% Plots
plot(real(wave));
%plot(fvec, fwave);



