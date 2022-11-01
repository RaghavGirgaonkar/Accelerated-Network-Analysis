%% Time of Signal based on Chirp Parameter Tau0 for both masses = 1.4 Solar Masses
% T_sig = Tau0 (in this case 54 seconds)
T_sig = 54;
%% Sampling Frequency
Fs = 2048;
%% Number of samples = num*Fs*T_sig
num = 3;
N = floor(num*T_sig*Fs);
timeVec = (0:N-1)*(1/Fs);
%% Min and Max Frequencies Hz
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

%% Create 2PN Waveform in Fourier Domain

fwavepos = waveform(fpos,ta,phase,fmin,m1,m2,0);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];
% fwave_t = fwave;
% 
% plot(real(fwave));
%% Set all indices to zero except for the ones corresponding to [fmin, fmax]
min_index  = floor(T_sig*fmin) + 1;
max_index = floor(T_sig*fmax) + 1;
% % 
fwave(1:min_index-1) = 0;
fwave(max_index+1: floor(N/2) + floor(N/2) - max_index + 1) = 0;
fwave(floor(N/2) + floor(N/2) - min_index + 3:end) = 0;
% plot(abs(fwave))
% plot(fftshift(abs(fwave)));
wave = ifft(fwave);
%% Plots
plot(timeVec,wave);
xlabel("Seconds");
title("2PN Waveform 30-700Hz")
% hold on;
% plot(timeVec,imag(wave));
%plot(fvec, fwave);



