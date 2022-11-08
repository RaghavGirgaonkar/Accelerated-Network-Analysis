%% Script to run matched filtering on 2PN Waveform
addpath ../MatchedFiltering/ 
%% Time of Signal based on Chirp Parameter Tau0 for both masses = 1.4 Solar Masses
% T_sig = Tau0 (in this case 54 seconds)
T_sig = 54;
initial_phase = 0;
%% Sampling Frequency
Fs = 2048;
%% Number of samples = num*Fs*T_sig
num = 3;
N = floor(num*T_sig*Fs);
timeVec = (0:N-1)*(1/Fs);
psd = ones(1,N);
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
datalen = N/Fs;
fpos = (0:floor(N/2))*(1/datalen);

%% Create q0 in Fourier Domain

fwavepos = waveform(fpos,ta,phase,fmin,fmax,m1,m2,datalen,initial_phase);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];

q0 = ifft(fwave);

%% Create waveform with certain SNR in Fourier Domain
snr = 10;

t = -80;

fwavepos = waveform(fpos,t,phase,fmin,fmax,m1,m2,datalen,initial_phase);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];

wave = ifft(fwave);

wave = snr*wave./norm(wave);

%% Add noise
n_iter = 1000;
SNRs = zeros(1,n_iter);

for j = 1:n_iter
final_wave = wave + randn(1,N);

[t_max, timesVec] = matchedfiltering(final_wave, q0, Fs, psd);

SNR = max(timesVec)/std(timesVec(1:10000));

SNRs(j) = SNR;

end

%% Plots
histogram(SNRs,50);
% plot(timeVec,timesVec, 'DisplayName','Likelihood Timeseries by Matched Filtering')
% plot(timeVec, final_wave,'DisplayName','Shifted Waveform + WGN');
% hold on;
% plot(timeVec,wave, 'DisplayName','Shifted Waveform with SNR = 10');
% % hold on;
% % plot(timeVec, q0);
% hold off;
legend
% xlabel("Seconds");
% title("2PN Waveform 30-700Hz")