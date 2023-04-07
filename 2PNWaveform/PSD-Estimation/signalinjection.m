function []=signalinjection(filename)
%% Script to inject 2PN signal in colored noise using the LIGO design sensitivity
%% Read JSON File 
addpath("../AllPSO/");
fname = filename;
str = fileread(fname);
filenames = jsondecode(str);
fname = filenames.signalparamfile;
str = fileread(fname);
params = jsondecode(str);
fname = filenames.psoparamfile;
str = fileread(fname);
pso = jsondecode(str);
fname = filenames.filenames;
str = fileread(fname);
files = jsondecode(str);
%Specify initial parameters
T_sig = params.signal.T_sig;
initial_phase = 0;
%% Sampling Frequency
Fs = params.sampling_freq;
%% Number of samples = num*Fs*T_sig
num = params.signal.num;
N = floor(num*T_sig*Fs);
timeVec = (0:N-1)*(1/Fs);
%% Min and Max Frequencies Hz
fmin = params.freq(1);
fmax = params.freq(2);
%% Positive Frequency Vector
datalen = N/Fs;
fpos = (0:floor(N/2))*(1/datalen);
%% Initial Time of Arrival and Coalescence Phase 
ta = params.ta;
phase = params.phase;
% Signal to noise ratio of the true signal
snr = params.snr;
%Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
% Mass parameters of the true signal
m1 = params.masses(1);
m2 = params.masses(2);

%% Get Noise timeseries and PSD vector

[colNoise, PSD] = LIGOnoise(T_sig, num, Fs);


%% Create Signal
wave = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, PSD);


%% Create final realization

data = wave + colNoise;

%% Plots (for testing)
figure;
hold on;
% plot(timeVec,data);
plot(timeVec, wave);
hold off;

%% Save in .mat file
% save('data-timeseries.mat','data');
% save('noise-timeseries.mat','colNoise');
end