function []=animate_matchedfiltering(filename)
%% Read JSON File 
% addpath("../../PSO/");
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
%Tau coeffs as phase parameters
if pso.type == "tau"
    m1 = m1*Msolar;
    m2 = m2*Msolar;
    M = m1 + m2;
    u = m1*m2/(m1 + m2);
    n = u/M;
    tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
    tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);
    type = 1;
    disp("Tau Space PSO");
else
    type = 0;
    disp("Mass Space PSO");
end
% Search range of phase coefficients
if type
    rmin = [params.rmin_tau(1), params.rmin_tau(2)];
    rmax = [params.rmax_tau(1), params.rmax_tau(2)];
else
    rmin = [params.rmin(1), params.rmin(2)];
    rmax = [params.rmax(1), params.rmax(2)];
end
% Number of independent PSO runs
nRuns = pso.nruns;
%% Do not change below
% Generate data realization
% wgn = 
% dataX = (0:(nSamples-1))/Fs;
% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
% rng('default');
%% Generate Noise
% [noise,PSD] = LIGOnoise(N,Fs, params.signal.noise);
%% Load One Sided PSDs
% E = load(files.psdfile);
% estPSDs = E.estS';
% sz = size(estPSDs);
% PSD = resample(estPSDs, floor(N/2) + 1, sz(2));
%% Generate Colored Noise
% noise = genColNoise(PSD, Fs, [fmin,fmax], params.signal.noise);
noise = randn(1,N);
PSD = ones(size(fpos));
% Generate 2PN signal
if type
     wave = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0,tau1p5,datalen, initial_phase, snr, N,PSD);
else
     wave = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, PSD);
end

%% Generate Final Signal

dataY = wave + noise;

dataX = timeVec;

params = struct('dataX', dataX,...
                  'fpos', fpos,...
                  'dataY', dataY,...
                  'frange', [fmin,fmax],...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'rmin',rmin,...
                  'rmax',rmax,...
                  'psd',PSD,...
                  'Fs',Fs);
mftimeseries = zeros(1,N);

%Innerproduct of a template at a certain time
% Define the filename for the gif
giffilename = 'matchedfiltering_2.gif';

time_lapse = 0.0001;
title_text = 'Generation of Matched-Filtering Timeseries';

% Create a new figure window
fig = figure;


for i = 1:32:N
q0 = gen2PNwaveform_tau(params.fpos, i/Fs, 0, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.psd);
q1 = gen2PNwaveform_tau(params.fpos, i/Fs, pi/2, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.psd);

I0 = innerprodpsd(dataY, q0, PSD);
I1 = innerprodpsd(dataY, q1, PSD);

I = sqrt(I0^2 + I1^2);

mftimeseries(i) = I;


hold on;
plot(timeVec, dataY); hold on;
plot(timeVec, wave, Color='red' ); hold on;
plot(timeVec, q0,Color='yellow'); hold on;
plot(timeVec(1:i), mftimeseries(1:i), Color='black'); hold on;
ylim([-4,12]);
xlabel("Time (s)");
ylabel("Strain h(t)");

title(title_text,'FontSize',16,'Color',[0 0 0]), hold on;

drawnow;
%Capture frame
% Capture the current frame and save it to the gif
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if i == 1
        imwrite(imind,cm,giffilename,'gif', 'Loopcount',Inf,'DelayTime',time_lapse);
    else
        imwrite(imind,cm,giffilename,'gif','WriteMode','append','DelayTime',time_lapse);
    end

    clf(fig);

end

end