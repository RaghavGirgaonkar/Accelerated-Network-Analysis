function []=fitValimg(filename)
%% Script to create an image of the Fitness Values in Tau Space
%% Read JSON File 
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



%% Grid for the Image
ngridpoints = params.fitnessimage.ngridpoints;

tau0_lin = linspace(params.fitnessimage.tau0_range(1), params.fitnessimage.tau0_range(2), ngridpoints);
tau1p5_lin = linspace(params.fitnessimage.tau1p5_range(1), params.fitnessimage.tau1p5_range(2), ngridpoints);

nsectors = params.fitnessimage.nsectors;
nindex = ngridpoints/nsectors;

fitvals = zeros(nindex);

% rng('default');
noise = load("noise_realizations.mat");
wgn = noise.wgn(params.signal.noise,1:end);

dataX = timeVec;

wave = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0,tau1p5,datalen, initial_phase, snr, N);
dataY = wave + wgn;
inParams = struct('dataX', dataX,...
                          'fpos', fpos,...
                          'dataY', dataY,...
                          'frange', [fmin,fmax],...
                          'datalen',datalen,...,
                          'initial_phase', initial_phase,...
                          'N', N,...
                          'rmin',rmin,...
                          'rmax',rmax);

startidx_tau0 = params.fitnessimage.startidx_tau0;
startidx_tau1p5 = params.fitnessimage.startidx_tau1p5;

parfor i = 1:nindex
    
    for j = 1:nindex

         fitvals(i,j) = -1*mfqc_tau([tau0_lin(startidx_tau0+i-1), tau1p5_lin(startidx_tau1p5+j-1)], inParams);
    end
     

end

save(files.fitvalimgfile,"fitvals");



end

