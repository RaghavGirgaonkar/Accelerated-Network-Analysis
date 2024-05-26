%% Script to demonstrate code optimization

%% Define Data Parameters
sampFreq = 4096;
datalen = 512;
N = datalen*sampFreq;
frange = [30,700];
fpos = (0:floor(N/2))*(1/datalen);
timeVec = (0:(datalen*sampFreq - 1))/sampFreq;

%% Create Colored Gaussian Noise
[noise, PSD] = LIGOnoise(N, sampFreq, 1, 'sample');
% S = load('realization_512s_inj138_fs4096_n1.mat');
% noise = S.dataY;
% PSD = S.dsstPSD;

%Create Total PSD and Transfer Function
negFStrt = 1-mod(N,2);
kNyq = floor(N/2)+1;

PSDtotal = [PSD, PSD((kNyq-negFStrt):-1:2)];
TFtotal = 1./sqrt(PSDtotal);


%Pre-calculate frequency magnitude and phase difference terms
 [A,avec, phaseDiff] = preprocessing(frange(1),frange(2),fpos, datalen, N);

%Create general normalization factor N_fW
AbysqrtPSD = A.*TFtotal;
innProd = (1/N)*(AbysqrtPSD)*AbysqrtPSD';
genNormfacSqr = real(innProd);
N_fW = 1/sqrt(genNormfacSqr);

%% Generate Signal to Inject
 chirptimes = [29.6373, 1.1045];
 ta = 138;
 initial_phase = 0;
 phase = 0;
 snr = 30;

 
 
 % Generate unnormalized signal vector of arbitrary SNR
fwavepos = waveform_tau(fpos,ta,phase,frange(1),frange(2),chirptimes(1),chirptimes(2),datalen,initial_phase, avec);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];
signal = ifft(A.*fwave);



%Normalize signal to preset SNR using expression for strain innerproduct
normFac = real(snr*(1/sqrt((1/(N*sampFreq))*(sum((fft(signal).*TFtotal).*conj(fft(signal).*TFtotal))))));
innProd = (1/N)*(fft(signal).*TFtotal)*(fft(signal).*TFtotal)';
normFacW = snr/sqrt(innProd);
% normFacW = real(snr*(1/sqrt((1/(N))*(sum((fft(signal).*TFtotal).*conj(fft(signal).*TFtotal))))));


signal = signal*normFac;
% signal = signal*normFacW;
% signal = signal*N_fW;

% whtndsignal = ifft(fft(signal).*TFtotal);

% strain = noise;

%Inject signal directly into strain
 strain = noise + signal;
%  strain = filtdata + signal;


%Highpassing (Optional)
% rolloff = 0.125;
% fmin = frange(1);
% segwin = strain.*tukeywin(length(strain),rolloff*sampFreq/length(strain))';
% seghpass = highpass(segwin, fmin, sampFreq, ImpulseResponse="iir",Steepness=0.95);
% strain = seghpass;

%Tukey Window data
rolloff = 0.5; %Roll-off in seconds
winstrain = strain.*tukeywin(length(strain), rolloff*sampFreq/N)';

% mf = (1/sampFreq)*ifft((fft(winstrain).*TFtotal).*conj(fft(qsignal).*TFtotal));



%% Whiten and Normalize Strain data

whtndstrain = ifft((1/sqrt(sampFreq))*fft(winstrain).*TFtotal);

% whtndstrain = whtndstrain + whtndsignal;

%% Perform Optimized Matched Filtering

%Create fftdatabyPSD vector
fftdatabyPSD = (A.*fft(whtndstrain).*TFtotal);


% N_fW = 1/sqrt((1/N)*sum((A.*TFtotal).*conj(A.*TFtotal)));

%Generate phase Fourier terms of quadrature templates q0,q1 using the same
%parameters used for signal injection
fwavepos = waveform_tau(fpos,0,phase,frange(1),frange(2),chirptimes(1),chirptimes(2),datalen,initial_phase, avec);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

phaseq0 = [fwavepos, fwaveneg];

phaseq1 = phaseq0.*phaseDiff; %Phase difference of pi/2

%Calculate M1 and M2
m1 = ifft(fftdatabyPSD.*conj(N_fW*phaseq0));
m2 = ifft(fftdatabyPSD.*conj(N_fW*phaseq1));

%Generate final timeseries
m = sqrt(m1.^2 + m2.^2);

%Calculate estimate Time of Arrival and SNR
[SNR,tindex] = max(m);
TOA = tindex/sampFreq;
% SNRs(i) = SNR;
% TOAs(i) = TOA;


% figure;
% histogram(SNRs,'NumBins',50);
% 
% 
% figure;
% histogram(TOAs,'NumBins',50);


% %Display parameters
dispstr = ['The estimated TOA = ',num2str(TOA),' and the estimated SNR = ',num2str(SNR)];
disp(dispstr);

%% Plots
%Plot m1 and m2 individually
figure;
plot(timeVec, m1);
xlabel('Time (s)');
title('M_1 Matched Filter Timeseries');

figure;
plot(timeVec, m2);
xlabel('Time (s)');
title('M_2 Matched Filter Timeseries');

figure;
plot(timeVec, m);
xlabel('Time (s)');
title('$\sqrt{M_1^2 + M_2^2}$ Matched Filter Timeseries','Interpreter','latex');













