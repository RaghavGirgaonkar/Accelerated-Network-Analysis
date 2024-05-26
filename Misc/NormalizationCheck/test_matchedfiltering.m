% Define Parameters
sampFreq = 4096;
datalen = 512;
N = datalen*sampFreq;
frange = [30,700];
fpos = (0:floor(N/2))*(1/datalen);
timeVec = (0:(datalen*sampFreq - 1))/sampFreq;

%Create Colored Gaussian Noise
[noise, PSD] = LIGOnoise(N, sampFreq, 1, 'sample');

%Create Total PSD
negFStrt = 1-mod(N,2);
kNyq = floor(N/2)+1;

PSDtotal = [PSD, PSD((kNyq-negFStrt):-1:2)];

%Generate signal normalized to specified SNR
 chirptimes = [29.6373, 1.1045];
 ta = 138;
 initial_phase = 0;
 phase = 0;
 snr = 30;

 [A,avec, phaseDiff] = preprocessing(frange(1),frange(2),fpos, datalen, N);

 signal = gen2PNtemplate_tau(fpos, ta, phase, frange(1), frange(2),chirptimes(1), chirptimes(2),datalen,initial_phase,snr, N, A, avec, PSDtotal);
 
    signal = signal*sqrt(sampFreq);
%  signal = signal/norm(signal);

 figure;
 plot(timeVec,signal);

 %Test 1: Inject signal directly into strain data

straindata = noise + signal;

  % Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', straindata,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_new(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test1'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test1'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test1'); xlabel('t');


%Test 2: Inject signal into strain data and whiten together without
%normalizing whitened strain to unit variance

[whtndseg,wstdv, TFtotal] = segdatacond(straindata, PSD, sampFreq, [1,20*sampFreq]);

% whtndseg = whtndseg*wstdv;

% Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', whtndseg,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'TFtotal',TFtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_whtnd_og(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test2'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test2'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test2'); xlabel('t');


%Test 3: Inject signal into strain data and whiten together with
%normalizing whitened strain to unit variance

[whtndseg,~, TFtotal] = segdatacond(straindata, PSD, sampFreq, [1,20*sampFreq]);


% Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', whtndseg,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'TFtotal',TFtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_whtnd_og(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test3'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test3'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test3'); xlabel('t');


%Test 4: Inject a whitened signal into a whitened data segment which has
%been normalized to unit variance.

%First generate normalized whiten data realization
[whtndseg,~, TFtotal] = segdatacond(noise, PSD, sampFreq, [1,20*sampFreq]);

 whtndseg = whtndseg;

%Generate whitened signal
whtndsig = ifft(fft(signal).*(TFtotal));

%Add whitened signal to normalized whitened data
whtndseg = whtndseg + whtndsig;

% Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', whtndseg,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'TFtotal',TFtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_whtnd_og(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test4'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test4'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test4'); xlabel('t');

%Test 5, add injection directly into strain data, but scaled by the stdv of
%the whitened data.

%Get stdv
[~,wstdv,~] = segdatacond(noise, PSD, sampFreq, [1,200*sampFreq]);

% Create signal and scale 

sigscaled = signal*wstdv;

%Add to strain

straindata = noise + sigscaled;

  % Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', straindata,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_new(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test5'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test5'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test5'); xlabel('t');



%Test 6, matched filtering on whitened data normalized to unit variance:


%First generate normalized whiten data realization
[whtndseg,~, TFtotal] = segdatacond(noise, PSD, sampFreq, [1,datalen*sampFreq]);

% Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', whtndseg,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'TFtotal',TFtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_whtnd_og(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test6'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test6'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test6'); xlabel('t');

%Test 7: Inject signal directly into the strain data scaled by the whtnd
%stdv

[~,wstdv, ~] = segdatacond(noise, PSD, sampFreq, [1,datalen*sampFreq]);

straindata = noise + wstdv*signal;

[whtndseg,~, TFtotal] = segdatacond(straindata, PSD, sampFreq, [1,20*sampFreq]);

 % Input Parameters Structure:
inParams = struct('fpos', fpos,...
                  'dataY', whtndseg,...
                  'frange', frange,...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'avec', avec,...
		          'T_sig', 54,...
                  'PSDtotal',PSDtotal,...
                  'TFtotal',TFtotal,...
                  'Fs',sampFreq);

[mf1,mf2,max_val,max_arg] = mfgw_tau_whtnd_og(chirptimes, inParams);
mftimeseries = sqrt(mf1(1:end - 54*sampFreq).^2 + mf2(1:end - 54*sampFreq).^2);

figure;
plot(timeVec,mf1); title('MF1 Test7'); xlabel('t');

figure;
plot(timeVec,mf2); title('MF2 Test7'); xlabel('t');

figure;
plot(timeVec(1:end - 54*sampFreq),mftimeseries); title('MFTimeseries Test7'); xlabel('t');


%Test 8: Repeat Test 7 100 times and get histogram of recovered SNR

snrs = [];
for i = 1:100
    disp(i);
    [noise, ~] = LIGOnoise(N, sampFreq, 1, 'sample');
    [~,wstdv, ~] = segdatacond(noise, PSD, sampFreq, [1,datalen*sampFreq]);

    straindata = noise + wstdv*signal;

    [whtndseg,~, TFtotal] = segdatacond(straindata, PSD, sampFreq, [1,20*sampFreq]);
    inParams.dataY = whtndseg;
    [~,~,max_val,~] = mfgw_tau_whtnd_og(chirptimes, inParams);
    snrs = [snrs, sqrt(-max_val)];
end

figure;
histogram(snrs, 50);
hold on;
xline(snr);
hold off;






