%% Script to generate a 2PN waveform normalized to a given SNR
function [signal]=gensignal()
    %Injection parameters
    masses = [4,300]; %Change to chirp-times as needed
    initial_phase = 0; 
    phase = 0; %Coalescence phase of the signal
    ta = 138; %TOA
    datalen = 512; %Length of signal in seconds
    frange = [30,700]; %Frequency range
    snr = 15; % Injection SNR
    sampFreq = 4096; %Sampling Frequency
    N = datalen*sampFreq; %Total number of samples
    fpos = (0:floor(datalen*sampFreq/2))*(1/datalen); %Positive Fourier frequency vector
    
    %PSD (Uncomment if using different PSD for normalization)
    %Get two-sided PSD from iLIGO sensitivities
    [~, PSD] = LIGOnoise(N,sampFreq,1,'sample');
    
    %Make total PSD vector
    negFStrt = 1-mod(N,2);
    kNyq = floor(N/2)+1;
    PSDtotal = [PSD, PSD((kNyq-negFStrt):-1:2)];
    
    
    %Preprocessing Vectors
    [A,avec, phaseDiff] = preprocessing(frange(1), frange(2), fpos, datalen, datalen*sampFreq);
    
    
    %Generate signal
    %Comment below line and uncomment chirp time line if using chirptimes
     signal = gen2PNtemplate_mass(fpos, ta, phase, frange(1), frange(2),masses(1), masses(2),datalen,initial_phase,snr,N, A, avec, PSDtotal);
    
    
    %Uncomment if using chirptimes to create signal (Make sure to then comment
    %above line)
    % signal = gen2PNtemplate_tau(fpos, ta, phase, frange(1), frange(2),masses(1), masses(2),datalen,initial_phase,snr,N, A, avec, PSDtotal);
    
    %Uncomment if using negative chirp times to create waveform
%     signal = gen2PNtemplate_tau_negative(fpos, ta, phase, frange(1), frange(2),masses(1), masses(2),datalen,initial_phase,snr,N, A, avec, PSDtotal);
end