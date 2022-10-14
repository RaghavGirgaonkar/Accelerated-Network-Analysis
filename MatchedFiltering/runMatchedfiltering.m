%Script to run matchedfiltering with WGN

%Specify initial parameters
%Time length of signal and total time (in seconds)
T_sig = 2;
T_total = 16;
%Desired SNR
A = 10;
%Parameters of the Quadratic Chirp
a = 10;
b = 3;
c = 3;
coeffs = [a,b,c];
%Sampling frequency (assumed 4 times that of the nyquist freq in specified QC)
nyquist_freq = 2*(coeffs(1)*T_sig + coeffs(2)*T_sig.^2 + coeffs(3)*T_sig.^3);
sampling_freq = 4*nyquist_freq;
sampling_interval = 1/sampling_freq;
%Time vectors for signal and total series
timeVecSig = 0:sampling_interval:T_sig;
timeVecTot = 0:sampling_interval:T_total;
%Number of samples
nsamples_sig = length(timeVecSig);
nsamples_tot = length(timeVecTot);
%Generate PSD vector (assumed constant in this special case)
psd = ones(1,nsamples_tot);
%Generate signal
signal = genqc(timeVecSig,A,coeffs);
final_signal = [signal, zeros(1,(nsamples_tot-nsamples_sig))];
%Plot signal 
%plot(timeVecTot, final_signal);
%Generate White Gaussian Noise
wgn = randn(1, nsamples_tot);
% plot(timeVecTot, wgn);
% hold on;
% plot(timeVecTot, final_signal, 'r');
%Shift the signal forward by 1 second
t = 10;
shifted_signal = [zeros(1,t*sampling_freq), signal, zeros(1, nsamples_tot - nsamples_sig - t*sampling_freq)];
% plot(timeVecTot, shifted_signal);
%Add noise to shifted signal to create total signal 
total_signal = shifted_signal + wgn;
% plot(timeVecTot, total_signal);
% hold on;
% plot(timeVecTot, shifted_signal, 'r');
%Run Matched Filtering
[ta, timesVec] = matchedfiltering(total_signal, final_signal, sampling_freq, psd);
plot(timeVecTot, timesVec);
title("Time Sequence of Matched Filtering")
xlabel("Seconds");
ylabel("Amplitude");





