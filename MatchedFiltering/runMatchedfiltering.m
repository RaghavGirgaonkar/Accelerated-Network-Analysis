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
signal = genqc(timeVecSig,A,coeffs,0);
final_signal = [signal, zeros(1,(nsamples_tot-nsamples_sig))];
%Shift the signal forward by t seconds
t = 10;
shifted_signal = [zeros(1,floor(t*sampling_freq)-1), signal, zeros(1, nsamples_tot - nsamples_sig - floor(t*sampling_freq)+1)];
% plot(timeVecTot, shifted_signal);
niter = 1000;
SNRs = zeros(1,niter);

for i = 1:niter
%Plot signal 
%plot(timeVecTot, final_signal);
%Generate White Gaussian Noise
wgn = randn(1, nsamples_tot);
% plot(timeVecTot, wgn);
% hold on;
% plot(timeVecTot, final_signal, 'r');
%Add noise to shifted signal to create total signal 
total_signal = shifted_signal + wgn;
% plot(timeVecTot, total_signal);
% hold on;
% plot(timeVecTot, shifted_signal, 'r');
%Run Matched Filtering
[ta, timesVec] = matchedfiltering(total_signal, final_signal, sampling_freq, psd);
%Find SNR
[max_val, max_sample] = max(real(timesVec));
stdv = std(timesVec(1:4*sampling_freq));
snr = max_val/stdv;
fprintf("SNR is %f\n", snr);
SNRs(i) = snr;
end
%Plot Histogram
histogram(SNRs,50);
%Plot
% plot(timeVecTot, timesVec);
% title("Time Sequence of Matched Filtering")
% xlabel("Seconds");
% ylabel("Amplitude");





