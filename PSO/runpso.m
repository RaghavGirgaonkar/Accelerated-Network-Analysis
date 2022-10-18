%% Minimize the fitness function CRCBQCFITFUNC using PSO
%Specify initial parameters
%Time length of signal and total time (in seconds)
T_sig = 2;
T_total = 10;
% Signal to noise ratio of the true signal
snr = 10;
% Phase coefficients parameters of the true signal
a1 = 10;
a2 = 3;
a3 = 3;
coeffs = [a1,a2,a3];
% Search range of phase coefficients
rmin = [1, 1, 1];
rmax = [15, 5, 5];
%Sampling frequency (assumed 4 times that of the nyquist freq in specified QC)
nyquist_freq = 2*(coeffs(1) + 2*coeffs(2)*T_sig + 3*coeffs(3)*T_sig.^2);
sampling_freq = 4*nyquist_freq;
sampling_interval = 1/sampling_freq;
%Time vectors for signal and total series
timeVecSig = (0:(sampling_freq*T_sig - 1))/sampling_freq;
timeVecTot = (0:(sampling_freq*T_total - 1))/sampling_freq;
%Number of samples
nsamples_sig = length(timeVecSig);
nsamples_tot = length(timeVecTot);
% Number of independent PSO runs
nRuns = 8;
%% Do not change below
% Generate data realization
% dataX = (0:(nSamples-1))/Fs;
% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
rng('default');
% Generate pure signal
sig = genqc(timeVecSig,snr,[a1,a2,a3],pi/4);
%Shift the signal forward by t seconds
t = 6;
shifted_signal = [zeros(1,floor(t*sampling_freq)-1), sig, zeros(1, nsamples_tot - nsamples_sig - floor(t*sampling_freq)+1)];
%Generate Final Signal
wgn = randn(1, nsamples_tot);
dataY = shifted_signal + wgn;
sizedataY = size(dataY)
dataX = timeVecTot;
% Input parameters for CRCBQCHRPPSO
inParams = struct('dataX', dataX,...
                  'dataY', dataY,...
                  'dataXSq',dataX.^2,...
                  'dataXCb',dataX.^3,...
                  'rmin',rmin,...
                  'rmax',rmax);
% CRCBQCHRPPSO runs PSO on the CRCBQCHRPFITFUNC fitness function. As an
% illustration of usage, we change one of the PSO parameters from its
% default value.
outStruct = crcbqcpso(inParams,struct('maxSteps',2000),nRuns, t, T_sig, sampling_freq);

%%
% Plots
figure;
hold on;
plot(dataX,dataY,'.');
plot(dataX,shifted_signal,'r');
% for lpruns = 1:nRuns
%     plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
% end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
legend('Data','Signal',...
       ['Estimated signal: ',num2str(nRuns),' runs'],...
       'Estimated signal: Best run');
disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                             '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                             '; a3=',num2str(outStruct.bestQcCoefs(3)),...
                             '; A = ',num2str(outStruct.bestAmp),...
                             '; \phi = ',num2str(outStruct.bestPhase),...
                             '; t_a = ',num2str(t)]);

