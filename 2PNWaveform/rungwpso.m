%% Minimize the fitness function PSOFITFUNC using PSO
%Specify initial parameters
T_sig = 54;
initial_phase = 0;
%% Sampling Frequency
Fs = 2048;
%% Number of samples = num*Fs*T_sig
num = 3;
N = floor(num*T_sig*Fs);
timeVec = (0:N-1)*(1/Fs);
%% Min and Max Frequencies Hz
fmin = 30;
fmax = 700;
%% Positive Frequency Vector
datalen = N/Fs;
fpos = (0:floor(N/2))*(1/datalen);
%% Initial Time of Arrival and Phase 
ta = 0;
phase = 0;
% Signal to noise ratio of the true signal
snr = 10;
% Phase coefficients parameters of the true signal
m1 = 1.4;
m2 = 1.4;
coeffs = [m1,m2];
% Search range of phase coefficients
rmin = [1, 1];
rmax = [10, 10];
% Number of independent PSO runs
nRuns = 2;
%% Do not change below
% Generate data realization
% dataX = (0:(nSamples-1))/Fs;
% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
% rng('default');
% N = 10;
% a1_errors = zeros(1,N);
% a2_errors = zeros(1,N);
% a3_errors = zeros(1,N);
% A_errors = zeros(1,N);
% phi_errors = zeros(1,N);
% ta_errors = zeros(1,N);
% Generate 2PN signal
wave = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N);
%Generate Final Signal
wgn = randn(1, N);
dataY = wave + wgn;
sizedataY = size(dataY);
dataX = timeVec;
% Input parameters for CRCBQCHRPPSO
inParams = struct('dataX', dataX,...
                  'fpos', fpos,...
                  'dataY', dataY,...
                  'frange', [fmin,fmax],...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'snr', snr,...
                  'N', N,...
                  'rmin',rmin,...
                  'rmax',rmax);
% CRCBQCHRPPSO runs PSO on the PSOFITFUNC fitness function. As an
% illustration of usage, we change one of the PSO parameters from its
% default value.
outStruct = crcbqcpso(inParams,struct('maxSteps',2000),nRuns,Fs);
% a1_errors(i) = outStruct.bestQcCoefs(1) - a1;
% a2_errors(i) = outStruct.bestQcCoefs(2) - a2;
% a3_errors(i) = outStruct.bestQcCoefs(3) - a3;
% A_errors(i) = outStruct.bestAmp - snr;
% phi_errors(i) = outStruct.bestPhase - pi/4;
% ta_errors(i) = outStruct.bestTime - t;
% end
%%Plots
% figure;
% plot(a1_errors);
% title('Estimated Errors for a1');
% ylabel('a1_estimated - a1_actual');
% 
% figure;
% plot(a2_errors);
% title('Estimated Errors for a2');
% ylabel('a2_estimated - a2_actual');
% 
% figure;
% plot(a3_errors);
% title('Estimated Errors for a1');
% ylabel('a3_estimated - a3_actual');
% 
% figure;
% plot(A_errors);
% title('Estimated Errors for A');
% ylabel('A_estimated - A_actual');
% 
% figure;
% plot(phi_errors);
% title('Estimated Errors for phi');
% ylabel('phi_estimated - phi_actual');
% 
% figure;
% plot(ta_errors);
% title('Estimated Errors for t_a');
% ylabel('ta_estimated - ta_actual');
% 
% 
% %%Histograms
% figure;
% histogram(a1_errors, 50);
% title('Histogram for Errors of a1');
% 
% figure;
% histogram(a2_errors, 50);
% title('Histogram for Errors of a2');
% 
% 
% figure;
% histogram(a3_errors, 50);
% title('Histogram for Errors of a3');
% 
% figure;
% histogram(A_errors, 50);
% title('Histogram for Errors of A');
% 
% figure;
% histogram(phi_errors, 50);
% title('Histogram for Errors of phi');
% 
% figure;
% histogram(ta_errors, 50);
% title('Histogram for Errors of ta');
%%
% Plots
figure;
hold on;
plot(dataX,dataY,'.');
plot(dataX,wave,'r');
for lpruns = 1:nRuns
      plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
legend('Data','Signal',...
        ['Estimated signal: ',num2str(nRuns),' runs'],...
        'Estimated signal: Best run');
% legend('Data','Signal',...
%        'Estimated signal: Best run');
disp(['Estimated parameters: m1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; m2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; A = ',num2str(outStruct.bestAmp),...
                             '; phi = ',num2str(outStruct.bestPhase),...
                              '; t_a = ',num2str(outStruct.bestTime)]);

