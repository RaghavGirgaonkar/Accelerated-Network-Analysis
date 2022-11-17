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
%% Initial Time of Arrival and Coalescence Phase 
ta = 82;
phase = 0;
% Signal to noise ratio of the true signal
snr = 10;
% Phase coefficients parameters of the true signal
m1 = 1.4;
m2 = 1.4;
% coeffs = [m1,m2];
% Search range of phase coefficients
rmin = [1.2, 1.2];
rmax = [10, 10];
% Number of independent PSO runs
nRuns = 8;
%% Do not change below
% Generate data realization
% dataX = (0:(nSamples-1))/Fs;
% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
% rng('default');

% Generate 2PN signal
wave = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N);
%Generate Final Signal
wgn = randn(1, N);
dataY = wave + wgn;
dataX = timeVec;
% Input parameters for CRCBQCHRPPSO
inParams = struct('dataX', dataX,...
                  'fpos', fpos,...
                  'dataY', dataY,...
                  'frange', [fmin,fmax],...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'rmin',rmin,...
                  'rmax',rmax);
% CRCBQCHRPPSO runs PSO on the PSOFITFUNC fitness function. As an
% illustration of usage, we change one of the PSO parameters from its
% default value.
maxSteps = 10;
outStruct = crcbqcpso(inParams,struct('maxSteps',maxSteps),nRuns,Fs);

save('/scratch/09197/raghav/outStruct.mat','outStruct');

% Plots
figure;
hold on;
plot(dataX,dataY,'.');
plot(dataX,wave,'r');
% for lpruns = 1:nRuns
%       plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
% end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
% legend('Data','Signal',...
%         ['Estimated signal: ',num2str(nRuns),' runs'],...
%         'Estimated signal: Best run');
legend('Data','Signal',...
        'Estimated signal: Best run');
saveas(gcf,"/scratch/09197/raghav/psoresults.pdf");
hold off;

figure;
iterVec = linspace(1,maxSteps,maxSteps);
hold on;
for lpruns = 1:nRuns
      plot(iterVec,outStruct.allRunsOutput(lpruns).allBestFit, 'DisplayName',num2str(lpruns));
end
title("Best Fitness Values for All Runs");
xlabel("Iteration");
ylabel("Best Fitness Value");
legend;
saveas(gcf,"/scratch/09197/raghav/bestfitness.pdf");
hold off;

figure;
hold on;
for lpruns = 1:nRuns
      rVec = s2rv(outStruct.allRunsOutput(lpruns).allBestLoc,inParams);
      plot(rVec(:,1),rVec(:,2),'DisplayName',num2str(lpruns));
end
scatter(m1,m2,140,'red','filled','D','DisplayName','Original Parameters');
title("Best Parameter Values for All Runs");
xlabel("m_1");
ylabel("m_2");
legend;
saveas(gcf,"/scratch/09197/raghav/bestloc.pdf");
hold off;



disp(['Estimated parameters: m1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; m2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; A = ',num2str(outStruct.bestAmp),...
                             '; phi = ',num2str(outStruct.bestPhase),...
                              '; t_a = ',num2str(outStruct.bestTime)]);

