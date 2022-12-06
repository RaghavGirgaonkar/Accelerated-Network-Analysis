function []=rungwpso(filename)
%% Read JSON File 
fname = filename;
str = fileread(fname);
params = jsondecode(str);
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
if params.pso.type == "tau"
    m1 = m1*Msolar;
    m2 = m2*Msolar;
    M = m1 + m2;
    u = m1*m2/(m1 + m2);
    n = u/M;
    tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
    tau1p5 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1))*(1/n)*((743/336)+ (11*n/4));
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
nRuns = params.pso.nruns;
%% Do not change below
% Generate data realization
% dataX = (0:(nSamples-1))/Fs;
% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
% rng('default');

% Generate 2PN signal
if type
    wave = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0,tau1p5,datalen, initial_phase, snr, N);
else
    wave = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N);
end
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
maxSteps = params.pso.maxSteps;
if type
    original_fitVal = -1*mfqc_tau([tau0, tau1p5], inParams);
    outStruct = crcbqcpso_tau(inParams,struct('maxSteps',maxSteps),nRuns,Fs);
    bestFitVal = -1*outStruct.bestFitness;
else
    original_fitVal = -1*mfqc([m1, m2], inParams);
    outStruct = crcbqcpso(inParams,struct('maxSteps',maxSteps),nRuns,Fs);
    bestFitVal = -1*outStruct.bestFitness;
end
% save('/scratch/09197/raghav/outStruct_tau.mat','outStruct');

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
saveas(gcf,params.files.psoresultplot);
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
saveas(gcf,params.files.bestfitplot);
hold off;

if type
    figure;
    hold on;
    for lpruns = 1:nRuns
          rVec = s2rv(outStruct.allRunsOutput(lpruns).allBestLoc,inParams);
          plot(rVec(:,1),rVec(:,2),'DisplayName',num2str(lpruns));
    end
    scatter(tau0,tau1p5,140,'red','filled','D','DisplayName','Original Parameters');
    title("Best Parameter Values for All Runs");
    xlabel("\tau_0");
    ylabel("\tau_{1.5}");
    legend;
    boundary_fig = boundary_plot;
    saveas(gcf,params.files.bestlocplot);
    hold off;

    

    t0 = outStruct.bestQcCoefs(1);
    t1p5 = outStruct.bestQcCoefs(2);
    est_M = (5/(32*fmin))*(t1p5/(pi*pi*t0))*(c^3/G);
    est_u = (1/(16*fmin*fmin))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(c^3/G);
    
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;
    
    
    disp(['Original parameters: tau0= ',num2str(tau0),...
                                  '; tau1p5= ',num2str(tau1p5),...
                                  '; m1= ', num2str(m1/Msolar),...
                                  '; m2= ', num2str(m2/Msolar),...
                                  '; A = ',num2str(snr),...
                                 '; phi = ',num2str(phase),...
                                  '; t_a = ',num2str(ta),...
                                  '; FitVal = ',num2str(original_fitVal)]);
    
    disp(['Estimated parameters: tau0=',num2str(outStruct.bestQcCoefs(1)),...
                                  '; tau1p5=',num2str(outStruct.bestQcCoefs(2)),...
                                  '; m1= ', num2str(est_m1/Msolar),...
                                  '; m2= ', num2str(est_m2/Msolar),...
                                  '; A = ',num2str(outStruct.bestAmp),...
                                 '; phi = ',num2str(outStruct.bestPhase),...
                                  '; t_a = ',num2str(outStruct.bestTime),...
                                  '; FitVal = ',num2str(bestFitVal)]);
else
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
    saveas(gcf,params.files.bestlocplot);
    hold off;

    disp(['Original parameters:  m1= ',num2str(m1),...
                                  '; m2= ',num2str(m2),...
                                  '; A = ',num2str(snr),...
                                 '; phi = ',num2str(phase),...
                                  '; t_a = ',num2str(ta),...
                                  '; FitVal = ',num2str(original_fitVal)]);

    disp(['Estimated parameters: m1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; m2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; A = ',num2str(outStruct.bestAmp),...
                             '; phi = ',num2str(outStruct.bestPhase),...
                              '; t_a = ',num2str(outStruct.bestTime),...
                              '; FitVal = ',num2str(bestFitVal)]);
end

fitnessVals = [original_fitVal, bestFitVal];
fileID = fopen(params.files.fitnessvaluesfile,'a');
fprintf(fileID,'%f %f\n',fitnessVals);
fclose(fileID);

end