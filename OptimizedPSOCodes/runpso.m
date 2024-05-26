function [] = runpso(segdata, whtnddata,TFtotal, filename, segnum,Seglen)
%Script to run Chirp-time or Mass Space PSO on a user specified datasegment and PSD.
%Refer to User Manual for usage

%Raghav Girgaonkar, July 2023

%% Read JSON Files 
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
%% Specify initial parameters
%T_sig_len = params.signal.T_sig_len;
T_sig = params.signal.T_sig;
initial_phase = 0;
%% Sampling Frequency
Fs = params.sampling_freq;
%% Number of samples = num*Fs*T_sig
%num = params.signal.num;
N = floor(Seglen*Fs);
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
%% Signal to noise ratio of the true signal
snr = params.snr;
%% Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
%% Mass parameters of the true signal
m1 = params.masses(1);
m2 = params.masses(2);
%% Tau coeffs as phase parameters
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
%% Search range of phase coefficients
if type
    rmin = [params.rmin_tau(1), params.rmin_tau(2)];
    rmax = [params.rmax_tau(1), params.rmax_tau(2)];
else
    rmin = [params.rmin(1), params.rmin(2)];
    rmax = [params.rmax(1), params.rmax(2)];
end
%% Number of independent PSO runs
nRuns = pso.nruns;

%% Pre-processing
[A,avec, phaseDiff] = preprocessing(fmin,fmax, fpos, datalen, N);

dataY = whtnddata;

AbysqrtPSD = A.*TFtotal;

%% Create General Normalization Factor
% Scalar factor of 1/N (check NoteonNormalizaton.pdf for more details)
dataLen = N;
innProd = (1/dataLen)*(AbysqrtPSD)*AbysqrtPSD';
genNormfacSqr = real(innProd);
genNormfac = 1/sqrt(genNormfacSqr);

%% Data Products 
fftdataY = fft(dataY);
fftdataY = fftdataY.*A;%Pre-multiply frequency magnitude vector A for optimization

%% Get FFT of data by total PSD
fftdataYbyPSD = fftdataY.*TFtotal;

dataX = timeVec;

%% Input Parameters:
inParams = struct('dataX', dataX,...
                  'fpos', fpos,...
                  'dataY', dataY,...
                  'fftdataYbyPSD', fftdataYbyPSD,...
                  'frange', [fmin,fmax],...
                  'datalen',datalen,...,
                  'initial_phase', initial_phase,...
                  'N', N,...
                  'A', A,...
                  'phaseDiff', phaseDiff,...
                  'normfac', genNormfac,...
                  'avec', avec,...
		          'T_sig', T_sig,...
                  'rmin',rmin,...
                  'rmax',rmax,...
                  'Fs',Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxSteps = pso.maxSteps;
if type
    if rmin(1) >= 0 && rmin(2) >= 0
    	original_fitVal = -1*mfgw_tau([tau0, tau1p5], inParams);
    	outStruct = crcbgwpso_tau(inParams,struct('maxSteps',maxSteps),nRuns,Fs);
    else
	original_fitVal = -1*mfgw_tau_negative([tau0, tau1p5], inParams);
    	outStruct = crcbgwpso_tau_negative(inParams,struct('maxSteps',maxSteps),nRuns,Fs);
    end
    bestFitVal = -1*outStruct.bestFitness;
else
    original_fitVal = -1*mfgw_mass([m1, m2], inParams);
    outStruct = crcbgwpso_mass(inParams,struct('maxSteps',maxSteps),nRuns,Fs);
    bestFitVal = -1*outStruct.bestFitness;
end

%% Optional: Uncomment to save output of PSO in a .mat file
% save(files.output_struct_location,'outStruct');

%% Plots, uncomment all saveas() commands to save figures in custom directories
figure;
hold on;
plot(dataX,segdata,'.');
% plot(dataX,wave,'r');
% for lpruns = 1:nRuns
%       plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
% end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
% legend('Data','Signal',...
%         ['Estimated signal: ',num2str(nRuns),' runs'],...
%         'Estimated signal: Best run');
legend('Data','Signal',...
        'Estimated signal: Best run');
% saveas(gcf,files.psoresultplot);
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
% saveas(gcf,files.bestfitplot);
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
    boundary_plot;
%     saveas(gcf,files.bestlocplot);
    hold off;

    

    t0 = outStruct.bestQcCoefs(1);
    t1p5 = outStruct.bestQcCoefs(2);
    est_ta = outStruct.bestTime;
    est_SNR = outStruct.bestAmp;
    est_M = (5/(32*fmin))*(t1p5/(pi*pi*t0))*(c^3/G);
    est_u = (1/(16*fmin*fmin))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(c^3/G);
    
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/(2*Msolar);
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/(2*Msolar);
    
    %% This will display parameters given through signal.json and PSO-estimated parameters
    %% Uncomment Original parameter display command if needed
%     disp(['Original parameters: tau0= ',num2str(tau0),...
%                                   '; tau1p5= ',num2str(tau1p5),...
%                                   '; m1= ', num2str(m1/Msolar),...
%                                   '; m2= ', num2str(m2/Msolar),...
%                                   '; A = ',num2str(snr),...
%                                  '; phi = ',num2str(phase),...
%                                   '; t_a = ',num2str(ta),...
%                                   '; FitVal = ',num2str(original_fitVal)]);
    
    disp(['Estimated parameters: tau0=',num2str(outStruct.bestQcCoefs(1)),...
                                  '; tau1p5=',num2str(outStruct.bestQcCoefs(2)),...
                                  '; m1= ', num2str(est_m1),...
                                  '; m2= ', num2str(est_m2),...
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
%     saveas(gcf,files.bestlocplot);
    hold off;
    
    %% This will display parameters given through signal.json and PSO-estimated parameters
    %% Uncomment Original parameter display command if needed
%     disp(['Original parameters:  m1= ',num2str(m1),...
%                                   '; m2= ',num2str(m2),...
%                                   '; A = ',num2str(snr),...
%                                  '; phi = ',num2str(phase),...
%                                   '; t_a = ',num2str(ta),...
%                                   '; FitVal = ',num2str(original_fitVal)]);
	
    est_m1 = outStruct.bestQcCoefs(1);
    est_m2 = outStruct.bestQcCoefs(2);
    est_ta = outStruct.bestTime;
    est_SNR = outStruct.bestAmp;
    disp(['Estimated parameters: m1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; m2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; A = ',num2str(outStruct.bestAmp),...
                             '; phi = ',num2str(outStruct.bestPhase),...
                              '; t_a = ',num2str(outStruct.bestTime),...
                              '; FitVal = ',num2str(bestFitVal)]);
end

%{
S = load('injectedsignals.mat'); sigsegnum = S.segnum; sigsnrs = S.injectedsnrs;

index = find(sigsegnum == segnum);

if index
    signalsnr = sigsnrs(index);
else
    signalsnr = 0;
end
%}


%candidates = [segnum, t0, t1p5, real(est_m1), imag(est_m1), real(est_m2), imag(est_m2), est_ta, est_SNR,signalsnr];
candidates = [segnum, t0, t1p5, real(est_m1), imag(est_m1), real(est_m2), imag(est_m2), est_ta, est_SNR];
%candidates = [segnum, est_m1, est_m2, est_ta, est_SNR];
fileID = fopen(files.candidatefile,'a');
fprintf(fileID,'%d\t%f\t%f\t%f + %fi\t\t%f + %fi\t\t%f\t\t%f\n',candidates);
%fprintf(fileID,'%d\t%f\t%f\t%f + %fi\t\t%f + %fi\t\t%f\t\t%f\t%f\n',candidates);
%fprintf(fileID,'%d\t%f\t%f\t%f\t%f\n',candidates);
fclose(fileID);
end


