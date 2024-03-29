function outResults = crcbgwpso_tau_negative(inParams,psoParams,nRuns, sampling_freq)
%Regression of 2PNWaveform using Chirp-Time Space PSO
% Input: inParams: Struct containing data and signal parameters
%        psoParams: Struct containing PSO parameters
%        nRuns: Number of PSO iterations
%        sampling_freq: Sampling Frequency of data
% Output: outResults: Struct containing the following parameters
% 'allRunsOutput': An N element struct array containing results from each PSO
%              run. The fields of this struct are:
%                 'fitVal': The fitness value.
%                 'qcCoefs': The coefficients [tau0, tau1.5].
%                 'estSig': The estimated signal.
%                 'totalFuncEvals': The total number of fitness
%                                   evaluations.
% 'bestRun': The best run.
% 'bestFitness': best fitness from the best run.
% 'bestSig' : The signal estimated in the best run.
% 'bestQcCoefs' : [tau0, tau1.5] found in the best run.

%Raghav Girgaonkar, April 2023

nSamples = length(inParams.dataX);

fHandle = @(x) psofitfunc_tau_negative(x,inParams);

params = inParams;

nDim = 2;
outStruct = struct('bestLocation',[],...
                   'bestFitness', [],...
                   'totalFuncEvals',[],...
                   'allBestFit',[],...
                   'allBestLoc',[]);
                    
outResults = struct('allRunsOutput',struct('fitVal', [],...
                                           'qcCoefs',zeros(1,2),...
                                           'estTa',[],...
                                           'estSig',zeros(1,nSamples),...
                                           'totalFuncEvals',[],...
                                           'allBestFit',zeros(1,psoParams.maxSteps),...
                                           'allBestLoc',zeros(nDim,psoParams.maxSteps)),...
                    'bestRun',[],...
                    'bestFitness',[],...
                    'bestSig', zeros(1,nSamples),...
                    'bestQcCoefs',zeros(1,2),...
                    'estAmp',[],...
                    'estPhase',[],...
                    'bestAmp',[],...
                    'bestPhase',[],...
                    'bestTime',[]);

%Allocate storage for outputs: results from all runs are stored
for lpruns = 1:nRuns
    outStruct(lpruns) = outStruct(1);
    outResults.allRunsOutput(lpruns) = outResults.allRunsOutput(1);
end
%Independent runs of PSO in parallel. Change 'parfor' to 'for' if the
%parallel computing toolbox is not available.
% fprintf("Running PSO\n");
parpool(nRuns);
parfor lpruns = 1:nRuns
    %Reset random number generator for each worker
    rng(lpruns);
    outStruct(lpruns)=crcbpso(fHandle,nDim,psoParams,2);
end

%Prepare output
fitVal = zeros(1,nRuns);
% sampling_interval = 1/sampling_freq;
%Time vectors for signal and total series
% timeVecSig = (0:(sampling_freq*T_sig - 1))/sampling_freq;
for lpruns = 1:nRuns 
    outResults.allRunsOutput(lpruns).allBestFit = outStruct(lpruns).allBestFit;
    outResults.allRunsOutput(lpruns).allBestLoc = outStruct(lpruns).allBestLoc;
    fitVal(lpruns) = outStruct(lpruns).bestFitness;
    outResults.allRunsOutput(lpruns).fitVal = fitVal(lpruns);
    [~,chirptimes,~] = fHandle(outStruct(lpruns).bestLocation);
%     Q = qcCoefs
%     index = ta_index
    [estAmp, estTa, estPhase] = getparamestimates_negative(chirptimes, params);
    outResults.allRunsOutput(lpruns).qcCoefs = chirptimes;
    %Calculate time using sampling freq and ta_index
%     estTa = ta_index/sampling_freq;
    
    outResults.allRunsOutput(lpruns).estTa = estTa;
    tau0 = chirptimes(1);
    tau1p5 = chirptimes(2);
%     phaseq0 = gen2PNwaveform_tau_negative(params.fpos, estTa, 0, params.frange(1), params.frange(2), tau0,...
%     tau1p5,params.datalen,0,1,params.N,params.avec, params.normfac);
%     fftq0 = phaseq0;
%     fftq1 = phaseq0.*params.phaseDiff;
%     %Estimated Phase
% %     yq0 = inParams.dataY*q0(:);
% %     yq1 = inParams.dataY*q1(:);
%     yq0 = innerprodpsd(fftq0, params.fftdataYbyPSD);
%     yq1 = innerprodpsd(fftq1, params.fftdataYbyPSD);
%     estPhase = atan2(yq1,yq0);
    outResults.allRunsOutput(lpruns).estPhase = estPhase;

    %Estimated Amplitude
%     estAmp = cos(estPhase)*yq0 + sin(estPhase)*yq1;
    outResults.allRunsOutput(lpruns).estAmp = estAmp;
    
    %Estimated Signal
%     estSigTemp = genqc(timeVecSig,1,qcCoefs,estPhase);
    estSigphase = gen2PNwaveform_tau_negative(params.fpos, estTa, estPhase, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,estAmp,params.N,params.avec, params.normfac);
    estSigfourier = (params.A).*estSigphase;
    estSig = ifft(estSigfourier);
%     estSigTemp_shifted = [zeros(1,floor(estTa*sampling_freq)-1), estSigTemp, zeros(1, nSamples - nSamplesSig - floor(estTa*sampling_freq)+1)];
%     estSig = estAmp*estSigTemp;
    outResults.allRunsOutput(lpruns).estSig = estSig;
    outResults.allRunsOutput(lpruns).totalFuncEvals = outStruct(lpruns).totalFuncEvals;
end
%Find the best run
[~,bestRun] = min(fitVal(:));
outResults.bestRun = bestRun;
outResults.bestFitness = outResults.allRunsOutput(bestRun).fitVal;
outResults.bestSig = outResults.allRunsOutput(bestRun).estSig;
outResults.bestAmp = outResults.allRunsOutput(bestRun).estAmp;
outResults.bestPhase = outResults.allRunsOutput(bestRun).estPhase;
outResults.bestQcCoefs = outResults.allRunsOutput(bestRun).qcCoefs;
outResults.bestTime = outResults.allRunsOutput(bestRun).estTa;

