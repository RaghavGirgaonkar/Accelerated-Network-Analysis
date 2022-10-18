function outResults = crcbqcpso(inParams,psoParams,nRuns, timeshift, T_sig,sampling_freq)
%Regression of quadratic chirp using PSO
%O = CRCBQCPPSO(I,P,N)
%I is the input struct with the fields given below.  P is the PSO parameter
%struct. Setting P to [] will invoke default parameters (see CRCBPSO). N is
%the number of independent PSO runs. The output is returned in the struct
%O. The fields of I are:
% 'dataY': The data vector (a uniformly sampled time series).
% 'dataX': The time stamps of the data samples.
% 'dataXSq': dataX.^2
% 'dataXCb': dataX.^3
% 'rmin', 'rmax': The minimum and maximum values of the three parameters
%                 a1, a2, a3 in the candidate signal:
%                 a1*dataX+a2*dataXSq+a3*dataXCb
%The fields of O are:
% 'allRunsOutput': An N element struct array containing results from each PSO
%              run. The fields of this struct are:
%                 'fitVal': The fitness value.
%                 'qcCoefs': The coefficients [a1, a2, a3].
%                 'estSig': The estimated signal.
%                 'totalFuncEvals': The total number of fitness
%                                   evaluations.
% 'bestRun': The best run.
% 'bestFitness': best fitness from the best run.
% 'bestSig' : The signal estimated in the best run.
% 'bestQcCoefs' : [a1, a2, a3] found in the best run.

%Soumya D. Mohanty, May 2018

nSamples = length(inParams.dataX);

nSamplesSig = sampling_freq*T_sig;

fHandle = @(x) psofitfunc(x,inParams);

nDim = 3;
outStruct = struct('bestLocation',[],...
                   'bestFitness', [],...
                   'totalFuncEvals',[]);
                    
outResults = struct('allRunsOutput',struct('fitVal', [],...
                                           'qcCoefs',zeros(1,3),...
                                           'estSig',zeros(1,nSamples),...
                                           'totalFuncEvals',[]),...
                    'bestRun',[],...
                    'bestFitness',[],...
                    'bestSig', zeros(1,nSamples),...
                    'bestQcCoefs',zeros(1,3),...
                    'estAmp',[],...
                    'estPhase',[],...
                    'bestAmp',[],...
                    'bestPhase',[]);

%Allocate storage for outputs: results from all runs are stored
for lpruns = 1:nRuns
    outStruct(lpruns) = outStruct(1);
    outResults.allRunsOutput(lpruns) = outResults.allRunsOutput(1);
end
%Independent runs of PSO in parallel. Change 'parfor' to 'for' if the
%parallel computing toolbox is not available.
parfor lpruns = 1:nRuns
    %Reset random number generator for each worker
    rng(lpruns);
    outStruct(lpruns)=crcbpso(fHandle,nDim,psoParams);
end

%Prepare output
fitVal = zeros(1,nRuns);
% sampling_interval = 1/sampling_freq;
%Time vectors for signal and total series
timeVecSig = (0:(sampling_freq*T_sig - 1))/sampling_freq;
for lpruns = 1:nRuns   
    fitVal(lpruns) = outStruct(lpruns).bestFitness;
    outResults.allRunsOutput(lpruns).fitVal = fitVal(lpruns);
    [~,qcCoefs] = fHandle(outStruct(lpruns).bestLocation);
    outResults.allRunsOutput(lpruns).qcCoefs = qcCoefs;
    estSigq0 = genqc(timeVecSig,1,qcCoefs,0);
    estSigq1 = genqc(timeVecSig,1,qcCoefs,pi/2);
    estSigq0_shifted = [zeros(1,floor(timeshift*sampling_freq)-1), estSigq0, zeros(1, nSamples - nSamplesSig - floor(timeshift*sampling_freq)+1)];
    estSigq1_shifted = [zeros(1,floor(timeshift*sampling_freq)-1), estSigq1, zeros(1, nSamples - nSamplesSig - floor(timeshift*sampling_freq)+1)];
%     sizeq0 = size(estSigq0_shifted)
    %Estimated Phase
    yq0 = inParams.dataY*estSigq0_shifted(:);
    yq1 = inParams.dataY*estSigq1_shifted(:);
    estPhase = atan(yq1/yq0);
    outResults.allRunsOutput(lpruns).estPhase = estPhase;
    estSigTemp = genqc(timeVecSig,1,qcCoefs,estPhase);
    estSigTemp_shifted = [zeros(1,floor(timeshift*sampling_freq)-1), estSigTemp, zeros(1, nSamples - nSamplesSig - floor(timeshift*sampling_freq)+1)];
    %Estimated Amplitude
    estAmp = inParams.dataY*estSigTemp_shifted(:);
    outResults.allRunsOutput(lpruns).estAmp = estAmp;
    estSig = estAmp*estSigTemp_shifted;
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

