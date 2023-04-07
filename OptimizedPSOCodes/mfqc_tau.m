function [mfVal, max_arg] = mfqc_tau(x,params)
%MatchedFiltering for Chirp time space PSO
%Raghav Girgaonkar, April 2023

%Generate normalized quadrature templates
tau0 = x(1);
tau1p5 = x(2);
phaseq0 = gen2PNwaveform_tau(params.fpos, 0, 0, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.avec, params.normfac);

fftq0 = phaseq0;

fftq1 = phaseq0.*params.phaseDiff;

%Compute fitness value after maximizing by matched filtering
mf1 = matchedfiltering(params.fftdataYbyPSD, fftq0);
mf2 = matchedfiltering(params.fftdataYbyPSD, fftq1);
[max_val, max_arg] = max(mf1(1:end - params.T_sig*params.Fs).^2 + mf2(1:end - params.T_sig*params.Fs).^2);
mfVal = -1*max_val;
