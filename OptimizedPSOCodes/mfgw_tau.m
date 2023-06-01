function [mfVal, max_arg] = mfgw_tau(x,params)
%MatchedFiltering for Chirp time space PSO
%Generates a combined matched-filtering timeseries from both quadrature templates and 
%returns index and value of the maximum of this series,
%Input: x = [tau0, tau1.5], vector containing chirp-time parameters for
%           creating quadrature templates
%       params: Struct containing signal parameters
%Output: mfVal: Maximum value of total matchedfiltering timeseries
%        max_arg: Index of maximum value

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
