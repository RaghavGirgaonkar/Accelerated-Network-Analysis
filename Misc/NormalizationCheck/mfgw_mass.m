function [mfVal, max_arg] = mfgw_mass(x,params)
%MatchedFiltering for Mass space PSO
%Generates a combined matched-filtering timeseries from both quadrature templates and 
%returns index and value of the maximum of this series,
%Input: x = [m1, m2], vector containing mass parameters for
%           creating quadrature templates
%       params: Struct containing signal parameters
%Output: mfVal: Maximum value of total matchedfiltering timeseries
%        max_arg: Index of maximum value


%Raghav Girgaonkar, April 2023
m1 = x(1);
m2 = x(2);
phaseq0 = gen2PNwaveform(params.fpos, 0, 0, params.frange(1), params.frange(2), m1,...
    m2,params.datalen,0,1,params.N,params.avec, params.normfac);

fftq0 = phaseq0;

fftq1 = phaseq0.*params.phaseDiff;


%Compute fitness value after maximizing by matched filtering
mf1 = matchedfiltering(params.fftdataYbyPSD, fftq0);
mf2 = matchedfiltering(params.fftdataYbyPSD, fftq1);
[max_val, max_arg] = max(mf1(1:end - params.T_sig*params.Fs).^2 + mf2(1:end - params.T_sig*params.Fs).^2);
mfVal = -1*max_val;
