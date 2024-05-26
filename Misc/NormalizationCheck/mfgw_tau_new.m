function [mf1, mf2, mfVal, max_arg] = mfgw_tau_new(x,params)
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
q0 = gen2PNtemplate_tau(params.fpos, 0, 0, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.A,params.avec, params.PSDtotal);
q1 = gen2PNtemplate_tau(params.fpos, 0, pi/2, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.A,params.avec, params.PSDtotal);

%Compute fitness value after maximizing by matched filtering
mf1 = mlftr(params.dataY, q0, params.PSDtotal);
mf2 = mlftr(params.dataY, q1, params.PSDtotal);
% mf1 = mlftr2(params.dataY, q0, params.TFtotal);
% mf2 = mlftr2(params.dataY, q1, params.TFtotal);
% [max_val, max_arg] = max(mf1(1:end - params.T_sig*params.Fs).^2 + mf2(1:end - params.T_sig*params.Fs).^2);
[max_val, max_arg] = max(mf1(1:end - params.T_sig*params.Fs).^2 + mf2(1:end - params.T_sig*params.Fs).^2);
mfVal = -1*max_val;
