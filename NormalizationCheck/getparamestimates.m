function [estAmp, estTa, estPhase] = getparamestimates(x, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Generate normalized quadrature templates
tau0 = x(1);
tau1p5 = x(2);
q0 = gen2PNtemplate_tau(params.fpos, 0, 0, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.A,params.avec, params.PSDtotal);
q1 = gen2PNtemplate_tau(params.fpos, 0, pi/2, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.A,params.avec, params.PSDtotal);

%Compute fitness value after maximizing by matched filtering
% mf1 = mlftr(params.dataY, q0, params.PSDtotal);
% mf2 = mlftr(params.dataY, q1, params.PSDtotal);
mf1 = mlftr2(params.dataY, q0, params.TFtotal);
mf2 = mlftr2(params.dataY, q1, params.TFtotal);
[max_val, max_arg] = max(mf1(1:end - params.T_sig*params.Fs).^2 + mf2(1:end - params.T_sig*params.Fs).^2);
mfVal = -1*max_val;
%Estimated SNR
estAmp = sqrt(mfVal);
%Estimated TOA:
estTa = max_arg/params.Fs;

yq0 = mf1(max_arg);
yq1 = mf2(max_arg);

estPhase = atan2(yq1,yq0);
end

