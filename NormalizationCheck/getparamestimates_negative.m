function [estAmp, estTa, estPhase] = getparamestimates_negative(chirptimes, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tau0 = chirptimes(1);
tau1p5 = chirptimes(2);
phaseq0 = gen2PNwaveform_tau_negative(params.fpos, 0, 0, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.avec, params.normfac);

fftq0 = phaseq0;

fftq1 = phaseq0.*params.phaseDiff;

%Compute fitness value after maximizing by matched filtering
mf1 = matchedfiltering(params.fftdataYbyPSD, fftq0);
mf2 = matchedfiltering(params.fftdataYbyPSD, fftq1);
[max_val, max_arg] = max(mf1(1:end - params.T_sig*params.Fs).^2 + mf2(1:end - params.T_sig*params.Fs).^2);
%Estimated SNR
estAmp = sqrt(max_val);
%Estimated TOA:
estTa = max_arg/params.Fs;

yq0 = mf1(max_arg);
yq1 = mf2(max_arg);

estPhase = atan2(yq1,yq0);
end

