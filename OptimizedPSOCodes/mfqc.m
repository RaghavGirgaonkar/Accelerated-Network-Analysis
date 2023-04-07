function [mfVal, max_arg] = mfqc(x,params)
%Generate normalized 2PN Waveform
% phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
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
