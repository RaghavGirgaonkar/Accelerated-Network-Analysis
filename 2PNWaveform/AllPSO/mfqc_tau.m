function [mfVal, max_arg] = mfqc_tau(x,params)
%Generate normalized 2PN Waveform
% phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
tau0 = x(1);
tau1p5 = x(2);
q0 = gen2PNwaveform_tau(params.fpos, 0, 0, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.psd);
q1 = gen2PNwaveform_tau(params.fpos, 0, pi/2, params.frange(1), params.frange(2), tau0,...
    tau1p5,params.datalen,0,1,params.N,params.psd);


%Compute fitness value after maximizing by matched filtering
mf1 = matchedfiltering(params.dataY, q0, params.Fs,params.psd);
mf2 = matchedfiltering(params.dataY, q1, params.Fs,params.psd);
[max_val, max_arg] = max(mf1.^2 + mf2.^2);
mfVal = -1*max_val;