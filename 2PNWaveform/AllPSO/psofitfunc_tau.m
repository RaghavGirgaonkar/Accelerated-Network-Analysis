function [fitVal,varargout] = psofitfunc_tau(xVec,params)


%rows: points
%columns: coordinates of a point
[nVecs,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nVecs,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
xVec(validPts,:) = s2rv(xVec(validPts,:),params);

for lpc = 1:nVecs
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        [fitVal(lpc), max_index] = mfqc(x, params);
    end
end

%Return max_index if requested
if nargout > 1
    varargout{1} = xVec;
    varargout{2}=max_index;
end

%Max of matchedfiltering series after maximizing over amplitude and phase
%parameters
function [mfVal, max_arg] = mfqc(x,params)
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