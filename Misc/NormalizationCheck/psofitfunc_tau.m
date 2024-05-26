function [fitVal,varargout] = psofitfunc_tau(xVec,params)
%Fitness function for Chirp-time Space PSO
%Input: xVec: normalized location vector of mass parameters
%       params: Struct containing signal and data parameters
%Output: fitval: Fitness value at location specified by xVec
%        varargout: Additional output arguments sucha as the index of max
%        value in matchedfiltering timeseries

%Raghav Girgaonkar, April 2023

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
        [fitVal(lpc), max_index] = mfgw_tau_whtnd(x, params);
    end
end

%Return max_index if requested
if nargout > 1
    varargout{1} = xVec;
    varargout{2}=max_index;
end

