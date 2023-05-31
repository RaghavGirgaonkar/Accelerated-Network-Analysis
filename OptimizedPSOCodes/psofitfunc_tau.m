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
        [fitVal(lpc), max_index] = mfgw_tau(x, params);
    end
end

%Return max_index if requested
if nargout > 1
    varargout{1} = xVec;
    varargout{2}=max_index;
end

