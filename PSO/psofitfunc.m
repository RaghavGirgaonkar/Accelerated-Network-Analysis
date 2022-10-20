function [fitVal,varargout] = psofitfunc(xVec,params)


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
if nargout > 2
    varargout{1} = xVec;
    varargout{2}=max_index;
end

%Max of matchedfiltering series after maximizing over amplitude and phase
%parameters
function [mfVal, max_arg] = mfqc(x,params)
%Generate normalized quadratic chirp
phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
q0 = sin(2*pi*phaseVec);
q0 = q0/norm(q0);
q1 = sin(2*pi*phaseVec + pi/2);
q1 = q1/norm(q1);

%Compute fitness value after maximizing by matched filtering
m1 = matchedfiltering(params.dataY, q0);
m2 = matchedfiltering(params.dataY, q1);
[max_val, max_arg] = max(m1.^2 + m2.^2);
mfVal = -1*max_val;