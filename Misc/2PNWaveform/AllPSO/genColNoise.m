function outNoise = genColNoise(PSDs, Fs,frange, noise_num)
% PSDs is assumed to be a row vector with PSDs for postive DFT frequencies
[~,y] = size(PSDs);
N = y+y-2;
datalen = N/Fs;%Total datalength in seconds

fvec = 0:(1/datalen):(Fs/2);

fmin = frange(1);
fmax = frange(2);

% Modifications
minidx = find(fvec==fmin);
maxidx = find(fvec==fmax);

Snmin = PSDs(minidx);
Snmax = PSDs(maxidx);

PSDs(1:minidx) = Snmin;
PSDs(maxidx:end) = Snmax;

% loglog(totfvec,TotalpsdVec);
 
%% Make colored Noise
fltrOrdr = 500;
% 
outNoise_t = statgaussnoisegen(N+fltrOrdr,[fvec(:),PSDs(:)],fltrOrdr,Fs, noise_num);

outNoise = outNoise_t(fltrOrdr+1:end);
% plot(timeVec,outNoise);




