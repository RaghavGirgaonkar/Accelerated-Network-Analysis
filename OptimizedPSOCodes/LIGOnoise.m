function [outNoise, interPSD] = LIGOnoise(T_sig, num, Fs,noise_num)
%Load PSD 
y = load('iLIGOSensitivity.txt','-ascii');
freqs = y(:,1);
sqrtPSD = y(:,2);


%Freq
% Fs = 2048;
% num=3;
% T_sig=54;

N = num*T_sig*Fs; %Total number of Time Samples

T = N/Fs;
timeVec = (0:(N-1))/Fs;

% fvec = freqs(1):(1/T):freqs(end);
fvec = 0:(1/T):(Fs/2);

%Interpolation
interPSD = interp1(freqs,sqrtPSD, fvec);

% kNyq = floor(N/2)+1;
% negFStrt = 1-mod(nSamples,2);
% psdVectotal = [interPSD,psdVals((kNyq-negFStrt):-1:2)];

% Modifications
minidx = find(fvec==50);
maxidx = find(fvec==700);

Sn50 = interPSD(minidx);
Sn700 = interPSD(maxidx);

interPSD(1:minidx) = Sn50;
interPSD(maxidx:end) = Sn700;

%Construct two-sided PSD
% kNyq = floor(N/2)+1;
% negFStrt = 1-mod(N,2);
% TotalpsdVec = [interPSD((kNyq-negFStrt):-1:2), interPSD];
% totfvec = [fvec((kNyq-negFStrt):-1:2), fvec];

% loglog(totfvec,TotalpsdVec);

% figure;
% plot(fvec,interPSD);
% xlabel('Frequency (Hz)');
% ylabel('PSD');
% xlim([0,2048]);
% title("Original PSD (Two-sided)")
% 
%% Make colored Noise
fltrOrdr = 500;
% 
outNoise_t = statgaussnoisegen(N+fltrOrdr,[fvec(:),interPSD(:)],fltrOrdr,Fs, noise_num);

outNoise = outNoise_t(fltrOrdr+1:end);
% plot(timeVec,outNoise);
% 
% %Estimate PSD using Welch's Method
% [pxx,f]=pwelch(outNoise, 256,[],[],Fs);
% figure;
% plot(f,pxx);
% xlabel('Frequency (Hz)');
% ylabel('PSD');
% xlim([0,2048]);
% title("Estimated PSD (One-sided)")
% % Plot the colored noise realization
% % figure;     
% % plot(timeVec,outNoise);
% % xlabel("Time");
% % ylabel("Magnitude");
% % title("Colored Noise");
% % figure;
% % histogram(outNoise);
% % title("Histogram of Colored Noise");



