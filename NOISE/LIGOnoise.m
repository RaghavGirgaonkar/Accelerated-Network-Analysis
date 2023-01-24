%Load PSD 
y = load('iLIGOSensitivity.txt','-ascii');
freqs = y(:,1);
sqrtPSD = y(:,2);


%Freq
Fs = 16384;
num=3;
T_sig=54;

N = num*T_sig*Fs; %Total number of Time Samples

T = N/Fs;
timeVec = (0:(N-1))/Fs;

% fvec = freqs(1):(1/T):freqs(end);
fvec = 0:(1/T):(Fs/2);

%Interpolation
interPSD = interp1(freqs,sqrtPSD, fvec);

% Modifications
minidx = find(fvec==50);
maxidx = find(fvec==700);

Sn50 = interPSD(minidx);
Sn700 = interPSD(maxidx);


interPSD(1:minidx) = Sn50;
interPSD(maxidx:end) = Sn700;

figure;
plot(fvec,interPSD);
xlabel('Frequency (Hz)');
ylabel('PSD');
% xlim([0,2048]);
title("Original PSD (Two-sided)")

%Make colored Noise
fltrOrdr = 500;

outNoise = statgaussnoisegen(N,[fvec(:),interPSD(:)],fltrOrdr,Fs);

outNoise_final = outNoise(fltrOrdr+1:end-fltrOrdr);

%Estimate PSD using Welch's Method
[pxx,f]=pwelch(outNoise_final, 1024,[],[],Fs);
figure;
loglog(f,pxx);
xlabel('Frequency (Hz)');
ylabel('PSD');
% xlim([0,2048]);
title("Estimated PSD (One-sided)")
% Plot the colored noise realization
% figure;     
% plot(timeVec,outNoise);
% xlabel("Time");
% ylabel("Magnitude");
% title("Colored Noise");
% figure;
% histogram(outNoise);
% title("Histogram of Colored Noise");



