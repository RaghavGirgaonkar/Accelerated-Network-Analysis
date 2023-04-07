%Load PSD 
y = load('iLIGOSensitivity.txt','-ascii');
freqs = y(:,1);
sqrtPSD = y(:,2);

P = sqrtPSD.^2;

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
interPSD = interp1(freqs,P, fvec);

% Modifications
minidx = find(fvec==50);
maxidx = find(fvec==700);

Sn50 = interPSD(minidx);
Sn700 = interPSD(maxidx);


interPSD(1:minidx-1) = Sn50;
interPSD(maxidx+1:end) = Sn700;

figure;
loglog(fvec,interPSD);
xlabel('Frequency (Hz)');
ylabel('PSD');
% xlim([0,2048]);
title("Original PSD (Two-sided)")

%Make colored Noise
fltrOrdr = 1000;

outNoise = statgaussnoisegen(N,[fvec(:),interPSD(:)],fltrOrdr,Fs);

outNoise_final = outNoise(fltrOrdr+1:end);

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
% plot(outNoise_final);
% xlabel("Time");
% ylabel("Magnitude");
% title("Colored Noise befor truncation");
% figure;     
plot(outNoise);
xlabel("Time");
ylabel("Magnitude");
title("Colored Noise");
% figure;
% histogram(outNoise);
% title("Histogram of Colored Noise");



