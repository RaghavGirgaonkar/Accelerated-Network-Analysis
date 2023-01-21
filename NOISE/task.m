%Load .txt file

data = load("testData.txt");

timestamps = data(:,1);
datavals = data(:,2);

figure;
plot(timestamps, datavals);
xlabel("Time");
ylabel("Magnitude");
title("Colored Noise");

figure;
[S,F,T] = spectrogram(datavals,128,[],[],1/timestamps(2));
imagesc(T,F,abs(S));
axis xy;
title("Spectrogram of Colored Noise");

sig_idx = 5/timestamps(2);

noise = datavals(1:sig_idx+1);
sampFreq = 1/timestamps(2);
S = size(timestamps);
nSamples = S(1);

%Estimate PSD
[pxx,f]=pwelch(noise, 256,[],[],sampFreq);
figure;
plot(f,pxx);
xlabel('Frequency (Hz)');
ylabel('PSD');
title("Estimated PSD (One-sided)");

freqVec = f';
psdVals = pxx;

%Generate Whitening Filter
% Design FIR filter with T(f)= 1/square root of target PSD
fltrOrdr = 500;
onebysqrtPSD = 1./sqrt(psdVals);
b = fir2(fltrOrdr,freqVec/(sampFreq/2),onebysqrtPSD');

whiteneddata = sqrt(sampFreq)*fftfilt(b,datavals);

[pxx,f]=pwelch(whiteneddata(2/timestamps(2): 6/timestamps(2)), 256,[],[],sampFreq);
figure;
plot(f,pxx);
xlabel('Frequency (Hz)');
ylabel('PSD');
title("Estimated PSD 2 to 6 (One-sided)");


figure;
plot(timestamps, whiteneddata);
xlabel("Time");
ylabel("Magnitude");
title("Whitened Noise");

figure;
[S,F,T] = spectrogram(whiteneddata,128,[],[],1/timestamps(2));
imagesc(T,F,abs(S));
axis xy;
title("Spectrogram of Whitened Data");




