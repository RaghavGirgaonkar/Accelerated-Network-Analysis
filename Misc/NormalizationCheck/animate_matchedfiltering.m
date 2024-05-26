%% Script to create an animation of matched filtering timeseries.

%% Define Data Parameters
sampFreq = 4096;
datalen = 128;
N = datalen*sampFreq;
frange = [30,700];
fpos = (0:floor(N/2))*(1/datalen);
timeVec = (0:(datalen*sampFreq - 1))/sampFreq;

%% Create Colored Gaussian Noise
[noise, PSD] = LIGOnoise(N, sampFreq, 1, 'sample');


%Create Total PSD and Transfer Function
negFStrt = 1-mod(N,2);
kNyq = floor(N/2)+1;

PSDtotal = [PSD, PSD((kNyq-negFStrt):-1:2)];
TFtotal = 1./sqrt(PSDtotal);


%Pre-calculate frequency magnitude and phase difference terms
 [A,avec, phaseDiff] = preprocessing(frange(1),frange(2),fpos, datalen, N);

%Create general normalization factor N_fW
AbysqrtPSD = A.*TFtotal;
innProd = (1/N)*(AbysqrtPSD)*AbysqrtPSD';
genNormfacSqr = real(innProd);
N_fW = 1/sqrt(genNormfacSqr);

%% Generate Signal to Inject
 chirptimes = [29.6373, 1.1045];
 ta = 34;
 initial_phase = 0;
 phase = 0;
 snr = 12;

signal = gen2PNtemplate_tau(fpos, ta, phase, frange(1), frange(2),chirptimes(1), chirptimes(2),datalen,initial_phase,snr, N, A, avec, PSDtotal);

signal = signal*sqrt(sampFreq);

figure;
plot(timeVec, noise);
hold on;
plot(timeVec, 30*signal);
xlabel('Time (s)');
ylabel('Strain h(t)');
axis tight;
ax = gca(); ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
%Save whitened and normalized signal and noise
%Window first
whtndsignal = ifft((1/sqrt(sampFreq))*fft(signal).*TFtotal);

%Inject signal into data
strain = noise + signal;

%Whiten and normalize strain
rolloff = 0.5;
fmin = frange(1);
strainwin = strain.*tukeywin(length(strain),rolloff*sampFreq/length(strain))';
whtndstrain = ifft((1/sqrt(sampFreq))*fft(strainwin).*TFtotal);


%Make GIF
% Define the filename for the gif
giffilename = 'matchedfiltering_animation.gif';
time_lapse = 0.0001;
title_text = 'Generation of Matched-Filtering Timeseries';

% Create a new figure window
fig = figure;
mftimeseries = zeros(1,N);
for i = 1:1*sampFreq:N
    q0 =  gen2PNtemplate_tau(fpos, i/sampFreq, 0, frange(1), frange(2),chirptimes(1), chirptimes(2),datalen,initial_phase,1, N, A, avec, PSDtotal);
    whtndq0 = ifft(fft(q0).*TFtotal);
    q1 =  gen2PNtemplate_tau(fpos, i/sampFreq, pi/2, frange(1), frange(2),chirptimes(1), chirptimes(2),datalen,initial_phase,1, N, A, avec, PSDtotal);
    whtndq1 = ifft(fft(q1).*TFtotal);

%     I0 = innerproduct(strain, q0, PSDtotal, sampFreq);
%     I1 = innerproduct(strain, q1, PSDtotal, sampFreq);

    I0 = real((1/N)*sum(fft(whtndstrain).*conj(fft(whtndq0))));
    I1 = real((1/N)*sum(fft(whtndstrain).*conj(fft(whtndq1))));
    
    I = sqrt(I0^2 + I1^2);
    
    mftimeseries(i) = I;
    
    
    hold on;
    plot(timeVec, whtndstrain); hold on;
    plot(timeVec, 5*whtndsignal, Color='red' ); hold on;
    plot(timeVec, 20*whtndq0,Color='yellow'); hold on;
    plot(timeVec(1:i), mftimeseries(1:i), Color='black',LineWidth=1.5); hold on;
    ylim([-4,12]);
    xlim([0,128]);
    xlabel("Time (s)");
    ylabel("Whitened Strain h(t)");
    ax = gca(); ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
%     axis tight;
    
    title(title_text,'FontSize',16,'Color',[0 0 0]), hold on;
    
    drawnow;
    %Capture frame
    % Capture the current frame and save it to the gif
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
    
        if i == 1
            imwrite(imind,cm,giffilename,'gif', 'Loopcount',Inf,'DelayTime',time_lapse);
        else
            imwrite(imind,cm,giffilename,'gif','WriteMode','append','DelayTime',time_lapse);
        end
    
        clf(fig);

end

% 
% figure; 
% plot(timeVec, noise);
% hold on;
% plot(timeVec, 10*signal);
% axis tight;
