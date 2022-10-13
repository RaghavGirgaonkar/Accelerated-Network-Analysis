function ta = matchedfiltering(y,q,sampling_freq,psd)
%function to perform FFT correlation based matched filtering using a signal
%y = signal + noise and q a unit normalised template starting at t = 0,
%assumes two-sided psd is provided
%returns ta, the best matching time of the template q in the signal y.

Y = fft(y);
Z = Y./psd;
Q = fft(q);

timesVec = ifft((Z.*Q));
plot(timesVec);
[max_val, max_sample] = max(timesVec);
fprintf('The max value of t_a is = %f at time = %f\n',max_val,sampling_freq*max_sample);
ta = sampling_freq*max_sample;