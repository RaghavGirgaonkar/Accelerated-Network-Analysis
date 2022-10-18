function timesVec = matchedfiltering(y,q)
%function to perform FFT correlation based matched filtering using a signal
%y = signal + noise and q a unit normalised template starting at t = 0,
%assumes two-sided psd is provided
%returns ta, the best matching time of the template q in the signal y and
%the value of the likelihood at that time.



Y = fft(y);
Z = Y;
Q = fft(q);


timesVec = ifft((Z.*conj(Q)));


% plot(real(timesVec));
% [max_val, max_sample] = max(real(timesVec));
% fprintf('The max value of t_a is = %f at time = %f\n',max_val,(1/sampling_freq)*max_sample);
% ta = (1/sampling_freq)*max_sample;