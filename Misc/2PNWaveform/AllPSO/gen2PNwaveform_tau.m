function wave = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0, tau1p5,datalen,initial_phase,snr,N,PSD)

fwavepos = waveform_tau(fpos,ta,phase,fmin,fmax,tau0,tau1p5,datalen,initial_phase);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];

wave = ifft(fwave);

sampFreq = N/datalen;

[wave,~] = normsig4psd(wave,sampFreq,PSD,snr);
% 
% wave = snr*wave/norm(wave);


