function wave = gen2PNwaveform(fpos, ta, phase, fmin, fmax,m1, m2,datalen,initial_phase,snr,N,avec, normfac)
%Returns normalized phase vector of waveform in the Fourier Domain

fwavepos = waveform(fpos,ta,phase,fmin,fmax,m1,m2,datalen,initial_phase, avec);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];

wave = snr*normfac*fwave;


