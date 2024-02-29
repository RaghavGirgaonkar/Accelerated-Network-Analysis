function wave = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0, tau1p5,datalen,initial_phase,snr,N, avec, normfac)
%Creates and Returns normalized phase vector of waveform in the Fourier Domain
%General expression for a 2PN waveform in a fourier domain is 
% W(f) = A(f)exp(-iPsi(f)), this function returns exp(-iPsi(f))
% Waveform generation happens through chirp-time values tau0, tau1.5 

fwavepos = waveform_tau(fpos,ta,phase,fmin,fmax,tau0,tau1p5,datalen,initial_phase, avec);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];

wave = snr*normfac*fwave;



