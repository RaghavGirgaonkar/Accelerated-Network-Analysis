function signal = gen2PNtemplate_mass(fpos, ta, phase, fmin, fmax,m1,m2,datalen,initial_phase,snr,N, A, avec, PSDtotal)
%Creates and Returns normalized phase vector of waveform in the Fourier Domain
%General expression for a 2PN waveform in a fourier domain is 
% W(f) = A(f)exp(-iPsi(f)), this function returns exp(-iPsi(f))
% Waveform generation happens through chirp-time values tau0, tau1.5 

fwavepos = waveform(fpos,ta,phase,fmin,fmax,m1,m2,datalen,initial_phase, avec);

if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end



fwave = [fwavepos, fwaveneg];
wf = ifft(A.*fwave);

%Normalize to unit 1
normfac = 1/sqrt(innerproduct(wf,wf,PSDtotal));

% Create final signal
signal = snr*normfac*wf;
