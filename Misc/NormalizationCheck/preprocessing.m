function [A,avec, phaseDiff] = preprocessing(fmin,fmax,fpos, datalen, N)
%PREPROCESSING Function to create vectors used for waveform generation
%   Input: fmin: Low frequency cutoff for waveform 
%          fmax: High frequency cutoff for waveform
%          fpos: Positive frequency vector
%          datalen: total length of data in seconds
%          N: Total number of samples in data segment
%   Output: A: Frequency magnitude vector 
%           avec: matrix containing alpha terms used for waveform
%           generation
%           phaseDiff: phase difference vector for quadrature templates

%Raghav Girgaonkar, May 2023


%Make Frequency Magnitude and Phase difference vectors
Apos = zeros(size(fpos));
phaseDiffpos = -1j*ones(size(fpos));

Apos(2:end) = fpos(2:end).^(-7/6);
%Modify Apos
min_index  = floor(datalen*fmin) + 1;
max_index = floor(datalen*fmax) + 1;
Apos(1:min_index-1) = 0;
Apos(max_index + 1: end) = 0;
%Make Aneg and Phasediffneg
if mod(N,2) == 0
    Aneg = conj(Apos(end-1:-1:2));
    phaseDiffneg = conj(phaseDiffpos(end-1:-1:2));
else
    Aneg = conj(Apos(end:-1:2));
    phaseDiffneg = conj(phaseDiffpos(end:-1:2));
end
%Make full A
A = [Apos, Aneg];
phaseDiff = [phaseDiffpos, phaseDiffneg];

a0fvec = ((fpos(2:end)./fmin).^(-5/3));
a1fvec = ((fpos(2:end)./fmin).^(-4/3));
a2fvec = ((fpos(2:end)./fmin).^(-3/3));
a3fvec = ((fpos(2:end)./fmin).^(-2/3));
a4fvec = ((fpos(2:end)./fmin).^(-1/3));

avec = [a0fvec ; a1fvec; a2fvec; a3fvec; a4fvec];
end

