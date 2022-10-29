function fwave = waveform(fvec, t, phase, N,fmin, m1, m2)
%Function to create Restricted 2PN Waveform in Fourier Domain
%Input is total frequency vector
%Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
%Calculate Mass Terms
m1 = m1*Msolar;
m2 = m2*Msolar;
M = (m1 + m2);
u = m1*m2/M;
n = u/M;
%Calculate Chirp Times
tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);

tau1 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1))*(1/n)*((743/336)+ (11*n/4));

tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);

tau2 = (5/(128*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1/3))*(1/n)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));

%Alpha Terms
alpha0 = 2*pi*fmin*(3*tau0/5);

alpha1 = 0;

alpha2 = 2*pi*fmin*tau1;

alpha3 = -2*pi*fmin*(3*tau1p5/2);

alpha4 = 2*pi*fmin*3*tau2;

alpha = [alpha0, alpha1, alpha2, alpha3, alpha4];

A = zeros(size(fvec));

A(2:end) = fvec(2:end).^(-7/6);

% A(1:10)

alphaTerm = zeros(size(fvec));

for k = 0:4
    alphaTerm(2:end) = alphaTerm(2:end) + alpha(k+1)*(fvec(2:end)./fmin).^((-5+k)/3);

end

%Construct the Negative Spectrum
% fneg = -1*fvec(end-1:-1:2);

% Aneg = A(end-1:-1:2)*-1;

% Aneg(end-9:end)

% alphaTermneg = alphaTerm(end-1:-1:2)*-1;

%Final Phase Term

Psi = 2*pi*t*fvec - phase - pi/4 + alphaTerm;

% Psineg = 2*pi*t*fneg - phase - pi/4 + alphaTermneg;

%Final Expression

fwavepos = A.*exp(-1*1j*Psi);
% fwaveneg = Aneg.*exp(-1*1i*Psineg);
if mod(N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];




