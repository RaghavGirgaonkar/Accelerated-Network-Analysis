function fwave = newtonian_waveform(fvec, t, phase, fmin, fmax,m1, m2,datalen,initial_phase)
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

chirp_mass = (u^3*M^2)^(1/5);

%Zi

zi = 34.54*((chirp_mass/Msolar)^(-5/3))*(fmin/40)^(-8/3);

%Calculate Chirp Times
% tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
% 
% tau1 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1))*(1/n)*((743/336)+ (11*n/4));
% 
% tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);
% 
% tau2 = (5/(128*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1/3))*(1/n)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));
% 
% %Alpha Terms
% alpha0 = 2*pi*fmin*(3*tau0/5);
% 
% alpha1 = 0;
% 
% alpha2 = 2*pi*fmin*tau1;
% 
% alpha3 = -2*pi*fmin*(3*tau1p5/2);
% 
% alpha4 = 2*pi*fmin*3*tau2;
% 
% alpha = [alpha0, alpha1, alpha2, alpha3, alpha4];

A = zeros(size(fvec));

A(2:end) = (fvec(2:end)./fmin).^(-7/6);



alphaTerm = zeros(size(fvec));

alphaTerm(2:end) = (1/5)*(8 - 3*(fvec(2:end)./fmin).^(-5/3) -5*(fvec(2:end)./fmin));
alphaTerm(1) = 8/5;

% 
% for k = 0:4
%     alphaTerm(2:end) = alphaTerm(2:end) + alpha(k+1)*(fvec(2:end)./fmin).^((-5+k)/3);
% 
% end


%Final Phase Term

Psi = -2*pi*t*fvec + phase + pi/4 + 2*pi*fmin*zi*alphaTerm + initial_phase;
% Psi = 2*pi*t*fvec - phase - pi/4;

%Final Expression

fwave = A.*exp(1j*Psi);
% fwave = A;

min_index  = floor(datalen*fmin) + 1;
max_index = floor(datalen*fmax) + 1;

fwave(1:min_index-1) = 0;
fwave(max_index + 1: end) = 0;
% fwaveneg = Aneg.*exp(-1*1i*Psineg);





