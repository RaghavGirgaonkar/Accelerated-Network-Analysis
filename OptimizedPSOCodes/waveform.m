function fwave = waveform(fvec, t, phase, fmin, fmax,m1, m2,datalen,initial_phase, avec)
%Function to create Restricted 2PN Waveform in Fourier Domain
%Input: fvec = Row vector of positive DTFT frequencies 
%       t = Time of arrival in seconds
%       phase = Phase of Coalescence
%       fmin, fmax = minimum and maximum frequency cutoffs for the waveform
%       m1, m2 = Component masses
%       datalen = Length of data realization in seconds
%       avec = Pre-calculated alpha term vector
%Output: fwave = Phase term of waveform in the Fourier domain

%Raghav Girgaonkar, April 2023

%% Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
%% Calculate Mass Terms
m1 = m1*Msolar;
m2 = m2*Msolar;
M = (m1 + m2);
u = m1*m2/M;
n = u/M;
%% Calculate Chirp Times
tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);

tau1 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1))*(1/n)*((743/336)+ (11*n/4));

tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);

tau2 = (5/(128*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1/3))*(1/n)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));

%% Alpha Terms
alpha0 = 2*pi*fmin*(3*tau0/5);

alpha1 = 0;

alpha2 = 2*pi*fmin*tau1;

alpha3 = -2*pi*fmin*(3*tau1p5/2);

alpha4 = 2*pi*fmin*3*tau2;

alphaTerm = zeros(size(fvec));

alphaTerm(2:end) = alpha0*avec(1,:)... 
+ alpha1*avec(2,:)...
+ alpha2*avec(3,:)... 
+ alpha3*avec(4,:)... 
+ alpha4*avec(5,:);

F = 2*pi*fvec*(tau0 + tau1 - tau1p5 + tau2);

%% Final Phase Term

Psi = 2*pi*t*fvec - phase - pi/4 + alphaTerm + F + initial_phase;


%% Final Expression
P = exp(-1j*Psi);
fwave = P;

%% Cut between fmin and fmax
min_index  = floor(datalen*fmin) + 1;
max_index = floor(datalen*fmax) + 1;

fwave(1:min_index-1) = 0;
fwave(max_index + 1: end) = 0;





