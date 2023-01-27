a = 1:10;
N = 100;
Fs = 1;
A = zeros( ...
    1,N);

A(40:60) = 1;

% for i = 1:6
%     sum = 0;
%     for j = 0:9
%         sum = sum + a(j+1)*exp(-2*pi*1i*(i-1)*j/10);
%     end
%     A(i) = sum;
% end
% 
% % Aneg = zeros(1,4);
% 
% Aneg = conj(A(5:-1:2));
% 
% finalA = [A, Aneg]

fftA = fft(A);

datalen = N/Fs;
fpos = (0:floor(N/2))*(1/datalen);
% fneg = fpos(end-1:-1:2)*-1;


fftA(1:floor(100/2)+1) = fftA(1:floor(100/2)+1).*exp(1*1j*2*pi*fpos*10);
fftA(floor(100/2)+2:end) = conj(fftA(floor(100/2):-1:2));

plot(A);
hold on;
b = ifft(fftA);
plot(real(b));

