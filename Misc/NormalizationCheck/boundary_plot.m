function [boundary_fig] = boundary_plot
% Creates and returns a Boundary Plot as a figure in the Chirp-time (Tau1.5 and Tau0) space
% Search range of component masses is 1.4 to 30 Solar Masses

%Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
%Frequency min max values
fmin = 30;
fmax = 700;
%Mass vectors
m_lin = linspace(1.4,30,10000);
m1 = m_lin;
m2 = m1;
%Calculate Mass Terms
m1 = m1*Msolar;
m2 = m2*Msolar;
M = (m1 + m2);
u = m1.*m2./M;
n = u./M;
%Calculate Chirp Times firt when m1 = m2
tau0_1 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3).^(-5/3)).*(1./n);

% tau1 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1))*(1/n)*((743/336)+ (11*n/4));

tau1p5_1 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3).^(-2/3)).*(1./n);

% tau2 = (5/(128*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1/3))*(1/n)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));

%Boundary Array for m2 = 1.4
%Mass vectors
m1 = m_lin;
m2 = m1(1);
%Calculate Mass Terms
m1 = m1*Msolar;
m2 = m2*Msolar;
M = (m1 + m2);
u = (m1*m2)./M;
n = u./M;
tau0_2 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3).^(-5/3)).*(1./n);
tau1p5_2 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3).^(-2/3)).*(1./n);

%Boundary Array for m1 = 30
%Mass vectors
m2 = m_lin;
m1 = m2(end);
%Calculate Mass Terms
m1 = m1*Msolar;
m2 = m2*Msolar;
M = (m1 + m2);
u = (m1*m2)./M;
n = u./M;
tau0_3 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3).^(-5/3)).*(1./n);
tau1p5_3 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3).^(-2/3)).*(1./n);

%Final Tau0 and Tau1.5 arrays

Tau0 = [tau0_1, tau0_2, tau0_3];

Tau1p5 = [tau1p5_1, tau1p5_2, tau1p5_3];

%Scatterplot
sz = 5;
c = 'black';
% boundary_fig = figure;
% hold on;
% t1 = scatter(tau0_1, tau1p5_1,sz, 'red', 'filled', 'HandleVisibility','off');
% hold on;
% t2 = scatter(tau0_2, tau1p5_2,sz, c, 'filled', 'HandleVisibility','off');
% hold on;
% t3 = scatter(tau0_3, tau1p5_3,sz, c, 'filled', 'HandleVisibility','off');
% hold off;
boundary_fig = scatter(Tau0, Tau1p5,sz,c,"filled", 'DisplayName','Tau Space Boundary');
% xlabel('\tau_0');
% ylabel('\tau_{1.5}');
% title("Boundary of Region of Physical Parameters in \tau_0 and \tau_{1.5} Space");
% saveas(gcf,"boundary-tauspace.pdf");