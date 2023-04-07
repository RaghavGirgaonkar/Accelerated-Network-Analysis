addpath('../AllPSO/');

ngridpoints = 2000;
nsectors = 16;

%Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
fmin = 30;

m1val = 3;
m2val = 27;

m1 = m1val*Msolar;
m2 = m2val*Msolar;
M = m1 + m2;
u = m1*m2/(m1 + m2);
n = u/M;
tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);

nindex = ngridpoints/sqrt(nsectors);

Fmatrix = zeros(ngridpoints);

filedir = '/Users/raghav/Documents/FitnessImg/Sig2MassFitImg/';

for i = 1:nindex:ngridpoints
    
    for j = 1:nindex:ngridpoints
            S = load([filedir 'files' num2str(j) num2str(i) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
            M = abs(S.fitvals);
            M = M';
            Fmatrix(i:i+nindex-1, j:j+nindex-1) = M;
    end
     

end

[i,j] = find(Fmatrix==max(max(Fmatrix)));

% x = linspace(0.1, 90, ngridpoints);
% y = linspace(0.1, 2, ngridpoints);
x = linspace(1.4, 30, ngridpoints);
y = linspace(1.4, 30, ngridpoints);

% % 
figure;
hold on;
imagesc(x, y, Fmatrix);
% SFP = surf(x,y,Fmatrix, Fmatrix);
% set(SFP,'LineStyle','none');

axis xy;
% xlim([0.1 90]);
% ylim([0.1 2]);
xlim([1.4 30]);
ylim([1.4 30]);
% xlabel("\tau_0");
% ylabel("\tau_{1.5}");
xlabel("m_1");
ylabel("m_2");
% title("Fitness Values: Tau Space Signal 3, SNR 10")
title("Fitness Values: Mass Space Signal 2, SNR 10");
% scatter(tau0,tau1p5,70,'red','filled','D','DisplayName','Original Parameters');
scatter(m1val,m2val,70,'red','filled','D','DisplayName','Original Parameters');
scatter(m2val,m1val,70,'red','filled','D','HandleVisibility','off');
% scatter(x(j),y(i),70,'green','filled','D','DisplayName','Maximum Fitness Value');
% boundary_plot;
colorbar;
legend;
hold off;
saveas(gcf, '~/Desktop/Sig2MassFitImg.png');
