addpath('../AllPSO/');

ngridpoints = 2000;
nsectors = 16;

nindex = ngridpoints/sqrt(nsectors);

Fmatrix = zeros(ngridpoints);

for i = 1:nindex:ngridpoints
    
    for j = 1:nindex:ngridpoints
            S = load(['/Users/raghav/Documents/TACC_Output/FitValImg/MaxVal/files' num2str(j) num2str(i) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
            M = abs(S.fitvals);
            M = M';
            Fmatrix(i:i+nindex-1, j:j+nindex-1) = M;
    end
     

end

[i,j] = find(Fmatrix==max(max(Fmatrix)));

x = linspace(0, 0.8, ngridpoints);
y = linspace(0, 0.1, ngridpoints);

tau0 = 3.9998;
tau1p5 = 0.52173;
% % 
figure;
hold on;
imagesc(x, y, Fmatrix);
% SFP = surf(x,y,Fmatrix, Fmatrix);
% set(SFP,'LineStyle','none');

axis xy;
xlim([0 0.8])
ylim([0 0.1])
xlabel("\tau_0");
ylabel("\tau_{1.5}");
title("Fitness Values: Tau Space with QC")
% scatter(tau0,tau1p5,70,'red','filled','D','DisplayName','Original Parameters');
scatter(x(j),y(i),70,'green','filled','D','DisplayName','Maximum Fitness Value');
% boundary_plot;
colorbar;
legend;
hold off;
