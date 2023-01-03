addpath('../AllPSO/');

ngridpoints = 2000;
nsectors = 16;

nindex = ngridpoints/sqrt(nsectors);

Fmatrix = zeros(ngridpoints);

for i = 1:nindex:ngridpoints
    
    for j = 1:nindex:ngridpoints
            S = load(['files' num2str(j) num2str(i) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
            M = S.fitvals;
            M = M';
            Fmatrix(i:i+nindex-1, j:j+nindex-1) = M;
    end
     

end

[i,j] = find(Fmatrix==max(max(Fmatrix)));

x = linspace(0.1, 70, ngridpoints);
y = linspace(0.1, 2, ngridpoints);

tau0 = 11.9050;
tau1p5 = 0.7764;

figure;
hold on;
imagesc(x, y, Fmatrix);
% SFP = surf(x,y,Fmatrix, Fmatrix);
% set(SFP,'LineStyle','none');

axis xy;
xlabel("\tau_0");
ylabel("\tau_{1.5}");
title("Fitness Values: Tau Space (Without Noise)")
% scatter(tau0,tau1p5,70,'red','filled','D','DisplayName','Original Parameters');
scatter(x(j),y(i),70,'green','filled','D','DisplayName','Maximum Fitness Value');
boundary_plot;
colorbar;
legend;
hold off;
