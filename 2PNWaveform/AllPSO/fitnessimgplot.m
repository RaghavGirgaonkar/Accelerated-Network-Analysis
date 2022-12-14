one = readmatrix("/Users/kqa493/Documents/1-1000_1-1000.txt");
two = readmatrix("/Users/kqa493/Documents/1001-2000_1-1000.txt");
three = readmatrix("/Users/kqa493/Documents/1-1000_1001-2000.txt");
four = readmatrix("/Users/kqa493/Documents/1001-2000_1001-2000.txt");

ngridpoints = 2000;

one = one';
% one = one(end:-1:1,1:end);
two = two';
% two = two(end:-1:1,1:end);
three = three';
% three = three(end:-1:1,1:end);
four = four';
% four = four(end:-1:1,1:end);

z = [one, two; three, four];

[i,j] = find(z==max(max(z)));

x = linspace(0.1, 70, ngridpoints);
y = linspace(0.1, 2, ngridpoints);

tau0 = 11.9050;
tau1p5 = 0.7764;

figure;
hold on;
imagesc(x, y, z);
axis xy;
xlabel("\tau_0");
ylabel("\tau_{1.5}");
title("Fitness Values: Tau Space")
scatter(tau0,tau1p5,70,'red','filled','D','DisplayName','Original Parameters');
scatter(x(j),y(i),70,'green','filled','D','DisplayName','Maximum Fitness Value');
boundary_plot;
colorbar;
legend;
hold off;
