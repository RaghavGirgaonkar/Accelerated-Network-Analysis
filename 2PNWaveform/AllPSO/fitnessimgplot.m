one = readmatrix("1-1000_1-1000.txt");
two = readmatrix("1001-2000_1-1000.txt");
three = readmatrix("1-1000_1001-2000.txt");
four = readmatrix("1001-2000_1001-2000.txt");

ngridpoints = 2000;

one = one(end:-1:1,1:end);
two = two(end:-1:1,1:end);
three = three(end:-1:1,1:end);
four = four(end:-1:1,1:end);

z = [three, four; one, two];



x = linspace(0.1, 70, ngridpoints);
y = linspace(0.1, 2, ngridpoints);

figure;
% hold on;
imagesc(x, y, z);
axis xy;
xlabel("\tau_0");
ylabel("\tau_{1.5}");
title("Fitness Values: Tau Space")
% boundary_plot;
colorbar;
% hold off;
