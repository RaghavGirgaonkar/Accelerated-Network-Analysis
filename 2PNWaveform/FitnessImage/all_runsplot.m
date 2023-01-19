%% Make Tau0 and Tau1.5 Evolution Plots for 50 noise realizations
addpath('../AllPSO/');
nRuns = 8;
maxSteps = 500;
iterVec = linspace(1,maxSteps,maxSteps);
rmin = [0, 0];
rmax = [90, 2];

inParams = struct('rmin',rmin,...
                  'rmax',rmax);

tau0_evolution = zeros(50,500);
tau1p5_evolution = zeros(50,500);
tau0 = 3.9998*ones(1,500);
tau1p5= 0.52173*ones(1,500);

figure;
hold on;
for j = 1:50
        S = load(['/Users/raghav/Documents/MassvTauPSO/TauPSO/SNR_12_m1_5_m2_9_ta_40_iter_2/outstruct_tau_' num2str(j) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
        bestRun = S.outStruct.bestRun;
        rVec = s2rv(S.outStruct.allRunsOutput(bestRun).allBestLoc,inParams);
        tau0_evolution(j,:) = rVec(:,1)';
        tau1p5_evolution(j,:) = rVec(:,2)';
        plot(iterVec,rVec(:,1), 'HandleVisibility','off');
end

mean_tau0 = mean(tau0_evolution);
mean_tau1p5 = mean(tau1p5_evolution);
plot(iterVec,mean_tau0,DisplayName='Mean', Color='black', LineWidth=5);
plot(iterVec,tau0,DisplayName='Original Value', Color='cyan', LineWidth=2);
xlabel('Iteration');
ylabel('\tau_0')
title('Best Run Evolution of \tau_0 over 50 Noise realizations');
legend;
hold off;

figure;
hold on;
for j = 1:50
        S = load(['/Users/raghav/Documents/MassvTauPSO/TauPSO/SNR_12_m1_5_m2_9_ta_40_iter_2/outstruct_tau_' num2str(j) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
        bestRun = S.outStruct.bestRun;
        rVec = s2rv(S.outStruct.allRunsOutput(bestRun).allBestLoc,inParams);
        plot(iterVec,rVec(:,2),'HandleVisibility','off');
end

plot(iterVec,mean_tau1p5,DisplayName='Mean', Color='black', LineWidth=5)
plot(iterVec,tau1p5,DisplayName='Original Value', Color='cyan', LineWidth=2);
xlabel('Iteration');
ylabel('\tau_{1.5}')
title('Best Run Evolution of \tau_{1.5} over 50 Noise realizations');
legend;
hold off;

%% Fitness Value Plots
all_fitval = zeros(50,500);
figure;
hold on;
for j = 1:50
        S = load(['/Users/raghav/Documents/MassvTauPSO/TauPSO/SNR_12_m1_5_m2_9_ta_40_iter_2/outstruct_tau_' num2str(j) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
        bestRun = S.outStruct.bestRun;
        rVec = S.outStruct.allRunsOutput(bestRun).allBestFit;
        all_fitval(j,:) = rVec;
        plot(iterVec,rVec, 'HandleVisibility','off');
end

mean_fit = mean(all_fitval);
plot(iterVec,mean_fit,DisplayName='Mean', Color='black', LineWidth=5);
% plot(iterVec,tau0,DisplayName='Original Value', Color='cyan', LineWidth=2);
xlabel('Iteration');
ylabel('Fitness Value')
title('Best Run Evolution of Fitness Value over 50 Noise realizations');
legend;
hold off;

