function []=psooutput(outStruct,maxSteps, nRuns, dataX, dataY, tau0, tau1p5, m1, m2, snr, ta, phase, original_fitVal, bestFitVal, c, G, type)
%OUTPUT Display PSO estimated parameters and make plots
%   Input: outStruct: Structure returned by PSO
%          maxSteps: Maximum iterations of PS)
%          nRuns: Number of independent PSO runs
%          dataX: time vector
%          dataY: data vector
%          tau0, tau1p5: Chirp time parameters (initial)
%          m1, m2: Component masses (initial)
%          snr: Initial specified snr
%          ta: Initial time of arrival
%          phase: Initial Phase
%          original_fitVal: Fitness value at the true parameter location
%          bestFitVal: Final best run fitness value found by PSO
%          c, G: Constants
%          type: Flag to specify Chirp time (Mass) space: 1 (0)
%   Output: 1. Overlay plot of data realization and estimated signal from
%              best run
%           2. Plot of Evolution of Fitness Values for all runs over iteration
%           3. Plot of Evolution of Best parameter values for all runs
%           4. Displays Original (if custom signal is injected) and
%              estimated PSO parameters 
%   Note: saveas() commands can be commented and uncommented based on need
%         of saving figures

% Raghav Girgaonkar, May 2023

figure;
hold on;
plot(dataX,dataY,'.');
% plot(dataX,wave,'r');
% for lpruns = 1:nRuns
%       plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
% end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
% legend('Data','Signal',...
%         ['Estimated signal: ',num2str(nRuns),' runs'],...
%         'Estimated signal: Best run');
legend('Data','Signal',...
        'Estimated signal: Best run');
% saveas(gcf,files.psoresultplot);
hold off;

figure;
iterVec = linspace(1,maxSteps,maxSteps);
hold on;
for lpruns = 1:nRuns
      plot(iterVec,outStruct.allRunsOutput(lpruns).allBestFit, 'DisplayName',num2str(lpruns));
end
title("Best Fitness Values for All Runs");
xlabel("Iteration");
ylabel("Best Fitness Value");
legend;
% saveas(gcf,files.bestfitplot);
hold off;

if type
    figure;
    hold on;
    for lpruns = 1:nRuns
          rVec = s2rv(outStruct.allRunsOutput(lpruns).allBestLoc,inParams);
          plot(rVec(:,1),rVec(:,2),'DisplayName',num2str(lpruns));
    end
    scatter(tau0,tau1p5,140,'red','filled','D','DisplayName','Original Parameters');
    title("Best Parameter Values for All Runs");
    xlabel("\tau_0");
    ylabel("\tau_{1.5}");
    legend;
    boundary_plot;
%     saveas(gcf,files.bestlocplot);
    hold off;

    

    t0 = outStruct.bestQcCoefs(1);
    t1p5 = outStruct.bestQcCoefs(2);
    est_M = (5/(32*fmin))*(t1p5/(pi*pi*t0))*(c^3/G);
    est_u = (1/(16*fmin*fmin))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(c^3/G);
    
    est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
    est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;
    
    %% This will display parameters given through signal.json and PSO-estimated parameters
    %% Uncomment Original parameter display command if needed
    disp(['Original parameters: tau0= ',num2str(tau0),...
                                  '; tau1p5= ',num2str(tau1p5),...
                                  '; m1= ', num2str(m1/Msolar),...
                                  '; m2= ', num2str(m2/Msolar),...
                                  '; A = ',num2str(snr),...
                                 '; phi = ',num2str(phase),...
                                  '; t_a = ',num2str(ta),...
                                  '; FitVal = ',num2str(original_fitVal)]);
    
    disp(['Estimated parameters: tau0=',num2str(outStruct.bestQcCoefs(1)),...
                                  '; tau1p5=',num2str(outStruct.bestQcCoefs(2)),...
                                  '; m1= ', num2str(est_m1/Msolar),...
                                  '; m2= ', num2str(est_m2/Msolar),...
                                  '; A = ',num2str(outStruct.bestAmp),...
                                 '; phi = ',num2str(outStruct.bestPhase),...
                                  '; t_a = ',num2str(outStruct.bestTime),...
                                  '; FitVal = ',num2str(bestFitVal)]);
else
    figure;
    hold on;
    for lpruns = 1:nRuns
          rVec = s2rv(outStruct.allRunsOutput(lpruns).allBestLoc,inParams);
          plot(rVec(:,1),rVec(:,2),'DisplayName',num2str(lpruns));
    end
    scatter(m1,m2,140,'red','filled','D','DisplayName','Original Parameters');
    title("Best Parameter Values for All Runs");
    xlabel("m_1");
    ylabel("m_2");
    legend;
%     saveas(gcf,files.bestlocplot);
    hold off;
    
    %% This will display parameters given through signal.json and PSO-estimated parameters
    %% Uncomment Original parameter display command if needed
    disp(['Original parameters:  m1= ',num2str(m1),...
                                  '; m2= ',num2str(m2),...
                                  '; A = ',num2str(snr),...
                                 '; phi = ',num2str(phase),...
                                  '; t_a = ',num2str(ta),...
                                  '; FitVal = ',num2str(original_fitVal)]);

    disp(['Estimated parameters: m1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; m2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; A = ',num2str(outStruct.bestAmp),...
                             '; phi = ',num2str(outStruct.bestPhase),...
                              '; t_a = ',num2str(outStruct.bestTime),...
                              '; FitVal = ',num2str(bestFitVal)]);
end

end
