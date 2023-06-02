% % 
% massf = [0.72, 0.46, 0.54, 0.82]; tauf = [1,1,1,1]; snr = [6,8,10,20];
% 
% figure; 
% hold on; 
% plot(snr, massf, '-x', DisplayName='Mass PSO', Color='red'); 
% plot(snr, tauf, '-x', DisplayName='Tau PSO', Color='black'); 
% ylim([0,1.2]); 
% xlabel('SNR'); 
% ylabel("$\mathcal{F}$",'interpreter','latex'); 
% legend; 
% hold off;

%%%%%%% Comparison in PSO-Found Parameter Locations %%%%%%%%%%%%%%%%%%
% n = 50;
% fmin = 30;
% c = 3*10^8;
% Msolar = 1.989*10^30;
% G = 6.6743*10^-11;
% 
% m1 = 3;
% m2 = 27;
% 
% tau_dir = '/Users/kqa493/Documents/Corrected_GenMassvTauPSO_Output/TauPSO/Gen_SNR_20_3_27_90/outstruct_tau_';
% mass_dir = '/Users/kqa493/Documents/Corrected_GenMassvTauPSO_Output/MassPSO/Gen_SNR_20_3_27_90/outstruct_mass_';
% 
% plot_dir = '/Users/kqa493/Desktop/MassvTau-Report/sig3_20_comp.png';
% 
% % mass_dir = [mass_File_dir 'outstruct_mass_'];
% % tau_dir = [tau_File_dir 'outstruct_tau_'];
% 
% %%% Find All Tau Locations
% 
% all_tau0 = zeros(1,50);
% all_tau1p5 = zeros(1,50);
% % all_tas = zeros(1,50);
% for j = 1:n
% %     if j ==24
% %         continue;
% %     end
%         S = load([tau_dir num2str(j) '.mat']);
% %             disp(['files' num2str(j) num2str(i) '.mat']);
%         a = S.outStruct.bestQcCoefs;
% %         t = S.outStruct.bestTime;
%         all_tau0(j) = a(1);
%         all_tau1p5(j) = a(2);
% %         all_tas(j) = t;
% %         scatter(a(1),a(2), 'HandleVisibility','off');
%        
% end
% 
% %%% Find all Mass Locations
% 
% all_m1 = zeros(1,50);
% all_m2 = zeros(1,50);
% % all_tas = zeros(1,50);
% 
% for j = 1:50
% 
%         S = load([mass_dir num2str(j) '.mat']);
% %             disp(['files' num2str(j) num2str(i) '.mat']);
%         a = S.outStruct.bestQcCoefs;
% %         t = S.outStruct.bestTime;
% %         if a(1) < a(2)
% %                 all_m1(j) = a(1);
% %                 all_m2(j) = a(2);
% %         else
% %                 all_m2(j) = a(1);
% %                 all_m1(j) = a(2);
% %         end
%         all_m1(j) = a(1);
%         all_m2(j) = a(2);
% %         all_tas(j) = t;
% %         scatter(a(1),a(2), 'HandleVisibility','off');
%        
% end
% 
% %%% Convert the Tau Found Locations onto Mass found locations 
% conv_m1 = zeros(1,50);
% conv_m2 = zeros(1,50);
%  
% for j = 1:50
%     t0 = all_tau0(j);
%     t1p5 = all_tau1p5(j);
%     est_M = (5/(32*fmin))*(t1p5/(pi*pi*t0))*(c^3/G);
%     est_u = (1/(16*fmin*fmin))*(5/(4*pi^4*t0*t1p5^2))^(1/3)*(c^3/G);
%     
%     est_m1 = (est_M - sqrt(est_M^2 - 4*est_u*est_M))/2;
%     est_m2 = (est_M + sqrt(est_M^2 - 4*est_u*est_M))/2;
% 
%     conv_m1(j) = real(est_m1/Msolar);
%     conv_m2(j) = real(est_m2/Msolar);
% end
% 
% 
% %% Scatter Plot
% 
% % figure;
% % hold on;
% % scatter(all_m1, all_m2, 'filled', DisplayName='Mass PSO Found Locations', MarkerFaceColor='red');
% % scatter(conv_m1, conv_m2, 'filled',DisplayName='Tau PSO Found Locations', MarkerFaceColor='black');
% % scatter(m1,m2,200,'pentagram',DisplayName='Original Parameters', MarkerFaceColor='cyan', MarkerEdgeColor='black');
% % scatter(m2,m1,200,'pentagram','HandleVisibility','off', MarkerFaceColor='cyan', MarkerEdgeColor='black');
% % xlabel('m_1');
% % ylabel('m_2');
% % legend;
% % hold off;
% 
% figure;
% hold on;
% scatterhist(all_m1, all_m2,'NBins',[50,50],'Color','red');
% scatterhist(conv_m1, conv_m2,'NBins',[50,50],'Color','black');
% % scatter(conv_m1, conv_m2, 'filled',DisplayName='Tau PSO Found Locations', MarkerFaceColor='black');
% % scatter(m1,m2,200,'pentagram',DisplayName='Original Parameters', MarkerFaceColor='cyan', MarkerEdgeColor='black');
% % scatter(m2,m1,200,'pentagram','HandleVisibility','off', MarkerFaceColor='cyan', MarkerEdgeColor='black');
% xlabel('m_1');
% ylabel('m_2');
% legend;
% hold off;

% saveas(gcf,plot_dir);



%%%%% Median Values for Parameters Found

tau_dir = '/Users/kqa493/Documents/Corrected_GenMassvTauPSO_Output/TauPSO/Gen_SNR_20_20_26_36/outstruct_tau_';
mass_dir = '/Users/kqa493/Documents/Corrected_GenMassvTauPSO_Output/MassPSO/Gen_SNR_20_20_26_36/outstruct_mass_';

n=50;

snr = 8;

%% Tau Space
all_tau0 = zeros(1,50);
all_tau1p5 = zeros(1,50);
all_tas = zeros(1,50);
all_snrs = zeros(1,50);
all_phases = zeros(1,50);
for j = 1:n
%     if j ==24
%         continue;
%     end
        S = load([tau_dir  num2str(j) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
        a = S.outStruct.bestQcCoefs;
        t = S.outStruct.bestTime;
        p = S.outStruct.bestPhase;
        A = S.outStruct.bestAmp;
        all_tau0(j) = a(1);
        all_tau1p5(j) = a(2);
        all_tas(j) = t;
        all_snrs(j) = A;
        all_phases(j) = p;       
end

%% Mass Space
all_m1 = zeros(1,50);
all_m2 = zeros(1,50);
all_mtas = zeros(1,50);
all_msnrs = zeros(1,50);
all_mphases = zeros(1,50);

for j = 1:50

        S = load([mass_dir  num2str(j) '.mat']);
%             disp(['files' num2str(j) num2str(i) '.mat']);
        a = S.outStruct.bestQcCoefs;
        t = S.outStruct.bestTime;
        p = S.outStruct.bestPhase;
        A = S.outStruct.bestAmp;
        if snr > 6
            if a(1) < a(2)
                all_m1(j) = a(1);
                all_m2(j) = a(2);
            else
                all_m2(j) = a(1);
                all_m1(j) = a(2);
            end
        else
            all_m1(j) = a(1);
            all_m2(j) = a(2);
        end
        all_mtas(j) = t;
        all_msnrs(j) = A;
        all_mphases(j) = p;
       
end

med_tau0 = median(all_tau0);
med_tau1p5 = median(all_tau1p5);
med_tau_ta = median(all_tas);
med_tau_snrs = median(all_snrs);
med_tau_phases = median(all_phases);


med_m1 = median(all_m1);
med_m2 = median(all_m2);
med_mass_ta = median(all_mtas);
med_mass_snrs = median(all_msnrs);
med_mass_phases = median(all_mphases);


disp('Tau PSO');
disp(['Median parameters: tau0= ',num2str(med_tau0),...
                                  '; tau1p5= ',num2str(med_tau1p5),...
                                  '; A = ',num2str(med_tau_snrs),...
                                 '; phi = ',num2str(med_tau_phases),...
                                  '; t_a = ',num2str(med_tau_ta)]);
disp(' ');
disp('Mass PSO');

disp(['Median parameters:  m1= ',num2str(med_m1),...
                                  '; m2= ',num2str(med_m2),...
                                  '; A = ',num2str(med_mass_snrs),...
                                 '; phi = ',num2str(med_mass_phases),...
                                  '; t_a = ',num2str(med_mass_ta)]);

























