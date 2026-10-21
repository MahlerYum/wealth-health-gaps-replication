%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Figures in Mahler and Yum, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Load variables (empirical data from stata)
cd('..\workdata')
health10 = readmatrix('data_figures.xlsx','Sheet','health10');
eff10_edu_health = readmatrix('data_figures.xlsx','Sheet','eff10_edu_health');
wealth10_50_edu = readmatrix('data_figures.xlsx','Sheet','wealth10_50_edu');
effwealth_HS = readmatrix('data_figures.xlsx','Sheet','effwealth_HS');
effwealth_CL = readmatrix('data_figures.xlsx','Sheet','effwealth_CL');
smoking = readmatrix('data_figures.xlsx','Sheet','smoking');
sport = readmatrix('data_figures.xlsx','Sheet','sport');
nutrition = readmatrix('data_figures.xlsx','Sheet','nutrition');
survival_rates = readmatrix('data_figures.xlsx','Sheet','survival_rates');
frac_work10 = readmatrix('data_figures.xlsx','Sheet','frac_work10');
frac_work_edu10 = readmatrix('data_figures.xlsx','Sheet','frac_work_edu10');
avrg_inc10_edu_health = readmatrix('data_figures.xlsx','Sheet','avrg_inc10_edu_health');
wealth10_dat = readmatrix('data_figures.xlsx','Sheet','wealth10_dat');
wealth10_25 = readmatrix('data_figures.xlsx','Sheet','wealth10_25');
wealth10_50 = readmatrix('data_figures.xlsx','Sheet','wealth10_50');
wealth10_75 = readmatrix('data_figures.xlsx','Sheet','wealth10_75');
wealth10_50_occ = readmatrix('data_figures.xlsx','Sheet','wealth10_50_occ');

%% Load data from model simulation
cd('..')
cd('fortran\benchmark')

frac_work10_mod(1:4,2) = readmatrix('Figure5model.txt','Range','1:4');
frac_work10_mod(1:4,1) = readmatrix('Figure5model.txt','Range','5:8');
frac_work_edu10_mod(1:4,1) = readmatrix('Figure5model.txt','Range','9:12');
frac_work_edu10_mod(1:4,2) = readmatrix('Figure5model.txt','Range','13:16');
avrg_inc10_edu_health_mod(1:4,1) = readmatrix('Figure6model.txt','Range','5:8');
avrg_inc10_edu_health_mod(1:4,2) = readmatrix('Figure6model.txt','Range','1:4');
avrg_inc10_edu_health_mod(1:4,3) = readmatrix('Figure6model.txt','Range','13:16');
avrg_inc10_edu_health_mod(1:4,4) = readmatrix('Figure6model.txt','Range','9:12');
eff10_edu_health_mod(1:6,1) = readmatrix('Figure7model.txt','Range','7:12');
eff10_edu_health_mod(1:6,3) = readmatrix('Figure7model.txt','Range','1:6');
eff10_edu_health_mod(1:6,2) = readmatrix('Figure7model.txt','Range','19:24');
eff10_edu_health_mod(1:6,4) = readmatrix('Figure7model.txt','Range','13:18');
health10_mod(1:6,1) = readmatrix('Figure8model.txt','Range','7:12');
health10_mod(1:6,2) = readmatrix('Figure8model.txt','Range','13:18');
wealth10_mod(1:6,2) = readmatrix('Figure8model.txt','Range','1:6')./100;
wealth10_50_mod(2:6,2) = readmatrix('Figure9model1.txt','Range','1:5');
wealth10_50_mod(2:6,1) = readmatrix('Figure9model1.txt','Range','6:10');

whgap25_per(2:5,1) = readmatrix('Figure10model1.txt','Range','1:4')./100;
whgap25_per(2:5,2) = readmatrix('Figure10model2.txt','Range','1:4')./100;       
whgap25_per(2:5,3) = readmatrix('Figure10model3.txt','Range','1:4')./100;
whgap25_per(2:5,4) = readmatrix('Figure11model2.txt','Range','1:4');
dat = importdata('FigureA4model.txt');
whgap25_per(2:5,8) = dat(1:4)./100;

whgap50_per(2:5,1) = readmatrix('Figure10model1.txt','Range','5:8')./100;
whgap50_per(2:5,2) = readmatrix('Figure10model2.txt','Range','5:8')./100;
whgap50_per(2:5,3) = readmatrix('Figure10model3.txt','Range','5:8')./100;
whgap50_per(2:5,4) = readmatrix('Figure11model2.txt','Range','5:8');
whgap50_per(2:5,5) = readmatrix('FigureA5model1.txt','Range','5:8');
whgap50_per(2:5,6) = readmatrix('FigureA5model2.txt','Range','5:8');
whgap50_per(2:5,7) = readmatrix('FigureA5model3.txt','Range','5:8');
whgap50_per(2:5,8) = readmatrix('FigureA4model.txt','Range','5:8')./100;

whgap75_per(2:5,1) = readmatrix('Figure10model1.txt','Range','9:12')./100;
whgap75_per(2:5,2) = readmatrix('Figure10model2.txt','Range','9:12')./100;
whgap75_per(2:5,3) = readmatrix('Figure10model3.txt','Range','9:12')./100;
whgap75_per(2:5,4) = readmatrix('Figure11model2.txt','Range','9:12');
whgap75_per(2:5,8) = readmatrix('FigureA4model.txt','Range','9:12')./100;

dat = importdata('FigureA7model.txt');
wealth10_25_mod(1:6,1) = dat(7:12);
wealth10_25_mod(1:6,2) = dat(1:6);
wealth10_50_mod(1:6,1) = dat(19:24);
wealth10_50_mod(1:6,2) = dat(13:18);
wealth10_75_mod(1:6,1) = dat(31:36);
wealth10_75_mod(1:6,2) = dat(25:30);

dat = importdata('FigureA8model.txt');
wealth10_50_edu_mod(2:6,1) = dat(6:10);
wealth10_50_edu_mod(2:6,2) = dat(1:5);
wealth10_50_edu_mod(2:6,3) = dat(16:20);
wealth10_50_edu_mod(2:6,4) = dat(11:15);

% From Exoneous Health Model
cd('..\exogenous')
wealth10_50_mod2(2:6,2) = readmatrix('Figure9model2.txt','Range','1:5');
wealth10_50_mod2(2:6,1) = readmatrix('Figure9model2.txt','Range','6:10');

% Save Data for use later
cd('..\..\workdata')
save('data_fig.mat')

cd('..')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 1
hh1 = figure
subplot(1,3,1)
xx1 = plot(1:6,health10(1:6,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,health10(1:6,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Share Unhealthy')
grid on
ylim([0 0.5]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
title('Health by Education');
ll = legend([xx1 xx2],'Non-College', 'College');

subplot(1,3,2)
xx1 = plot(1:6,eff10_edu_health(1:6,1),'Color',[1 0 0],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,eff10_edu_health(1:6,3),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
hold off
grid on
ylim([0.55 0.85])
xlabel('Age Group')
ylabel('Average Health Effort')
title('Effort among Non-College')
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
ll = legend([xx2 xx1],'Healthy', 'Unhealthy');
set(ll,'Location','northwest');

subplot(1,3,3)
xx1 = plot(1:6,eff10_edu_health(1:6,2),'Color',[1 0 0],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,eff10_edu_health(1:6,4),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
hold off
grid on
ylim([0.55 0.85])
xlabel('Age Group')
ylabel('Average Health Effort')
title('Effort among College')
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
ll = legend([xx2 xx1],'Healthy', 'Unhealthy');
set(ll,'Location','southeast');

set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(gcf,'Position',[50 50 1300 500]);
saveas(gcf,'figures\Figure1','epsc2')

%% Figure 2
hh2 = figure;
subplot(1,2,1)
xx1 = semilogy(1:6,wealth10_50_edu(1:6,1),':rx','LineWidth',2,'MarkerSize',12);
hold on
xx2 = semilogy(1:6,wealth10_50_edu(1:6,2),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
title('Non-College');
subplot(1,2,2)
xx1 = semilogy(1:6,wealth10_50_edu(1:6,3),':rx','LineWidth',2,'MarkerSize',12);    
hold on
xx2 = semilogy(1:6,wealth10_50_edu(1:6,4),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
title('College');
ll = legend([xx2 xx1],'Healthy', 'Unhealthy');
set(ll,'Location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\Figure2','epsc2')

%% Figure 3
figure
subplot(1,2,1)
d11 = plot(1:4,effwealth_HS(1,:),'LineWidth',2,'Marker','x','MarkerSize',10);
hold on
d21 = plot(1:4,effwealth_HS(2,:),'LineStyle','-.','LineWidth',2,'Marker','s');
d31 = plot(1:4,effwealth_HS(3,:),'LineStyle','--','LineWidth',2,'Marker','o');
grid on
ylim([0.55 0.85])
xlabel('Wealth Quartiles')
ylabel('Average Health Effort')
title('Non-College')
xticks([1 2 3 4])
xticklabels({'1st', '2nd', '3rd', '4th'});
subplot(1,2,2)
d1 = plot(1:4,effwealth_CL(1,:),'LineWidth',2,'Marker','x','MarkerSize',10);
hold on
d2 = plot(1:4,effwealth_CL(2,:),'LineStyle','-.','LineWidth',2,'Marker','s');
d3 = plot(1:4,effwealth_CL(3,:),'LineStyle','--','LineWidth',2,'Marker','o');
grid on
ylim([0.55 0.85])
xlabel('Wealth Quartiles')
ylabel('Average Health Effort')
title('College')
xticks([1 2 3 4])
xticklabels({'1st', '2nd', '3rd', '4th'});
ll = legend([d1 d2 d3],'Age: 25-44', 'Age: 45-64', 'Age: 65-84');
set(ll,'Location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,'figures\Figure3','epsc2')


%% Figure 4
hh4 = figure;
subplot(1,2,1)
d1 = plot(1:38,survival_rates(:,2),'r','LineWidth',2);
hold on
d2 = plot(1:38,survival_rates(:,1),'LineStyle','-.','Color',[0 0.5 0],'LineWidth',2);
xlabel('Age')
ylabel('Survival Prob.')
xticks([1 6 11 16 21 26 31 36])
xticklabels({'25', '35', '45', '55',...
    '65', '75',	'85', '95'});
xlim([1 38]);
grid on
ll = legend([d2 d1],'Healthy', 'Unealthy');
set(ll,'Location','southwest');
title('Non-College')
subplot(1,2,2)
d1 = plot(1:38,survival_rates(:,4),'r','LineWidth',2);
hold on
d2 = plot(1:38,survival_rates(:,3),'LineStyle','-.','Color',[0 0.5 0],'LineWidth',2);
xlabel('Age')
ylabel('Survival Prob.')
xticks([1 6 11 16 21 26 31 36])
xticklabels({'25', '35', '45', '55',...
    '65', '75',	'85', '95'});
xlim([1 38]);
grid on
title('College')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,'figures\Figure4','epsc2')


%% Figure 5
hh5 = figure;
subplot(1,2,1)
dd1 = plot(1:4,frac_work10(:,1),':rx','LineWidth',2,'Marker','x','MarkerSize',10);
hold on
dd2 = plot(1:4,frac_work10(:,2),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10);
dd3 = plot(1:4,frac_work10_mod(:,1),'-.r','LineWidth',2,'Marker','s');     % Data
dd4 = plot(1:4,frac_work10_mod(:,2),'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s');        % Data
ylim([0 1]);
grid on
xlabel('Age Group')
ylabel('Share')
xticks([1 2 3 4 ])
xticklabels({'25-34', '35-44', '45-54', '55-64'});
xlim([1 4]);
ll = legend([dd4 dd3 dd2 dd1],'Model: Healthy', 'Model: Unhealthy', 'Data: Healthy', 'Data: Unhealthy');
set(ll,'Location','southwest');
title('Employment by Health')
subplot(1,2,2)
xx1 = plot(1:4,frac_work_edu10(:,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':','LineWidth',2,'Marker','x','MarkerSize',10);
hold on
xx2 = plot(1:4,frac_work_edu10(:,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','-.','Marker','x','LineWidth',2,'MarkerSize',10);
xx3 = plot(1:4,frac_work_edu10_mod(:,1),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2,'Marker','s');     % Data
xx4 = plot(1:4,frac_work_edu10_mod(:,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','LineWidth',2,'Marker','s');        % Data
xlabel('Age Group')
ylabel('Share')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'25-34', '35-44', '45-54', '55-64'});
title('Employment by Education')
ll = legend([xx3 xx4 xx1 xx2], 'Model: Non-CL', 'Model: CL','Data: Non-CL', 'Data: CL');
set(ll,'Location','southwest');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[100 100 1200 500]);
saveas(gcf,'figures\Figure5','epsc2')


%% Figure 6
hh6 = figure;
subplot(1,2,1)
dd1 = plot(1:4,avrg_inc10_edu_health(:,1)./(1000),':rx','LineWidth',2,'Marker','x','MarkerSize',10)
hold on
dd2 = plot(1:4,avrg_inc10_edu_health(:,2)./(1000),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10)
dd3 = plot(1:4,avrg_inc10_edu_health_mod(:,1),'-.r','LineWidth',2,'Marker','s')
dd4 = plot(1:4,avrg_inc10_edu_health_mod(:,2),'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s')%      % Data
hold off
ylim([20 130])
grid on
xlabel('Age Group')
ylabel('Labor Income in 2015 EUR (Ths.)')
xticks([1 2 3 4 ])
xticklabels({'25-34', '35-44', '45-54', '55-64'});
xlim([1 4]);
ll = legend([dd4 dd3 dd2 dd1],'Model: Healthy', 'Model: Unhealthy', 'Data: Healthy', 'Data: Unhealthy');
set(ll,'Location','northwest');
title('Non-College')
subplot(1,2,2)
xx1 = plot(1:4,avrg_inc10_edu_health(:,3)./(1000),':rx','LineWidth',2,'Marker','x','MarkerSize',10)
hold on
xx2 = plot(1:4,avrg_inc10_edu_health(:,4)./(1000),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10)
xx3 = plot(1:4,avrg_inc10_edu_health_mod(:,3),'-.r','LineWidth',2,'Marker','s')
xx4 = plot(1:4,avrg_inc10_edu_health_mod(:,4),'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s')
ylim([20 130])
xlabel('Age Group')
ylabel('Labor Income in 2015 EUR (Ths.)')
grid on
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'25-34', '35-44', '45-54', '55-64'});
title('College')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[100 100 1200 500]);
saveas(gcf,'figures\Figure6','epsc2')


%% Figure 7
hh7 = figure;
subplot(1,2,1)
xx1 = plot(1:6,eff10_edu_health(1:6,1),'Color',[1 0 0],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',10);
hold on
xx2 = plot(1:6,eff10_edu_health(1:6,3),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10);
xx3 = plot(1:6,eff10_edu_health_mod(1:6,1),'Color',[1 0 0],'LineStyle','-.','Marker','s','LineWidth',2);
xx4 = plot(1:6,eff10_edu_health_mod(1:6,3),'Color',[0 0.5 0],'LineStyle','-','Marker','s','LineWidth',2);
grid on
ylim([0.55 0.85])
xlabel('Age Group')
ylabel('Average Health Effort')
title('Non-College')
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
ll = legend([xx4 xx3 xx2 xx1],'Model: Healthy', 'Model: Unhealthy', 'Data: Healthy', 'Data: Unhealthy');
set(ll,'Location','northwest');

subplot(1,2,2)
xx1 = plot(1:6,eff10_edu_health(1:6,2),'Color',[1 0 0],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',10);
hold on
xx2 = plot(1:6,eff10_edu_health(1:6,4),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10);
xx3 = plot(1:6,eff10_edu_health_mod(1:6,2),'Color',[1 0 0],'LineStyle','-.','Marker','s','LineWidth',2);
xx4 = plot(1:6,eff10_edu_health_mod(1:6,4),'Color',[0 0.5 0],'LineStyle','-','Marker','s','LineWidth',2);
grid on
ylim([0.55 0.85])
xlabel('Age Group')
ylabel('Average Health Effort')
title('College')
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\Figure7','epsc2')


%% Figure 8
hh8 = figure;
subplot(1,2,1)
dd2 = semilogy(1:6,wealth10_dat(1:6,2).*1000,'Color',[0 0 0],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',10);
hold on
dd5 = semilogy(1:6,wealth10_mod(1:6,2).*1000,'Color',[0 0.4470 0.7410],'LineStyle','-','Marker','s','LineWidth',2);
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34','35-44', '45-54', '55-64', '65-74', '75-84'});
ll = legend([dd5 dd2],'Model','Data');
set(ll,'Location','northwest');
title('Median Wealth');

subplot(1,2,2)
xx1 = plot(1:6,health10(1:6,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',10);
hold on
xx2 = plot(1:6,health10(1:6,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10);
xx3 = plot(1:6,health10_mod(1:6,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(1:6,health10_mod(1:6,2),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel('Share Unhealthy')
grid on
ylim([0 0.5]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
title('Health Evolution by Education');
ll = legend([xx3 xx4 xx1 xx2],'Model: Non-CL', 'Model: CL','Data: Non-CL', 'Data: CL');
set(ll,'Location','northwest');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[100 100 1200 500]);
saveas(gcf,'figures\Figure8','epsc2')


%% Figure 9
tiledlayout(1,2)
nexttile
xx1 = semilogy(1:5,wealth10_50(2:6,1).*1000,':rx','LineWidth',2,'MarkerSize',10);
hold on
xx2 = semilogy(1:5,wealth10_50(2:6,2).*1000,'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10);
xx3 = semilogy(1:5,wealth10_50_mod(2:6,1).*1000,'-.r','LineWidth',3,'Marker','s');
xx4 = semilogy(1:5,wealth10_50_mod(2:6,2).*1000,'Color',[0 0.5 0],'LineStyle','-','LineWidth',3,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'35-44', '45-54', '55-64', '65-74', '75-84'});
ll = legend([xx4 xx3 xx2 xx1 ],'Model: Healthy', 'Model: Unhealthy','Data: Healthy', 'Data: Unhealthy');
set(ll,'Location','southeast');
title('Baseline Model vs. Data');

% Model with exogenous health
nexttile
xx1 = semilogy(1:5,wealth10_50(2:6,1).*1000,':rx','LineWidth',2,'MarkerSize',10);
hold on
xx2 = semilogy(1:5,wealth10_50(2:6,2).*1000,'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',10);
xx3 = semilogy(1:5,wealth10_50_mod2(2:6,1).*1000,'-.r','LineWidth',3,'Marker','s');
xx4 = semilogy(1:5,wealth10_50_mod2(2:6,2).*1000,'Color',[0 0.5 0],'LineStyle','-','LineWidth',3,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'35-44', '45-54', '55-64', '65-74', '75-84'});
title('Exogenous Health Model vs. Data');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1100 500]);
saveas(gcf,'figures\Figure9','epsc2')


%% Figure 10
xaxis = 1:4;
e = 5;
hh10 = figure;
subplot(1,3,1)
xx2 = plot(xaxis,whgap25_per(2:e,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap25_per(2:e,2),'Color',[0.6350 0.0780 0.1840],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(xaxis,whgap25_per(2:e,3),'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',2,'Marker','s');
ylabel(' Wealth Gap rel. to Healthy ')
xlabel('Age Group')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('25th Prctl. Wealth');

subplot(1,3,2)
xx2 = plot(xaxis,whgap50_per(2:e,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap50_per(2:e,2),'Color',[0.6350 0.0780 0.1840],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(xaxis,whgap50_per(2:e,3),'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('50th Prctl. Wealth');

subplot(1,3,3)
xx2 = plot(xaxis,whgap75_per(2:e,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap75_per(2:e,2),'Color',[0.6350 0.0780 0.1840],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(xaxis,whgap75_per(2:e,3),'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('75th Prctl. Wealth');
ll = legend([xx2 xx3 xx4], 'Baseline Model', 'Earnings Channel', 'Savings Channel');
set(ll,'Location','northeast');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\Figure10','epsc2')


%% Figure 11
hh11 = figure;
subplot(1,3,1)
xx2 = plot(xaxis,whgap25_per(2:e,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap25_per(2:e,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','-.','LineWidth',2,'Marker','s');
ylabel(' Wealth Gap rel. to Healthy ')
xlabel('Age Group')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
ll = legend([xx2 xx3], 'Baseline Model', 'Equal Efforts');
set(ll,'Location','southwest');
title('25th Prctl. Wealth');

subplot(1,3,2)
xx2 = plot(xaxis,whgap50_per(2:e,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap50_per(2:e,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','-.','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('50th Prctl. Wealth');

subplot(1,3,3)
xx2 = plot(xaxis,whgap75_per(2:e,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap75_per(2:e,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','-.','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('75th Prctl. Wealth');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\Figure11','epsc2')


%% Figure A.3
hhA3 = figure;
subplot(3,1,1)
xx1 = plot(1:6,smoking(1:6,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,smoking(1:6,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Index')
grid on
ylim([-0.6 0.6]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75+'});
title('Smoking');
ll = legend([xx1 xx2],'Non-College', 'College');
set(ll,'Location','southeast');
subplot(3,1,2)
xx1 = plot(1:6,sport(1:6,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,sport(1:6,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Index')
grid on
ylim([-0.6 0.6]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75+'});
title('Sport and Exercise');
subplot(3,1,3)
xx1 = plot(1:6,nutrition(1:6,1),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':','Marker','x','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,nutrition(1:6,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Index')
grid on
ylim([-0.6 0.6]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75+'});
title('Healthy Nutrition');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[10 10 700 700]);
saveas(gcf,'figures\FigureA3','epsc2')


%% Figure A.4
hhA4 = figure;
subplot(1,3,1)
xx2 = plot(xaxis,whgap25_per(2:5,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx4 = plot(xaxis,whgap25_per(2:5,8),'Color',[0.4660 0.6740 0.1880],'LineStyle','-.','LineWidth',2,'Marker','d');
ylabel(' Wealth Gap rel. to Healthy ')
xlabel('Age Group')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('25th Prctl. Wealth');

subplot(1,3,2)
xx2 = plot(xaxis,whgap50_per(2:5,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx4 = plot(xaxis,whgap50_per(2:5,8),'Color',[0.4660 0.6740 0.1880],'LineStyle','-.','LineWidth',2,'Marker','d');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('50th Prctl. Wealth');

subplot(1,3,3)
xx2 = plot(xaxis,whgap75_per(2:5,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx4 = plot(xaxis,whgap75_per(2:5,8),'Color',[0.4660 0.6740 0.1880],'LineStyle','-.','LineWidth',2,'Marker','d');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
ll = legend([xx2 xx4], 'Baseline Model', 'All Channels');
set(ll,'Location','northwest');
title('75th Prctl. Wealth');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\FigureA4','epsc2')


%% Figure A.5
hhA5 = figure;
subplot(1,3,1)
xx2 = plot(xaxis,whgap50_per(2:5,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap50_per(2:5,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(xaxis,whgap50_per(2:5,5),'Color',[0.8500 0.3250 0.0980],'LineStyle','-.','LineWidth',2,'Marker','d');
ylabel(' Wealth Gap rel. to Healthy ')
xlabel('Age Group')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
ll = legend([xx2 xx3 xx4], 'Baseline Model', 'Equal Effort', 'Eq. Effort: Age 25-44');
set(ll,'Location','southwest');
title('50th Prctl. Wealth');

subplot(1,3,2)
xx2 = plot(xaxis,whgap50_per(2:5,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap50_per(2:5,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(xaxis,whgap50_per(2:5,6),'Color',[0.8500 0.3250 0.0980],'LineStyle','-.','LineWidth',2,'Marker','d');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
ll = legend([xx4], 'Eq. Effort: Age 45-64');
set(ll,'Location','southwest');
title('50th Prctl. Wealth');

subplot(1,3,3)
xx2 = plot(xaxis,whgap50_per(2:5,1),'Color',[0 0.4470 0.7410],'LineStyle','--','Marker','s','LineWidth',2);
hold on
xx3 = plot(xaxis,whgap50_per(2:5,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','-.','LineWidth',2,'Marker','s');
xx4 = plot(xaxis,whgap50_per(2:5,7),'Color',[0.8500 0.3250 0.0980],'LineStyle','-.','LineWidth',2,'Marker','d');
xlabel('Age Group')
ylabel(' Wealth Gap rel. to Healthy ')
grid on
ylim([0 1]);
xlim([1 4]);
xticks([1 2 3 4])
xticklabels({'35-44', '45-54', '55-64', '65-74'});
title('50th Prctl. Wealth');
ll = legend([xx4], 'Eq. Effort: Age 65+');
set(ll,'Location','southwest');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\FigureA5','epsc2')


%% Figure A.6
hhA6 = figure;
subplot(1,2,1)
xx1 = semilogy(1:5,wealth10_50_occ(1:5,3).*1000,':rx','LineWidth',2,'MarkerSize',12);
hold on
xx2 = semilogy(1:5,wealth10_50_occ(1:5,4).*1000,'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74'});
title('Manual Occupations');
subplot(1,2,2)
xx1 = semilogy(1:5,wealth10_50_occ(1:5,1).*1000,':rx','LineWidth',2,'MarkerSize',12);    
hold on
xx2 = semilogy(1:5,wealth10_50_occ(1:5,2).*1000,'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74'});
title('Non-Manual Occupations');
ll = legend([xx2 xx1],'Healthy', 'Unhealthy');
set(ll,'Location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\FigureA6','epsc2')


%% Figure A.7
hhA7 = figure;
subplot(1,3,1)
xx1 = plot(1:6,wealth10_25(1:6,1),':rx','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,wealth10_25(1:6,2),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xx3 = plot(1:6,wealth10_25_mod(1:6,1),'-.r','LineWidth',2,'Marker','s');
xx4 = plot(1:6,wealth10_25_mod(1:6,2),'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([0 220]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
ll = legend([xx4 xx3 xx2 xx1],'Model: Healthy', 'Model: Unhealthy', 'Data: Healthy','Data: Unhealthy' );
set(ll,'Location','northwest');
title('25th Prctl. Wealth');

subplot(1,3,2)
xx1 = plot(1:6,wealth10_50(1:6,1),':rx','LineWidth',2,'MarkerSize',12);      % Maybe x2
hold on
xx2 = plot(1:6,wealth10_50(1:6,2),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xx3 = plot(1:6,wealth10_50_mod(1:6,1),'-.r','LineWidth',2,'Marker','s');       % flipped bad with good health
xx4 = plot(1:6,wealth10_50_mod(1:6,2),'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([0 220]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
title('50th Prctl. Wealth');

subplot(1,3,3)
xx1 = plot(1:6,wealth10_75(1:6,1),':rx','LineWidth',2,'MarkerSize',12);
hold on
xx2 = plot(1:6,wealth10_75(1:6,2),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xx3 = plot(1:6,wealth10_75_mod(1:6,1),'-.r','LineWidth',2,'Marker','s');
xx4 = plot(1:6,wealth10_75_mod(1:6,2),'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([0 220]);
xlim([1 6]);
xticks([1 2 3 4 5 6])
xticklabels({'25-34', '35-44', '45-54', '55-64', '65-74', '75-84'});
title('75th Prctl. Wealth');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(gcf,'Position',[50 50 1300 600]);
saveas(gcf,'figures\FigureA7','epsc2')


%% Figure A.8
hhA8 = figure;
subplot(1,2,1)
xx1 = semilogy(1:5,wealth10_50_edu(2:6,1),':rx','LineWidth',2,'MarkerSize',12);      % Maybe x2
hold on
xx2 = semilogy(1:5,wealth10_50_edu(2:6,2),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xx3 = semilogy(1:5,wealth10_50_edu_mod(2:6,1).*1000,'-r','LineWidth',2,'Marker','s');       
xx4 = semilogy(1:5,wealth10_50_edu_mod(2:6,2).*1000,'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'35-44', '45-54', '55-64', '65-74', '75-84'});
title('Non-College');
subplot(1,2,2)
xx1 = semilogy(1:5,wealth10_50_edu(2:6,3),':rx','LineWidth',2,'MarkerSize',12);      % Maybe x2
hold on
xx2 = semilogy(1:5,wealth10_50_edu(2:6,4),'Color',[0 0.5 0],'LineStyle','--','Marker','x','LineWidth',2,'MarkerSize',12);
xx3 = semilogy(1:5,wealth10_50_edu_mod(2:6,3).*1000,'-r','LineWidth',2,'Marker','s');       
xx4 = semilogy(1:5,wealth10_50_edu_mod(2:6,4).*1000,'Color',[0 0.5 0],'LineStyle','-','LineWidth',2,'Marker','s');
xlabel('Age Group')
ylabel('Wealth in 2015 EUR (Ths.)')
grid on
ylim([1250 200000])
yticks([2500 5000 10000 20000 40000 80000 160000])
yticklabels({'2.5','5','10','20','40','80','160',})
xlim([1 5]);
xticks([1 2 3 4 5])
xticklabels({'35-44', '45-54', '55-64', '65-74', '75-84'});
title('College');
ll = legend([xx4 xx3 xx2 xx1],'Model: Healthy', 'Model: Unhealthy', 'Data: Healthy', 'Data: Unhealthy');
set(ll,'Location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(gcf,'Position',[50 50 1200 500]);
saveas(gcf,'figures\FigureA8','epsc2')
