%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Tables in Mahler and Yum, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Load data used to create figures
cd('..\workdata')
load data_fig.mat

%% Load variables (empirical data from stata
eff_adj = readmatrix('data_tables.xlsx','Sheet','eff_adj');


%% Load data from model simulation
cd('..')
cd('fortran\benchmark')

eff_adj_mod = importdata('Table3model.txt');

dat = importdata('Table4model1.txt');
eff_weq = dat(1:16);
he_weq = dat(17:32);

dat_exp = importdata('Table4model2.txt');
eff_weq_exp = dat_exp(1:16);
he_weq_exp = dat_exp(17:32);

whgap25_per(2:5,9) = dat_exp(33:36);
whgap50_per(2:5,9) = dat_exp(37:40);
whgap75_per(2:5,9) = dat_exp(41:44);

tabA11 = importdata('TableA11model.txt');

dat_exp = importdata('TableA12model1.txt')./100;
whgap25_per(2:5,10) = dat_exp(1:4);
whgap50_per(2:5,10) = dat_exp(5:8);
whgap75_per(2:5,10) = dat_exp(9:12);
dat_exp2 = importdata('TableA12model2.txt')./100;
whgap25_per(2:5,11) = dat_exp2(1:4);
whgap50_per(2:5,11) = dat_exp2(5:8);
whgap75_per(2:5,11) = dat_exp2(9:12);
dat_exp3 = importdata('TableA12model3.txt')./100;
whgap25_per(2:5,12) = dat_exp3(1:4);
whgap50_per(2:5,12) = dat_exp3(5:8);
whgap75_per(2:5,12) = dat_exp3(9:12);
dat_exp4 = importdata('TableA12model4.txt')./100;
whgap25_per(2:5,13) = dat_exp4(1:4);
whgap50_per(2:5,13) = dat_exp4(5:8);
whgap75_per(2:5,13) = dat_exp4(9:12);

dat_exp5 = importdata('TableA13model1.txt');
whgap25_per(2:5,14) = dat_exp5(1:4);
whgap50_per(2:5,14) = dat_exp5(5:8);
whgap75_per(2:5,14) = dat_exp5(9:12);

dat_exp6 = importdata('TableA13model2.txt');
whgap25_per(2:5,15) = dat_exp6(1:4);
whgap50_per(2:5,15) = dat_exp6(5:8);
whgap75_per(2:5,15) = dat_exp6(9:12);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('..')
cd('..\mfiles')

%% Table 3
tab3(1:3,1) = eff_adj_mod(1:3);
tab3(1:3,2) = eff_adj(1:3);
tab3(1:3,3) = eff_adj_mod(4:6);
tab3(1:3,4) = eff_adj(4:6);
tab3(1:3,5) = eff_adj_mod(7:9);
tab3(1:3,6) = eff_adj(7:9);
tab3(1:3,7) = eff_adj_mod(10:12);
tab3(1:3,8) = eff_adj(10:12);
tab3 = round(tab3,2);

colnames = {'Model (+10%)' 'Data (+10%)' 'Model (+20%)' 'Data (+20%)' 'Model (-10%)' 'Data (-10%)' 'Model (-20%)' 'Data (-20%)'};
Tab3  = array2table(tab3, 'VariableNames', colnames);
Rows = {'25-44';'45-64';'65-84'};
Tab3.Properties.RowNames = Rows; 
table2latex(Tab3,'..\tables\Table3');


%% Table 4
eff_weq_rel = round ( ((eff_weq_exp - eff_weq)./eff_weq ) .* 100,1);
tab4(1:4,1:4) = reshape(eff_weq_rel,4,4);
he_weq_rel = round ( ((he_weq_exp - he_weq)./he_weq ) .* 100,1);
tab4(1:4,5:8) = reshape(he_weq_rel,4,4);
tab4(1:4,9) = round ( ((whgap25_per(2:5,9) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,1);
tab4(1:4,10) = round ( ((whgap50_per(2:5,9) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,1);
tab4(1:4,11) = round ( ((whgap75_per(2:5,9) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,1);
colnames = {'Av. Effort (1st)' 'Av. Effort (2nd)' 'Av. Effort (3rd)' 'Av. Effort (4th)' ...
    'Health (1st)' 'Health (2nd)' 'Health (3rd)' 'Health (4th)' 'WH-Gap (25th)' 'WH-Gap (50th)' 'WH-Gap (75th)' };
Tab4  = array2table(tab4, 'VariableNames', colnames);
Rows = {'35-44';'45-54';'55-64';'65-74'};
Tab4.Properties.RowNames = Rows; 
table2latex(Tab4,'..\tables\Table4');


%% Table A.11
tabA11 = round(tabA11,1);
colnames = {'Contribution (in %)' };
TabA11  = array2table(tabA11, 'VariableNames', colnames);
Rows = {'Lifetime Earnings';'Wealth at 65-66';'No. of Healthy Years';'Share of Healthy Years'};
TabA11.Properties.RowNames = Rows; 
table2latex(TabA11,'..\tables\TableA11');

%% Table A.12
tab12(1:4,1) = round ( (-1).* ((whgap25_per(2:5,2) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab12(1:4,2) = round ( (-1).* ((whgap50_per(2:5,2) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab12(1:4,3) = round ( (-1).* ((whgap75_per(2:5,2) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab12(1:4,4) = round ( (-1).* ((whgap25_per(2:5,10) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab12(1:4,5) = round ( (-1).* ((whgap50_per(2:5,10) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab12(1:4,6) = round ( (-1).* ((whgap75_per(2:5,10) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab12(1:4,7) = round ( (-1).* ((whgap25_per(2:5,11) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab12(1:4,8) = round ( (-1).* ((whgap50_per(2:5,11) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab12(1:4,9) = round ( (-1).* ((whgap75_per(2:5,11) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab12(5:8,1) = round ( (-1).* ((whgap25_per(2:5,3) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab12(5:8,2) = round ( (-1).* ((whgap50_per(2:5,3) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab12(5:8,3) = round ( (-1).* ((whgap75_per(2:5,3) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab12(5:8,4) = round ( (-1).* ((whgap25_per(2:5,12) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab12(5:8,5) = round ( (-1).* ((whgap50_per(2:5,12) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab12(5:8,6) = round ( (-1).* ((whgap75_per(2:5,12) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab12(5:8,7) = round ( (-1).* ((whgap25_per(2:5,13) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab12(5:8,8) = round ( (-1).* ((whgap50_per(2:5,13) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab12(5:8,9) = round ( (-1).* ((whgap75_per(2:5,13) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

colnames = {'Total (25th)' 'Total (50th)' 'Total (75th)' 'Wage Loss/Length (25th)' 'Wage Loss/Length (50th)' 'Wage Loss/Length (75th)' ...
    'Disutil./Quality (25th)' 'Disutil./Quality  (50th)' 'Disutil./Quality (75th)' };
TabA12  = array2table(tab12, 'VariableNames', colnames);
Rows = {'Earnings Channel: 35-44';'Earnings Channel: 45-54';'Earnings Channel: 55-64';'Earnings Channel: 65-74';...
    'Savings Channel: 35-44';'Savings Channel: 45-54';'Savings Channel: 55-64';'Savings Channel: 65-74'};
TabA12.Properties.RowNames = Rows; 
table2latex(TabA12,'..\tables\TableA12');


%% Table A.13
tab13(1:4,1) = round ( (-1).* ((whgap25_per(2:5,4) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab13(1:4,2) = round ( (-1).* ((whgap50_per(2:5,4) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab13(1:4,3) = round ( (-1).* ((whgap75_per(2:5,4) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab13(1:4,4) = round ( (-1).* ((whgap25_per(2:5,14) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab13(1:4,5) = round ( (-1).* ((whgap50_per(2:5,14) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab13(1:4,6) = round ( (-1).* ((whgap75_per(2:5,14) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

tab13(1:4,7) = round ( (-1).* ((whgap25_per(2:5,15) - whgap25_per(2:5,1))./whgap25_per(2:5,1) ) .* 100 ,0);
tab13(1:4,8) = round ( (-1).* ((whgap50_per(2:5,15) - whgap50_per(2:5,1))./whgap50_per(2:5,1) ) .* 100 ,0);
tab13(1:4,9) = round ( (-1).* ((whgap75_per(2:5,15) - whgap75_per(2:5,1))./whgap75_per(2:5,1) ) .* 100 ,0);

colnames = {'Total (25th)' 'Total (50th)' 'Total (75th)' 'Direct (25th)' 'Direct (50th)' 'Direct (75th)' ...
    'Indirect (25th)' 'Indirect  (50th)' 'Indirect (75th)' };
TabA13  = array2table(tab13, 'VariableNames', colnames);
Rows = {'35-44';'45-54';'55-64';'65-74'};
TabA13.Properties.RowNames = Rows; 
table2latex(TabA13,'..\tables\TableA13');

