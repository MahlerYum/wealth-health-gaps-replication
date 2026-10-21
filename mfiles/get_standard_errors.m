%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Standard Errors of Estimates in Mahler and Yum, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Load Empirical moments and standard errors
cd('..\workdata')

%% From Stata Output
% Vector of noisy data moment estimates and their estimated standard errors
mu_hat = importdata('mu_hat.txt');
moment_se = importdata('moment_se.txt');

%% From Fortran Output (or Jacobian from Excel calculation)
% Vector of initial parameter estimates from minimization using
% efficient weighting matrix 
theta_init = importdata('theta_init.txt');

% Vector of model-implied moments resulting from first estimation
h_theta = importdata('h_theta.txt');

% Jacobian matrix coming from forward perturbation of each parameter
G_fwd = importdata('G_fwd.txt');


%% Compute one-step estimators and worst-case standard errors using methodolgy in Plagborg-Moller & Cocci (2021) 
% Also use part of their code

cd('..\mfiles')

moment_jacob = G_fwd';
[p,k]  = size(moment_jacob);
GpG = moment_jacob'*moment_jacob;
param_ses = zeros(k,1);
est_onestep= zeros(k,1);
moment_loadings = zeros(p,k);

for t=1:k
transf_jacob = zeros(1,k);
transf_jacob(1,t) = 1;

% Set up median regression as described in paper
    Y = moment_jacob*(GpG\(transf_jacob'));
    moment_jacob_perp = null(moment_jacob');
    X = -moment_jacob_perp;
    
    
    % Run median regression
        z = qreg(moment_se.*Y, moment_se.*X, 0.5);
        moment_loadings(:,t) = Y-X*z;
        param_ses(t,1) = moment_se'*abs(moment_loadings(:,t));
        
     % Compute 1-step efficient estimator for that parameter
     est_onestep(t,1) = theta_init(t,1) + moment_loadings(:,t)'*(mu_hat - h_theta);
     
end

%% Create Table 2
tab2 = zeros(k,2);
tab2(1:9,:) = [theta_init(1:9),param_ses(1:9)];
tab2(10,:) = [theta_init(41),param_ses(41)];
tab2(11:12,:) = [theta_init(35:36),param_ses(35:36)];
tab2(13,:) = [theta_init(42),param_ses(42)];
tab2(14:15,:) = [theta_init(37:38),param_ses(37:38)];
tab2(16:17,:) = [theta_init(39:40),param_ses(39:40)];
tab2(18:33,:) = [theta_init(10:25),param_ses(10:25)];
tab2(34,:) = [theta_init(34),param_ses(34)];
tab2(35:36,:) = [theta_init(26:27),param_ses(26:27)];
tab2(35:36,:) = [theta_init(26:27),param_ses(26:27)];
tab2(37,:) = [theta_init(33),param_ses(33)];
tab2(38:39,:) = [theta_init(30:31),param_ses(30:31)];
tab2(40:41,:) = [theta_init(28:29),param_ses(28:29)];
tab2(42,:) = [theta_init(32),param_ses(32)];

colnames = {'Estimate' 'Standard Error' };
Tab2  = array2table(tab2, 'VariableNames', colnames);
table2latex(Tab2,'..\tables\Table2')
