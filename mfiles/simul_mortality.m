%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation to gauge health technology parameters mortality implications
%%% (and construct Column 2 of Table A10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Load variables (empirical data from stata)
cd('..\workdata')
types_70 = readmatrix('data_tables.xlsx','Sheet','types_70');
types_7084 = readmatrix('data_tables.xlsx','Sheet','types_7084');
eff_dec = readmatrix('data_tables.xlsx','Sheet','eff_dec');
survival_rates = readmatrix('data_figures.xlsx','Sheet','survival_rates');

%% Perform simulation 
ns = 100000;
nt = 5;

age = zeros(ns,1);
alive = zeros(ns,nt);
health = zeros(ns,nt);
educ = zeros(ns,1);
health_type = zeros(ns,1);

total = sum(sum(types_70));
type_vec = [types_70(1,1)/total sum(types_70(1,:))/total (sum(types_70(1,:))+types_70(2,1))/total 1];    % cumulative distribution over college and health types at age 70

surv_hs = survival_rates(24:38,1:2);
surv_hs(16,1:2) = surv_hs(15,1:2);        % Assume constant survival rate for >100 year olds
surv_hs(17,1:2) = surv_hs(15,1:2);

surv_cl = survival_rates(24:38,3:4);
surv_cl(16,1:2) = surv_cl(15,1:2);        % Assume constant survival rate for >100 year olds
surv_cl(17,1:2) = surv_cl(15,1:2);

alive(:,1) = 1;

% Take initial health and health_type correlation as in the data at age
% 70-84
for is = 1:ns
    x = rand;
    y = rand;
    if (x<= type_vec(1)) 
        educ(is) = 0;
        health_type(is) =0;
        if (y<= types_7084(1))
            health(is,1) = 1;
        else
            health(is,1) = 0;
        end
    elseif (x<= type_vec(2))
        educ(is) = 0;
        health_type(is) =1;
            if (y<= types_7084(2))
            health(is) = 1;
        else
            health(is) = 0;
        end
    elseif (x<= type_vec(3))
        educ(is) = 1;
        health_type(is) =0;
        if (y<= types_7084(3))
            health(is,1) = 1;
        else
            health(is,1) = 0;
        end
    else
        educ(is) = 1;
        health_type(is) =1;
        if (y<= types_7084(4))
            health(is,1) = 1;
        else
            health(is,1) = 0;
        end
    end
    age(is,1) =  round(unifrnd(1,7));

end


% Set which effort level to use
% eff = 0.71; % average effort
% eff = 0.89; % top 10% effort
% eff = 0.59; % bottom 10% effort

for id = 1:3
        eff = eff_dec(id);
% Simulate survival using estimated health technology parameters
for is = 1:ns
    a = age(is,1);
    for t=1:nt-1
        
        if (alive(is,t) == 1) 
        
        % Get probability of being in good health next period (Parameter
        % estimates from Table A.6)
        if (a<=3) 
        prh = exp(-0.906 -1.138 + eff*0.693 + eff*0.734 + health(is,t)*2.311 ...
            + health_type(is)*0.632 + educ(is) * 0.238) / (1 + exp(-0.906 -1.138 ...
            + eff*0.693 + eff*0.734 + health(is,t)*2.311 + health_type(is)*0.632 + educ(is) * 0.238));
        else
        prh = exp(-0.906 -1.586 + eff*0.693 + eff*0.734 + health(is,t)*2.311 ...
            + health_type(is)*0.632 + educ(is) * 0.238) / (1 + exp(-0.906 -1.586 ...
            + eff*0.693 + eff*0.734 + health(is,t)*2.311 + health_type(is)*0.632 + educ(is) * 0.238));
        end
        
        % Get realization if healthy next period
        x = rand;
        if (x<= prh)
            health(is,t+1) = 1;
        else
            health(is,t+1) = 0;
        end
        
        % Get realization if alive next period
        y = rand;
        if educ(is) == 0
            if (health(is,t) == 0)
                svp = surv_hs(t+a,2);
            else
                svp = surv_hs(t+a,1);
            end
        else
            if (health(is,t) == 0)
                svp = surv_cl(t+a,2);
            else
                svp = surv_cl(t+a,1);
            end
        end
        
        if (y<=svp)
            alive(is,t+1) = 1;
        else
            alive(is,t+1) = 0;
        end
        
        else
            alive(is,t+1) = 0;
        end
            
    end
end

mortal_rate(id) = 1-mean(alive(:,end));

end
mortal_rate = mortal_rate';

cd('..\mfiles')
colnames = {'Mortality Rate (10y) implied by Simulation' };
TabA10  = table(mortal_rate, 'VariableNames', colnames);
Rows = {'Average Lifestyle';'Healthy Lifestyle';'Unhealthy Lifestyle'};
TabA10.Properties.RowNames = Rows; 
table2latex(TabA10,'..\tables\TableA10');