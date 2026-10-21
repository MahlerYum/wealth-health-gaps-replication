%%% Estimating Productivity Process Parameters (AR(1) coefficient)
%%% Following STY (2004, JME)
clear all

cd('..\workdata')
% Import combined residuals from fixed effects wage regressions (after
% residualizing against cohort and education)
RegResid = importdata('prod_data.csv');
cd('..\mfiles')
Data = RegResid.data;
Data(isnan(Data))=0;
na = size(Data,2);

ECov = zeros(na,na);
for a = 1:na
    for b = 1:na
        x = Data(:,a).*[Data(:,a) ~= 0];
        xn = Data(:,b).*[Data(:,b) ~= 0];
        ECov(a,b) = (x'*xn) / sum(x~=0); 
    end
end

EVar = diag(ECov);

% If you want to include temporary shock, need to identify this separately
% following process in SYT (2004)

% Calculate the sum of residuals over three years
SResid = zeros(length(Data),na-2);
for j = 1:na-2
    for ii = 1:length(Data)
        if (Data(ii,j) ~=0 & Data(ii,j+1) ~= 0 & Data(ii,j+2)~=0)
            SResid(ii,j) = Data(ii,j) + Data(ii,j+1) + Data(ii,j+2);
        else
            SResid(ii,j) =0;
        end
    end
end
% Calculate the variance of these sums for each age group
SResid(SResid ==0) = NaN;
SVar = var(SResid,0,'omitnan')';

% SVar now includes na-2 additional moments that can be included in Mhat
% (Pick the first one (i.e. at age 25) but include the whole variance
% vector to identify rho)
Mhat = [EVar(1:na) ; SVar(1)];

fun = @(y) AR1GMM(y(1),y(2),y(3),y(4),Mhat,na);
y_0 = [0.0 0.0 0.9 0.0];
param = fminsearch(fun,y_0,optimset('TolX',1e-7,'TolFun',1e-7,'Display','iter','FunValCheck','on','MaxFunEvals',10000,'MaxIter',2000));

sigma2_perm = param(1);
sigma2_pers = param(2);
rho = param(3);
sigma2_trans = param(4);

% Save estimated AR(1) coefficient in text file
T = table(rho,'VariableNames',{'AR(1) Coefficient'});
cd('..\tables')
writetable(T,'AR(1)Coefficient.txt')


