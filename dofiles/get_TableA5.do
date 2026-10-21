**********************************************************
******* Calculate Empirical Data Moments along with their Standard Errors to create Table A.5 (in slightly different order)
**********************************************************

cd "$path\workdata"
use workdata04-18.dta, clear

* Create matrix to store all empirical moments required for the MSM estimation
matrix moments = J(64,2,.)		
matrix colnames moments = Value SE
matrix rownames moments = "Empl_Hlty_25" "Empl_Hlty_35" "Empl_Hlty_45" "Empl_Hlty_55" "Empl_UHlty_25" "Empl_UHlty_35" "Empl_UHlty_45" "Empl_UHlty_55" "Effort_NCl_Hlty_25" "Effort_NCl_Hlty_35" "Effort_NCl_Hlty_45" "Effort_NCl_Hlty_55" "Effort_NCl_Hlty_65" "Effort_NCl_Hlty_75" "Effort_NCl_UHlty_25" "Effort_NCl_UHlty_35" "Effort_NCl_UHlty_45" "Effort_NCl_UHlty_55" "Effort_NCl_UHlty_65" "Effort_NCl_UHlty_75" "Effort_Cl_Hlty_25" "Effort_Cl_Hlty_35" "Effort_Cl_Hlty_45" "Effort_Cl_Hlty_55" "Effort_Cl_Hlty_65" "Effort_Cl_Hlty_75" "Effort_Cl_UHlty_25" "Effort_Cl_UHlty_35" "Effort_Cl_UHlty_45" "Effort_Cl_UHlty_55" "Effort_Cl_UHlty_65" "Effort_Cl_UHlty_75" "Med_Wealth_over_LabInc_25" "Med_Wealth_over_LabInc_35" "Med_Wealth_over_LabInc_45" "Med_Wealth_over_LabInc_55" "Med_Wealth_over_LabInc_65" "Med_Wealth_over_LabInc_75" "Educ_Gradient_in_Empl" "Non-Adj_Share_25"  "Non-Adj_Share_45"  "Non-Adj_Share_65" "VSL_multiple" " StdDev_Effort" "Cons_Ratio" "LabInc_NCl_Hlty_25" "LabInc_NCl_Hlty_35" "LabInc_NCl_Hlty_45" "LabInc_NCl_Hlty_55" "LabInc_NCl_UHlty_25" "LabInc_NCl_UHlty_35" "LabInc_NCl_UHlty_45" "LabInc_NCl_UHlty_55" "LabInc_Cl_Hlty_25" "LabInc_Cl_Hlty_35" "LabInc_Cl_Hlty_45" "LabInc_Cl_Hlty_55" "LabInc_Cl_UHlty_25" "LabInc_Cl_UHlty_35" "LabInc_Cl_UHlty_45" "LabInc_Cl_UHlty_55" "Var(Log_LabInc)" "Pension_Repl_Rate" "Wealth_Gini"


local row = 1
local reps_bs = 500		// no. of bootstrap repetitions
local rseed = 12345		// seed for random number generation


**********************************************************

*** Employment by 10-year age group and health category
tabstat working2y if age_cat10 <65 & health_cat == 1, by(age_cat10) stats(mean semean) save
tabstatmat workstat_h
forvalues r = 1/4{
mat moments[`row',1]=workstat_h[`r',1]
mat moments[`row',2]=workstat_h[`r',2]
local row = `row'+1
}
tabstat working2y if age_cat10 <65 & health_cat == 0, by(age_cat10) stats(mean semean) save
tabstatmat workstat_uh
forvalues r = 1/4{
mat moments[`row',1]=workstat_uh[`r',1]
mat moments[`row',2]=workstat_uh[`r',2]
local row = `row'+1
}

*** Average Health Effort by 10-year age group, college, and health category
tabstat eff_norm if college == 0 & health_cat == 1, by(age_cat10) stats(mean semean) save
tabstatmat effstat_ncl_h
forvalues r = 1/6{
mat moments[`row',1]=effstat_ncl_h[`r',1]
mat moments[`row',2]=effstat_ncl_h[`r',2]
local row = `row'+1
}
tabstat eff_norm if college == 0 & health_cat == 0, by(age_cat10) stats(mean semean) save
tabstatmat effstat_ncl_uh
forvalues r = 1/6{
mat moments[`row',1]=effstat_ncl_uh[`r',1]
mat moments[`row',2]=effstat_ncl_uh[`r',2]
local row = `row'+1
}
tabstat eff_norm if college == 1 & health_cat == 1, by(age_cat10) stats(mean semean) save
tabstatmat effstat_cl_h
forvalues r = 1/6{
mat moments[`row',1]=effstat_cl_h[`r',1]
mat moments[`row',2]=effstat_cl_h[`r',2]
local row = `row'+1
}
tabstat eff_norm if college == 1 & health_cat == 0, by(age_cat10) stats(mean semean) save
tabstatmat effstat_cl_uh
forvalues r = 1/6{
mat moments[`row',1]=effstat_cl_uh[`r',1]
mat moments[`row',2]=effstat_cl_uh[`r',2]
local row = `row'+1
}

*** Median wealth by age group
forvalues ac = 25(10)75{
	bootstrap p50=r(p50), reps(`reps_bs') seed(`rseed'): summarize net_wealth0_norm if age_cat10 == `ac', d
	mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1
}

*** Education Gradient in Employment
capture program drop myemplratio
program myemplratio, rclass
tempname mean1
qui sum working2y if college == 1 & age_cat10<65
scalar `mean1' = r(mean)
qui sum working2y if college == 0 & age_cat10<65
return scalar ratio = `mean1'/r(mean)
end
// Bootstrap Standard Error of this Ratio
bootstrap ratio=r(ratio), reps(`reps_bs') seed(`rseed'): myemplratio
mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1

*** Share of non-adjusters in effort by long age group
tsset pid syear
gen l2eff_norm = l2.eff_norm
cap program drop myinaction
program myinaction, rclass
cap drop delta abs_delta inac
gen delta = eff_norm - l2eff_norm
gen abs_delta = abs(delta)
gen inac = 1 if abs_delta ==0
replace inac = 0 if abs_delta > 0 & abs_delta !=.
qui summarize inac if long_age_cat == 1, meanonly
return scalar mean1 = r(mean)
qui summarize inac if long_age_cat == 2, meanonly
return scalar mean2= r(mean)
qui summarize inac if long_age_cat == 3, meanonly
return scalar mean3 = r(mean)
end

bootstrap mean1=r(mean1) mean2=r(mean2) mean3=r(mean3), reps(`reps_bs') seed(`rseed'): myinaction
forvalues cc=1/3{
mat moments[`row',1]=r(table)[1,`cc']
mat moments[`row',2]=r(table)[2,`cc']
local row = `row'+1
}

*** Value of a statistical life year multiple of consumption
/* This value is calculated starting from VSL estimates for EU countries as reported in Table 6.1 in OECD (2012). The procedure is as follows (see also the explanation on page 20 in the paper):

1. Take the average VSL estimate as reported in table 6.1 of OECD (2012) and transform into 2018 Euros using a CPI of 1.29 and an exchange rate between $ and € of 0.8475 (Source: https://www.usinflationcalculator.com/ and https://www.exchangerates.org.uk/USD-EUR-spot-exchange-rates-history-2018.html)

2. Transform this VSL into a value of a statistical life year of an average person by using the fact that the average life expctancy of an average German man in 2018 (who was 44.4 years old) was 34.8 years and an assumed annual discount factor of 1\% (Source: https://www.statistischebibliothek.de/mir/servlets/MCRFileNodeServlet/DEHeft_derivate_00056831/5126204187004a.pdf)

3. Divide this VSLY by annual average per-capita consumption in Germany in 2018 (20420 EUR according to https://www.destatis.de/EN/Press/2020/03/PE20_089_811.html)
*/
scalar vsl_05 = 4704038
scalar vsl_05_se = 329474

scalar vsl_18_EUR = vsl_05 * 1.29 * 0.8475
scalar vsl_18_EUR_se = vsl_05_se * 1.29 *0.8475

scalar le_disc = 0
forvalues yy=1/35{
	scalar le_disc = le_disc + 0.99^(`yy'-1)
}

scalar vsly = vsl_18_EUR / le_disc
scalar vsly_se = vsl_18_EUR_se / le_disc

mat moments[`row',1]=vsly/20420
mat moments[`row',2]=vsly_se/20420
local row = `row'+1

*** Standard deviation of effort
bootstrap sd=r(sd), reps(`reps_bs') seed(`rseed'): summarize eff_norm if age <85
mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1

*** Ratio of average consumption of healthy to unhealthy
* Since consumption is reported at the HH level but our model is at the individual level, we use one-person households to calculate the share in annual consumption of a healthy versus an unhealthy individual
capture program drop myconsratio
program myconsratio, rclass
tempname mean1
qui sum cons_exp_annl if d11106 == 1 & health_cat == 1
scalar `mean1' = r(mean)
qui sum cons_exp_annl if d11106 == 1 & health_cat == 0
return scalar ratio = `mean1'/r(mean)
end
bootstrap ratio=r(ratio), reps(`reps_bs') seed(`rseed'): myconsratio
mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1

*** Average labor income by 10-year age cat, health, and education (including zeros) in Thousands
tabstat labinc2y if age_cat10 <65 & health_cat == 1 & college == 0, by(age_cat10) stats(mean semean) save
tabstatmat labinstat_h_ncl
forvalues r = 1/4{
mat moments[`row',1]=labinstat_h_ncl[`r',1]/1000
mat moments[`row',2]=labinstat_h_ncl[`r',2]/1000
local row = `row'+1
}
tabstat labinc2y if age_cat10 <65 & health_cat == 0 & college == 0, by(age_cat10) stats(mean semean) save
tabstatmat labinstat_uh_ncl
forvalues r = 1/4{
mat moments[`row',1]=labinstat_uh_ncl[`r',1]/1000
mat moments[`row',2]=labinstat_uh_ncl[`r',2]/1000
local row = `row'+1
}
tabstat labinc2y if age_cat10 <65 & health_cat == 1 & college == 1, by(age_cat10) stats(mean semean) save
tabstatmat labinstat_h_cl
forvalues r = 1/4{
mat moments[`row',1]=labinstat_h_cl[`r',1]/1000
mat moments[`row',2]=labinstat_h_cl[`r',2]/1000
local row = `row'+1
}
tabstat labinc2y if age_cat10 <65 & health_cat == 0 & college == 1, by(age_cat10) stats(mean semean) save
tabstatmat labinstat_uh_cl
forvalues r = 1/4{
mat moments[`row',1]=labinstat_uh_cl[`r',1]/1000
mat moments[`row',2]=labinstat_uh_cl[`r',2]/1000
local row = `row'+1
}

*** Variance of log labor income (cond. on working) in working age population by education
bootstrap var=r(Var), reps(`reps_bs') seed(`rseed'): summarize loglabinc2y_workers
mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1

*** Ratio of Pension Benefits to average Labor Income
capture program drop mypensratio
program mypensratio, rclass
tempname mean1
qui sum labinc if age>25 & age<65
scalar `mean1' = r(mean)
qui sum pension_ind 
return scalar ratio = r(mean)/`mean1'
end
bootstrap ratio=r(ratio), reps(`reps_bs') seed(`rseed'): mypensratio
mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1

*** Gini coefficient of wealth
capture program drop mygini
program mygini, rclass
qui inequal7 net_wealth0 if age >=25, returnscalars
return scalar gini = r(gini)
end
bootstrap gini=r(gini), reps(`reps_bs') seed(`rseed'): mygini
mat moments[`row',1]=r(table)[1,1]
mat moments[`row',2]=r(table)[2,1]
local row = `row'+1



************************************************
*** Export Moment Table A5


matrix list moments
outtable using "$path\tables\TableA5", mat(moments) replace format (%9.3f) caption("Empirical Moments and Standard Errors") nobox center