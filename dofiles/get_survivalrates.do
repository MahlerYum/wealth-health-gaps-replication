
****** Estimate survival rates

cd "$path/workdata/"
use survival_dat.dta, clear

*** Age
gen age = d11101
gen age2 = age^2
gen age3 = age^3

*** Education
gen college =.
replace college = 1 if d11108 == 3
replace college = 0 if d11108 == 1 | d11108 == 2

*** Health (using SRHS and same classification as in De Nardi et al '23, French '05)
gen health_cat=.
replace health_cat = 0 if ple0008 == 4 | ple0008 == 5
replace health_cat = 1 if ple0008 >=1 & ple0008 <=3

*** Survival Rates
gen death2y =.
replace death2y = 0 if todjahr == -2
replace death2y = 1 if todjahr == syear + 1 | todjahr == syear + 2
gen surv2y = 1-death2y


*** Fit Probit model on age cubic conditional on health and education
probit surv2y age age2 age3 if health_cat == 0  & college == 0
predict pi2y_surv_bad_HS

probit surv2y age age2 age3 if health_cat == 0 & college == 1
predict pi2y_surv_bad_CL

probit surv2y age age2 age3 if health_cat == 1  & college == 0
predict pi2y_surv_good_HS

probit surv2y age age2 age3 if health_cat == 1  & college == 1
predict pi2y_surv_good_CL

* Save
egen age_cat2 = cut(age), at(25(2)101)
mat survival_rates = J(38,4,.)
local rr = 1
forvalues a=25(2)99{
	qui sum pi2y_surv_good_HS if age_cat2 == `a'
	mat survival_rates[`rr',1] = r(mean)
	qui sum pi2y_surv_bad_HS if age_cat2 == `a'
	mat survival_rates[`rr',2] = r(mean)
	qui sum pi2y_surv_good_CL if age_cat2 == `a'
	mat survival_rates[`rr',3] = r(mean)
	qui sum pi2y_surv_bad_CL if age_cat2 == `a'
	mat survival_rates[`rr',4] = r(mean)
	
	local rr = `rr'+1
}

cd "$path\workdata"
putexcel set data_figures, sheet(survival_rates) modify
putexcel A1 = mat(survival_rates)

