**********************************************************
******* Construction of Health Effort using PCA of Smoking, Sport, Nutrition (for weights)
**********************************************************

* 1. Construct and standardize all health behaviors
* a. Number of Cigarettes smoked in a day (capped at 50)
gen ncig_new =  ple0086_v2
replace ncig_new = . if ple0086_v2 <=-3 | ple0086_v2 == -1
replace ncig_new = 0 if ncig_new ==-2
replace ncig_new = 50 if ple0086_v2 >50 & ple0086_v2 !=.

qui su ncig_new
gen ncig_new_sd = (ncig_new - r(mean) ) / r(sd)
replace ncig_new_sd = (-1)*ncig_new_sd		// such that high value corresponds to high effort

* b. Frequency of Sport or Exercise (from 1: almost never or never to 4: at least once a week)
gen sport = m11104 if m11104 >0
qui su sport
gen sport_new_sd = (sport - r(mean) ) / r(sd)

* c. Health-conscious nutrition (from 1: very conscious to 4: not at all conscious)
gen nut_new = ple0095
replace nut_new = . if ple0095 <=-1
qui su nut_new
gen nut_new_sd = (nut_new - r(mean) ) / r(sd)
replace nut_new_sd = (-1)*nut_new_sd	// such that high value corresponds to high effort

* 2. Use lagged sport because the question is asked in odd-years only whereas the other effort variables are asked in even years only
tsset pid syear
gen InvNcig0 = ncig_new_sd
gen Lsport0 = l.sport_new_sd
gen Nut0 = nut_new_sd


* 3. Residualize effort variables to strip them from variation coming from age, health, education, marital status, insurance status, work status, labor income, or wealth.
* We do this because we are interested in estimating the contribution of each health behavior to an underlying factor that explains most of the variation in the three behaviors.
* An alternative approach (as in Kling et al, 2007) would be to calculate a healthy behavior index by just equally weighing the three individual behaviors
* As it turns out, our weights are roughly the same anyway.

qui reg Lsport0 health_cat age age2 yschool i.marital i.insurance i.work_st labinc net_wealth, vce(robust)
predict Lsport0_resid, residuals

qui reg InvNcig0 health_cat age age2 yschool i.marital i.insurance i.work_st labinc net_wealth, vce(robust)
predict InvNcig0_resid, residuals

qui reg Nut0 health_cat age age2 yschool i.marital i.insurance i.work_st labinc net_wealth, vce(robust)
predict Nut0_resid, residuals

* 4. Use a principal component analysis to extract common "health behavior" factor as well as the loading of each behavior on this factor
corr Lsport0_resid InvNcig0_resid Nut0_resid 
pca Lsport0_resid InvNcig0_resid Nut0_resid 

* 5. Use factor loadings as weights on sum of three behaviors (as in Cole, Kim, Krueger, 2019, RES)
gen eff = e(L)[1,1]* Lsport0 + e(L)[3,1] * Nut0 + e(L)[2,1] * InvNcig0

* 6. Normalization of effort variable to be between 0 and 1
qui su eff
gen eff_norm = (eff - (r(min))) / (r(max)-r(min))
su eff_norm



