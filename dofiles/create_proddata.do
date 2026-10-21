**********************************************************
******* Generate data to export for productivity process estimation (AR(1) coefficient) by estimatin a fixed effects regression on wages of workers 
******

cd "$path\workdata"
use workdata04-18.dta, clear

* Keep only fulltime workers (i.e. workers that work 2 consecutive years)
drop if working2y !=1

* Estimate a fixed effects regression like De Nardi et al. 2023 (equation 1) but on hourly wages 
* Generate wages as in LMPQS (Lagakos et al. 2018, JPE)
gen wage2y = labinc2y / hours2y

replace wage2y = . if wage2y == 0
sort syear 
by syear: egen p_high = pctile(wage2y), p(99)
by syear: egen p_low = pctile(wage2y), p(1)

replace wage2y = . if wage2y>p_high | wage<p_low

* Take logs
generate lwage = log(wage2y)
* Create age and health interactions dummies (De Nardi et al. 2022, equation 1)
egen group_agehealth = group(age health_cat)

* Estimate fixed effects model  but on log wages instead of earnings 
xtreg lwage i.group_agehealth, fe

* Get combined residual
predict combresid, ue

* Residualize against year of birth (to get rid of cohort effects) and college
reg combresid i.ple0010_h college if ple0010_h >0
predict combresid_age if ple0010_h >0, resid

* Reshape and export the datasset for wage process shock esitmation 
keep if combresid_age !=.
keep combresid_age age pid
reshape wide combresid_age, i(pid) j(age)

cd "$path\workdata"
export delimited combresid_age* using "prod_data", replace