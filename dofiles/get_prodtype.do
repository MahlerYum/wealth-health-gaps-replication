**********************************************************
******* Construct individual fixed Productivity Types
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

* Get only individual fixed effect
predict prod_fe, u

* Residualize against year of birth (to get rid of cohort effects) and college
reg prod_fe ple0010_h college if ple0010_h >0
predict prod_fe_resid if ple0010_h >0, residuals

* Average per individual (as in Low and Pistaferri)
bysort pid: egen avr_prod_fe = mean(prod_fe_resid)

* Save as work data
keep prod_fe_resid avr_prod_fe age pid syear
cd "$path\workdata"
save work_prodtypes.dta, replace