**********************************************************
*******  Estimation of Relation between Wealth, Health and Effort (Table 1)
******

cd "$path\workdata"
use workdata04-18.dta, clear

xtset pid syear

* Impute effort for every year
gen odd = mod(syear,2)
gen leff = l.eff_norm if odd == 1
gen feff = f.eff_norm if odd == 1
egen eff_odd = rowmean(leff feff) // now have eff in 2005 and 2007 ..., 2015
gen eff_comp = eff_norm if odd == 0
replace eff_comp = eff_odd if odd == 1
gen wealth_tenth = net_wealth_full /10000

* Create average patience by individual from self-assessd patience (only asked in 2008, 2013, 2018)
gen patience = plh0253
replace patience =. if patience <0	
sort pid syear
by pid: egen avr_patience = mean(patience)	
tab avr_patience, g(patdummy)

* Labor income inlusive zeros for retirees
gen labinc0 = labinc
replace labinc0 = 0 if labinc0 == . & working == 0
replace labinc0 = 0 if age >= 65

reg health_cat wealth_tenth eff_comp age age2 l.health_cat yschool labinc0 hours i.pglfs i.gender i.marital_st i.insurance d11107 i.syear patdummy*, vce(robust)
gen byte used = e(sample)
gen sample = 0
replace sample = 1 if (used == 1 & f.health_cat !=. & f2.health_cat !=. & f3.health_cat !=.)

matrix est_wealth_health = J(9,4,.)
matrix colnames est_wealth_health = "$k=0$" "$k=1$" "$k=2$" "$k=3$"
matrix rownames est_wealth_health = "Wealth" "(SE)" "Wealth" "(SE)" "Effort" "(SE)" "(1) R^2" "(2) R^2" "N"


forvalues k=0/3{
qui reg f`k'.health_cat wealth_tenth age age2 l.health_cat yschool labinc0 hours i.pglfs i.gender i.marital_st i.insurance d11107 patdummy* i.syear if sample == 1, vce(robust)
mat est_wealth_health[1,`k'+1] = e(b)[1,1]*1000
mat est_wealth_health[2,`k'+1] = sqrt(e(V)[1,1])*1000
mat est_wealth_health[7,`k'+1] = e(r2)

qui reg f`k'.health_cat wealth_tenth eff_comp age age2 l.health_cat yschool labinc0 hours i.pglfs i.gender i.marital_st i.insurance d11107 i.syear patdummy* if sample == 1, vce(robust)
mat est_wealth_health[3,`k'+1] = e(b)[1,1]*1000
mat est_wealth_health[4,`k'+1] = sqrt(e(V)[1,1])*1000
mat est_wealth_health[8,`k'+1] = e(r2)
mat est_wealth_health[5,`k'+1] = e(b)[1,2]
mat est_wealth_health[6,`k'+1] = sqrt(e(V)[2,2])
mat est_wealth_health[9,`k'+1] = e(N)
}

outtable using "$path\tables\Table1", mat(est_wealth_health) replace format (%9.3f) caption("Effect of Wealth on Current and Future Health, with and without Effort") nobox center

