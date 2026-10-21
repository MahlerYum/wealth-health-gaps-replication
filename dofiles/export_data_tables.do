*********************************************************
*******  Export Data to create Tables (using Matlab)
******

cd "$path\workdata"
use workdata04-18.dta, clear

*** Save all moments in excel file
//putexcel set data_tables

*** Effort Adjustments 
mat eff_adj = J(12,1,.)

* Effort Adjustments 
tsset pid syear
gen delta_eff = eff_norm - l2.eff_norm
gen delta_eff_rel = eff_norm/l2.eff_norm

gen pos_spike10 = 1 if delta_eff_rel >1.1 & delta_eff_rel !=.
replace pos_spike10 = 0 if delta_eff_rel <=1.1

gen neg_spike10 = 1 if delta_eff_rel <0.9 & delta_eff_rel !=.
replace neg_spike10 = 0 if delta_eff_rel >=0.9 & delta_eff_rel !=.

gen pos_spike20 = 1 if delta_eff_rel >1.2 & delta_eff_rel !=.
replace pos_spike20 = 0 if delta_eff_rel <=1.2

gen neg_spike20 = 1 if delta_eff_rel <0.8 & delta_eff_rel !=.
replace neg_spike20 = 0 if delta_eff_rel >=0.8 & delta_eff_rel !=.

forvalues aa = 1/3{
		qui sum pos_spike10 if long_age_cat == `aa'
		mat eff_adj[`aa',1] = r(mean)
		qui sum pos_spike20 if long_age_cat == `aa'
		mat eff_adj[`aa'+3,1] = r(mean)
		qui sum neg_spike10 if long_age_cat == `aa'
		mat eff_adj[`aa'+6,1] = r(mean)
		qui sum neg_spike20 if long_age_cat == `aa'
		mat eff_adj[`aa'+9,1] = r(mean)
}

putexcel set data_tables, sheet(eff_adj, replace) modify
putexcel A1 = mat(eff_adj)