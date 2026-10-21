**********************************************************
*******  Export Data to create Figures (using Matlab)
******

cd "$path\workdata"
use workdata04-18.dta, clear

*** Save all moments in excel file
//putexcel set data_figures

*** Average Health by Age and Education
mat health10 = J(7,2,.)

forvalues aa = 25(10)85{
	forvalues cc =0/1{
	local rr = `aa'/10 - 1.5
		qui sum health_cat if age_cat10 == `aa' & college == `cc'
		mat health10[`rr',`cc'+1] = 1-r(mean)
	} 
}

putexcel set data_figures, sheet(health10, replace) modify
putexcel A1 = mat(health10)

*** Average Health Effort by Age, Health, and Education
mat eff10_edu_health = J(7,4,.)

forvalues aa = 25(10)85{
	forvalues cc =0/1{
		forvalues hh = 0/1{
	local rr = `aa'/10 - 1.5
		qui sum eff_norm if age_cat10 == `aa' & college == `cc' & health_cat == `hh'
		if `hh'== 0 {
			mat eff10_edu_health[`rr',`cc'+1] = r(mean)
		}
		else {
			mat eff10_edu_health[`rr',`cc'+3] = r(mean)
		}
		}
	} 
}

putexcel set data_figures, sheet(eff10_edu_health, replace) modify
putexcel A1 = mat(eff10_edu_health)

*** Median Wealth by Age, Health, and Education
mat wealth10_50_edu = J(7,4,.)
bysort age_cat10 college health_cat: egen we10_50_edu_health = median(net_wealth0)

forvalues aa = 25(10)85{
	forvalues cc =0/1{
		forvalues hh = 0/1{
	local rr = `aa'/10 - 1.5
		qui sum we10_50_edu_health if age_cat10 == `aa' & college == `cc' & health_cat == `hh'
		if `cc'== 0 {
			mat wealth10_50_edu[`rr',`hh'+1] = r(mean)
		}
		else {
			mat wealth10_50_edu[`rr',`hh'+3] = r(mean)
		}
		}
	} 
}

putexcel set data_figures, sheet(wealth10_50_edu, replace) modify
putexcel A1 = mat(wealth10_50_edu)

*** Average Effort by Education-specific Wealth Quartiles
gen wealth_quant = . 
forvalues a=1/3 {
    forvalues c=0/1{
quietly _pctile net_wealth0 if long_age_cat == `a' & college == `c', n(4)
replace wealth_quant = 1 if net_wealth0 <= r(r1) & long_age_cat ==`a' & college == `c'
replace wealth_quant = 2 if net_wealth0 > r(r1) & net_wealth0 <= r(r2) & long_age_cat ==`a' & college == `c'
replace wealth_quant = 3 if net_wealth0 > r(r2) & net_wealth0 <= r(r3) & long_age_cat ==`a' & college == `c'
replace wealth_quant = 4 if net_wealth0 > r(r3) & net_wealth0 !=. & long_age_cat == `a' & college == `c'
}
}

mat effwealth_HS = J(3,4,.)
mat effwealth_CL = J(3,4,.)

forvalues a=1/3{
	forvalues wq = 1/4{
		qui sum eff_norm if long_age_cat == `a' & wealth_quant == `wq' & college == 0
		mat effwealth_HS[`a',`wq']=r(mean)
		qui sum eff_norm if long_age_cat == `a' & wealth_quant == `wq' & college == 1
		mat effwealth_CL[`a',`wq']=r(mean)
	}
}

putexcel set data_figures, sheet(effwealth_HS, replace) modify
putexcel A1 = mat(effwealth_HS)
putexcel set data_figures, sheet(effwealth_CL, replace) modify
putexcel A1 = mat(effwealth_CL)

*** Average of each component of effort by age and education
gen temp_age = age_cat10
replace temp_age = 75 if temp_age == 85

mat smoking = J(6,2,.)
mat sport = J(6,2,.)
mat nutrition = J(6,2,.)
forvalues aa = 25(10)75{
	local rr = `aa'/10 - 1.5
	forvalues cc =0/1{
		qui sum InvNcig0 if temp_age == `aa' & college == `cc'
		mat smoking[`rr',`cc'+1]=r(mean)
				qui sum Lsport0 if temp_age == `aa' & college == `cc'
		mat sport[`rr',`cc'+1]=r(mean)
				qui sum Nut0 if temp_age == `aa' & college == `cc'
		mat nutrition[`rr',`cc'+1]=r(mean)
	}
}

putexcel set data_figures, sheet(smoking, replace) modify
putexcel A1 = mat(smoking)
putexcel set data_figures, sheet(sport, replace) modify
putexcel A1 = mat(sport)
putexcel set data_figures, sheet(nutrition, replace) modify
putexcel A1 = mat(nutrition)

*** Labor Supply Moments
mat frac_work10 = J(4,2,.)
forvalues a=25(10)55{
	local rr = `a'/10-1.5
	forvalues hh = 0/1{
		qui sum working2y if age_cat10 == `a' & health_cat == `hh'
		mat frac_work10[`rr',`hh'+1]=r(mean)
	}
}

putexcel set data_figures, sheet(frac_work10, replace) modify
putexcel A1 = mat(frac_work10)

mat frac_work_edu10 = J(4,2,.)
forvalues a=25(10)55{
	local rr = `a'/10-1.5
	forvalues cc = 0/1{
		qui sum working2y if age_cat10 == `a' & college == `cc'
		mat frac_work_edu10[`rr',`cc'+1]=r(mean)
	}
}

putexcel set data_figures, sheet(frac_work_edu10, replace) modify
putexcel A1 = mat(frac_work_edu10)

*** Labor Income
mat avrg_inc10_edu_health = J(4,4,.)
forvalues a=25(10)55{
	local rr = `a'/10-1.5
	forvalues cc = 0/1{
		forvalues hh = 0/1{		
		qui sum labinc2y if age_cat10 == `a' & college == `cc' & health_cat == `hh'
		if `cc'==0{
			mat avrg_inc10_edu_health[`rr',`hh'+1]=r(mean)
		}
		else {
		mat avrg_inc10_edu_health[`rr',`hh'+3]=r(mean)	
		}	
		}
	}
}

putexcel set data_figures, sheet(avrg_inc10_edu_health, replace) modify
putexcel A1 = mat(avrg_inc10_edu_health)

*** Wealth Moments
* Unconditional Quartiles
//gen temp_age = age_cat10
//replace temp_age = 75 if age_cat10 == 85
bysort temp_age: egen we10_50_free = median(net_wealth0)
gen we10_25_free = .
gen we10_75_free = .
forvalues a=25(10)75{
	    qui su net_wealth0 if temp_age == `a', d
		replace we10_25_free = r(p25) if temp_age == `a'
		replace we10_75_free = r(p75) if temp_age == `a'
}

mat wealth10_dat = J(6,3,.)
forvalues a=25(10)75{
	local rr = `a'/10-1.5
	qui sum we10_25_free if temp_age == `a'
	mat wealth10_dat[`rr',1] = r(mean) / 1000
	qui sum we10_50_free if temp_age == `a'
	mat wealth10_dat[`rr',2] = r(mean) / 1000
	qui sum we10_75_free if temp_age == `a'
	mat wealth10_dat[`rr',3] = r(mean) / 1000
}

putexcel set data_figures, sheet(wealth10_dat, replace) modify
putexcel A1 = mat(wealth10_dat)

* Conditional Quartiles
bysort age_cat10 health_cat: egen we10_50_0 = median(net_wealth0)
gen we10_25_0 = .
gen we10_75_0 = .
forvalues a=25(10)85{
	forvalues h=0/1{
	    qui su net_wealth0 if age_cat10 == `a' & health_cat == `h', d
		replace we10_25_0 = r(p25) if age_cat10 == `a' & health_cat == `h'
		replace we10_75_0 = r(p75) if age_cat10 == `a' & health_cat == `h'
	}
}

mat wealth10_25 = J(7,2,.)
mat wealth10_50 = J(7,2,.)
mat wealth10_75 = J(7,2,.)

forvalues a=25(10)85{
	local rr = `a'/10-1.5
	forvalues hh = 0/1{
	qui sum we10_25_0 if age_cat10 == `a' & health_cat == `hh'
	mat wealth10_25[`rr',`hh'+1] = r(mean) / 1000
	qui sum we10_50_0 if age_cat10 == `a' & health_cat == `hh'
	mat wealth10_50[`rr',`hh'+1] = r(mean) / 1000
	qui sum we10_75_0 if age_cat10 == `a' & health_cat == `hh'
	mat wealth10_75[`rr',`hh'+1] = r(mean) / 1000
	}
}

putexcel set data_figures, sheet(wealth10_25, replace) modify
putexcel A1 = mat(wealth10_25)

putexcel set data_figures, sheet(wealth10_50, replace) modify
putexcel A1 = mat(wealth10_50)

putexcel set data_figures, sheet(wealth10_75, replace) modify
putexcel A1 = mat(wealth10_75)

* Conditional on Occupation (Manual vs. Non-Manual)
* Occupation according to 2 digit ISCO classification 
recode e11105 (1110/1319 = 1 "legislators, senior officials and managers") (2111/2470 = 2 "professionals") (3111/3480 = 3 "technicians and ass. professionals") (4111/4223 = 4 "clerks") (5111/5220 = 5 "service and sales workers") (6100/6154 = 6 "skl. agriculture") (7111/7442 = 7 "craft and trade workers") (8111/8340 = 8 "plant and machine operators") (9111/9330 = 9 "elementary occ.") (100 = 10 "soldier"), generate(isco88_2dig)
replace  isco88_2dig = . if isco88_2dig <=0 

* Manual Occupations include argricultre, craft and trade, plant and machine operators, and other elementary occupations. Non-Manual occupactions are all other occupations
gen manual = 0 if isco88_2dig >=1 & isco88_2dig !=.
replace manual = 1 if isco88_2dig >=6 & isco88_2dig <=9

bysort age_cat10 manual health_cat: egen medweocc = median(net_wealth0)

* Save
mat wealth10_50_occ = J(7,4,.)

forvalues a=25(10)85{
	local rr = `a'/10-1.5
	forvalues mm = 0/1{
		if `mm' == 0 {
			scalar cc = 1
		}
		else {
			scalar cc = 3
		}
	forvalues hh = 0/1{

	qui sum medweocc if age_cat10 == `a' & health_cat == `hh' & manual == `mm'
	mat wealth10_50_occ[`rr',`hh'+cc] = r(mean) / 1000

	}
	}
}

putexcel set data_figures, sheet(wealth10_50_occ, replace) modify
putexcel A1 = mat(wealth10_50_occ)


*** Save
putexcel save