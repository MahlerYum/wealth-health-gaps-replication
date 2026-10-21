**********************************************************
******* Create Figures A1, A2 and Table A2
******

cd "$path\workdata"
use workdata04-18.dta, clear

*** Table A2
* Correlation between binary health and other health measures: SRHS on 5-point scale, frailty index, and pcs and mcs scores
gen srhs = m11126 if m11126 >=0 // order is flipped (i.e. from 1: very good to 5: bad)
replace srhs = -srhs
replace pcs =. if pcs <0
replace mcs =. if mcs <0

* Generate Frailty index by using info on disease diagnostics plus diffulties (same strategy as in Hosseini, Kopecky, and Zhao, 2023)
foreach var in ple0012 ple0011 ple0013 ple0014 ple0015 ple0016 ple0017 ple0018 ple0019 ple0020 ple0021 ple0022 ple0023{
		gen `var'_norm = `var'
	replace `var'_norm = . if `var'_norm <-2 | `var'_norm >1
	replace `var'_norm = 0 if `var'_norm == -2
}
egen ind_frailty = rowmiss(ple0012_norm - ple0023_norm)
egen frailty_new = rowtotal(ple0012_norm - ple0023_norm) if ind_frailty == 0

* Standardize each measure
foreach var in health_cat srhs frailty_new pcs mcs{
	egen `var'_std = std(`var')
}

* Create matrix to store all empirical moments required for the MSM estimation
matrix correlations = J(5,5,.)		

tsset pid syear
corr health_cat_std srhs_std f.frailty_new_std pcs_std mcs_std

mat correlations = r(C)
matrix colnames correlations = "Binary_Health" "5-point_SRHS" "Frailty" "PCS" "MCS" 
matrix rownames correlations = "Binary_Health" "5-point_SRHS" "Frailty" "PCS" "MCS"

outtable using "$path\tables\TableA2", mat(correlations) replace format (%9.2f) caption("Correlations across Different Health Measures") nobox center

*** Figures A1 and A2 Make Graphs (Appendix Figures A.1 + A.2)
 label define agegroup 25 "25-34" 35 "35-44" 45 "45-54" 55 "55-64" 65 "65-74" 75 "75-84" 85 "85+"
label valu age_cat10 agegroup

graph box pcs, over(age_cat10) nooutside title("Physical Condition Index") note("") 
graph export "$path\figures\FigureA1pcs.eps", replace

graph box mcs, over(age_cat10) nooutside title("Mental Condition Index") note("")
graph export "$path\figures\FigureA1mcs.eps", replace

bysort age_cat10: egen mean_frailty = mean(frailty_new)
twoway (line mean_frailty age_cat10, xtitle("Age Group") ytitle("Frailty Index") xlabel(25 "25-34" 35 "35-44" 45 "45-54" 55 "55-64" 65 "65-74" 75 "75-84" 85 "85+"))
graph export "$path\figures\FigureA2.eps", replace
