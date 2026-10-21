**********************************************************
******* Construct fixed Health Types
******

local rseed = 12345		// seed for random number generation

* Use Kmeans clustering on set of health-related moments as detaild in Appendix A6

* First Step: Residualize all info against age, education, gender, Labor Force Status, Labor Income, Marital Status, Health Insurance Type
* Variables are: No of annual doctor visits, self-rated health status (5-point scale), inpatient nights in hospital, PCS, MCS (only every two years), BMI (only every 2 years)

* Record as missing if negative values
local hmoments "m11127 m11126 m11102 pcs mcs bmi"
foreach var of varlist `hmoments'{
	replace `var'=. if `var' <0
}
replace m11125 =. if m11125 <0

* Residualize all variables against age, education, gender, insurance type, birth year and satisfaction with health 
	foreach var of varlist `hmoments' {
		qui reg `var' age age2 college gender i.insurance i.birthy m11125 if `var' >= 0 
		predict `var'_hresid if `var' >=0, residuals
	}

* Second Step: Standardizing residuals to give every variable the chance to be equally important
	foreach var of varlist m11127_hresid-bmi_hresid{
			egen z_`var' = std(`var')
	}

* Generate one average score per individual 
	foreach var of varlist z_m11127_hresid-z_bmi_hresid{
		bysort pid: egen avr`var' = mean(`var')
	}

* Consider only one observation per individual
	foreach var of varlist avrz_m11127_hresid-avrz_bmi_hresid{
	sort pid `var'
	by pid: replace `var' =. if `var' == `var'[1] & _n > 1 
	}

* Third Step: Use kmeans clustering to define clusters of individuals
local listz "avrz_m11127_hresid avrz_m11126_hresid avrz_m11102_hresid avrz_pcs_hresid avrz_mcs_hresid avrz_bmi_hresid"
corr `listz'

* Loop over cluster size
foreach var of varlist `listz'{
	sum `var'
}
cap drop ind_vec
egen ind_vec = rowmiss(`listz')
cap drop cs*
forvalues k = 1(1)5 {
	cluster kmeans `listz' if ind_vec == 0, k(`k') start(random(`rseed')) name(cs`k')
}

* Fourth Step: Calculate and compare the within sum of squares of each cluster solution to find optimal number of clusters as recommended
* WSS matrix
matrix WSS = J(5,5,.)
matrix colnames WSS = k WSS log(WSS) eta-squared PRE
* WSS for each clustering
forvalues k = 1(1)5 {
scalar ws`k' = 0
foreach v of varlist `listz'{
quietly anova `v'  cs`k' 
scalar ws`k'  = ws`k' + e(rss)
}
matrix WSS[`k', 1] = `k'
matrix WSS[`k', 2] = ws`k'
matrix WSS[`k', 3] = log(ws`k'')
matrix WSS[`k', 4] = 1 - ws`k'/WSS[1,2]
matrix WSS[`k', 5] = (WSS[`k'-1,2] - ws`k'')/WSS[`k'-1,2]
}

* Look for "elbow" in WSS 
mat list WSS // largest drop going from 1 to 2 

* Fifth Step: Make the health types based on the 2 clusters and check for age patterns etc.
sum(health_cat) if cs2 == 1 // Check that 0 health type is bad health
scalar m1 = r(mean)
sum health_cat if cs2 == 2
if r(mean) > m1 {
	gen health_type = cs2 - 1
}
else {
	gen health_type = cs2*(-1) + 2
}

* For every individual, health type should be constant over age
sort pid health_type
quietly by pid: replace health_type = health_type[1] if health_type == . & _n >1

* Check age pattern of health type (=> no obvious age pattern)
tab age_cat10, sum(health_type)
tab health_type if age<65 // In working age population, 1/3 are of low health type and 2/3 of high health type
