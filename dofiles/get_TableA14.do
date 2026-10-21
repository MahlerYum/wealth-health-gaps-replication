**********************************************************
******* Calculate initial distribution over fixed states (Table A14)
******

*** 1. Calculate and save productivity fixed effects from wage regression
do "$path\dofiles\get_prodtype.do"

* Merge these with work data
cd "$path\workdata"
use workdata04-18.dta, clear

//drop _merge
merge 1:1 pid syear using work_prodtypes.dta
drop if _merge==2

/*** 2. Use data from time preference experiment in 2006 (Money now or later experiment with 2 list for each participant) to calculate preferences
drop _merge
merge 1:1 pid syear using "$path\rawdata\timepref.dta"
drop if _merge==2
do "$path\dofiles\get_preftype.do"
*/

*** 3. Calculate fixed health types
do "$path\dofiles\get_healthtype.do"

*** 4. Make two productivity types by splitting the distribution of fixed effect estimates in half (following Low and Pistaferri, 2015 AER, we assume that those who never work are of the low productivity type)
local condition "(college !=.) & (health_type !=.) & (age>=25) & (age<65) & (working2y !=.) & (avrz_interest !=.)"
bysort pid: egen avr_work = mean(working2y) if `condition'

gen avr_prod_fe_imp = avr_prod_fe if `condition'
sum avr_prod_fe if `condition'
replace avr_prod_fe_imp = r(min) if avr_work == 0 & `condition'
sum avr_prod_fe_imp

* Condition on self-assessd patience (asked in 2008, 2013, 2018, so create average)
gen patience = plh0253
replace patience =. if patience <0		// zero means very unpatient, 10 very patient
sort pid syear
by pid: egen avr_patience = mean(patience)

local condition2 "(college !=.) & (avr_patience!=.) & (health_type !=.) & (age>=25) & (age<65) & (working2y !=.) & (avrz_interest !=.) & (avr_prod_fe_imp !=.)"

gen prod_type =.
sum avr_prod_fe_imp if `condition2', d
replace prod_type = 0 if avr_prod_fe_imp <r(p50) & `condition2' 
replace prod_type = 1 if avr_prod_fe_imp >=r(p50) & `condition2' 

*** 5. Make two time preference types also based on a 50:50 split of the implied interest rate distribution
gen patience_type = .
su avrz_interest if `condition2', d
replace patience_type = 1 if avrz_interest <r(p50) & `condition2' 
replace patience_type = 0 if avrz_interest >=r(p50) & `condition2' 

*** 6. Get distribution over fixed types (in working age population)
table (patience_type prod_type) (health_type) (college)  if `condition2' 

*** 7. Save in table together with average health and effort among 25-30 year-olds
mat init_distr=J(6,8,.)
matrix colnames init_distr = "$\beta=0_\theta=0_\eta=0$" "$\beta=0_\theta=0_\eta=1$" "$\beta=0_\theta=1_\eta=0$" "$\beta=0_\theta=1_\eta=1$" "$\beta=1_\theta=0_\eta=0$" "$\beta=1_\theta=0_\eta=1$" "$\beta=1_\theta=1_\eta=0$" "$\beta=1_\theta=1_\eta=1$"
matrix rownames init_distr = "Prob_Mass_($e=0$)" "Avg_Health_($e=0$)" "Avg_Effort_($e=0$)" "Prob_Mass_($e=1$)" "Avg_Health_($e=1$)" "Avg_Effort_($e=1$)" 
count if patience_type !=. & health_type !=. & prod_type !=. & college !=.
scalar total_n = r(N)
qui tab prod_type health_type if patience_type == 0 & college == 0, matcell(M)
mat init_distr[1,1]=M[1,1]/total_n
mat init_distr[1,2]=M[1,2]/total_n
mat init_distr[1,3]=M[2,1]/total_n
mat init_distr[1,4]=M[2,2]/total_n
qui tab prod_type health_type if patience_type == 1 & college == 0, matcell(M)
mat init_distr[1,5]=M[1,1]/total_n
mat init_distr[1,6]=M[1,2]/total_n
mat init_distr[1,7]=M[2,1]/total_n
mat init_distr[1,8]=M[2,2]/total_n
qui tab prod_type health_type if patience_type == 0 & college == 1, matcell(M)
mat init_distr[4,1]=M[1,1]/total_n
mat init_distr[4,2]=M[1,2]/total_n
mat init_distr[4,3]=M[2,1]/total_n
mat init_distr[4,4]=M[2,2]/total_n
qui tab prod_type health_type if patience_type == 1 & college == 1, matcell(M)
mat init_distr[4,5]=M[1,1]/total_n
mat init_distr[4,6]=M[1,2]/total_n
mat init_distr[4,7]=M[2,1]/total_n
mat init_distr[4,8]=M[2,2]/total_n

qui sum health_cat if health_type == 0 & age >=25 & age <= 30 & college == 0
scalar m1 = r(mean)
qui sum health_cat if health_type == 1 & age >=25 & age <= 30 & college == 0
scalar m2 = r(mean)
qui sum health_cat if health_type == 0 & age >=25 & age <= 30 & college == 1
scalar m3 = r(mean)
qui sum health_cat if health_type == 1 & age >=25 & age <= 30 & college == 1
scalar m4 = r(mean)
qui sum eff_norm if health_type == 0 & age >=25 & age <= 30 & college == 0
scalar m5 = r(mean)
qui sum eff_norm if health_type == 1 & age >=25 & age <= 30 & college == 0
scalar m6 = r(mean)
qui sum eff_norm if health_type == 0 & age >=25 & age <= 30 & college == 1
scalar m7 = r(mean)
qui sum eff_norm if health_type == 1 & age >=25 & age <= 30 & college == 1
scalar m8 = r(mean)

forvalues cc = 1(2)7{
	mat init_distr[2,`cc']=m1
	mat init_distr[2,`cc'+1]=m2
	mat init_distr[5,`cc']=m3
	mat init_distr[5,`cc'+1]=m4
	mat init_distr[3,`cc']=m5
	mat init_distr[3,`cc'+1]=m6
	mat init_distr[6,`cc']=m7
	mat init_distr[6,`cc'+1]=m8
}

*** Export
outtable using "$path\tables\TableA14", mat(init_distr) replace format (%9.3f) caption("Initial Distribution") nobox center

************************************************
*** For Mortality Simulation (Table A10) export relevant data

tab college health_type if age == 70, matcell(types_70)
putexcel set data_tables, sheet(types_70, replace) modify
putexcel A1 = mat(types_70)

matrix types_7084 = J(4,1,.)
qui sum health_cat if health_type == 0 & college == 0 & age>=70 & age <=84
mat types_7084[1,1] = r(mean)
qui sum health_cat if health_type == 1 & college == 0 & age>=70 & age <=84
mat types_7084[2,1] = r(mean)
qui sum health_cat if health_type == 0 & college == 1 & age>=70 & age <=84
mat types_7084[3,1] = r(mean)
qui sum health_cat if health_type == 1 & college == 1 & age>=70 & age <=84
mat types_7084[4,1] = r(mean)

putexcel set data_tables, sheet(types_7084, replace) modify
putexcel A1 = mat(types_7084)

matrix eff_dec = J(3,1,.)
qui sum eff_norm, d
mat eff_dec[1,1]=r(mean)
mat eff_dec[2,1]=r(p90)
mat eff_dec[3,1]=r(p10)

putexcel set data_tables, sheet(eff_dec, replace) modify
putexcel A1 = mat(eff_dec)