**********************************************************
*******  Estimation of Effect of Health on Work and Labor Income (Table A4)
******

cd "$path\workdata"
use workdata04-18.dta, clear

* Drop non-working age?

* Create working variable from labor force status info
gen working2 = 1 if pglfs == 11
replace working2 = 0 if pglfs !=11
replace working2 = 0 if labinc <= 5400 & age>=25 & age <65   
replace working2 = 0 if age >= 65 

* Labor income and hours conditional on working
gen labinc_work = labinc if working == 1
gen hours_work = hours if working == 1
replace labinc_work =. if hours_work ==.
replace hours_work =. if labinc_work ==.
gen loglabinc = log(labinc_work) 
gen loghours = log(hours_work) 

* Number of annual doctor visits and inpatient nights at hospital
gen doctor_visits = m11127 if m11127 >=0
gen hosp_nights = m11102 if m11102>=0

* Occupation according to 2 digit ISCO classification 
recode e11105 (1110/1319 = 1 "legislators, senior officials and managers") (2111/2470 = 2 "professionals") (3111/3480 = 3 "technicians and ass. professionals") (4111/4223 = 4 "clerks") (5111/5220 = 5 "service and sales workers") (6100/6154 = 6 "skl. agriculture") (7111/7442 = 7 "craft and trade workers") (8111/8340 = 8 "plant and machine operators") (9111/9330 = 9 "elementary occ.") (100 = 10 "soldier"), generate(isco88_2dig)
replace  isco88_2dig = . if isco88_2dig <=0 


* Sample for extensive margin labor supply
tsset pid syear
gen sample_work = 1 if working2 !=. & age !=. & age2 !=. & yschool !=. & gender !=. & marital_st !=. & d11107 !=. & insurance !=. & syear !=. & l.working2 !=. & l2.working2 !=. & health_cat !=. & doctor_visits !=. & hosp_nights !=. & age >=25 & age<65

* Sample for intensive margin labor supply and labor income conditional on working
gen sample_hours = 1 if loghours !=. & age !=. & age2 !=. & yschool !=. & gender !=. & marital_st !=. & d11107 !=. & insurance !=. & syear !=. & l.loghours !=. & l2.loghours !=. & health_cat !=. & doctor_visits !=. & hosp_nights !=. & isco88_2dig !=. & age >=25 & age<65

* Create dummies for all categorical variables
tab marital, g(marry)
tab isco88_2dig, g(occ)
replace insurance = insurance - 1

keep if age < 65
mat tabA4=J(3,3,.)

preserve

* Individual FE IV Estimates
* Need to demean all variables (omit time constant ones)
keep if sample_work == 1
xtdata working working2 age age2 yschool insurance health_cat doctor_visits hosp_nights syear d11107 marry*, fe clear i(pid)

gen syear2 = round(syear)	// to get lags of working 
xtset pid syear2
ivregress gmm working2 age age2 d11107 insurance l.working2 l2.working2 (health_cat= doctor_visits hosp_nights) marry* i.syear2, vce(robust)

mat tabA4[1,1]=e(b)[1,1]
mat tabA4[2,1]=sqrt(e(V)[1,1])
mat tabA4[3,1]=e(N)

restore

* Same for hours and labor income conditional on work (include occupation dummy)
keep if working == 1 & sample_hours == 1
xtdata loglabinc loghours age age2 yschool insurance health_cat doctor_visits hosp_nights syear d11107 marry* occ*, fe clear i(pid)
gen syear2 = round(syear)
xtset pid syear2

ivregress gmm loglabinc age age2 d11107 insurance l.loglabinc l2.loglabinc (health_cat= doctor_visits hosp_nights) occ* i.syear2, vce(robust)
mat tabA4[1,2]=e(b)[1,1]
mat tabA4[2,2]=sqrt(e(V)[1,1])
mat tabA4[3,2]=e(N)

ivregress gmm loghours age age2 d11107 insurance l.loghours l2.loghours (health_cat= doctor_visits hosp_nights) occ* i.syear2, vce(robust)
mat tabA4[1,3]=e(b)[1,1]
mat tabA4[2,3]=sqrt(e(V)[1,1])
mat tabA4[3,3]=e(N)

matrix colnames tabA4 = "work" "log(income|work)" "log(hours|work)"
matrix rownames tabA4 = "Health" "SE" "N"

outtable using "$path\tables\TableA4", mat(tabA4) replace format (%9.3f) caption("Effect of Health on Work Status, Labor Income, and Hours Worked") nobox center






