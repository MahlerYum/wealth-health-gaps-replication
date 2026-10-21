**********************************************************
******* Create Working Data Set by Merging different Raw Data sets from SOEP-Core, v35, doi:10.5684/soep-core.v35
******* and performing some sample selection and variable generation
**********************************************************

****** 1. Raw Data Merge

* Change directory to where the raw data is stored
cd "$path/rawdata/"

* Use raw data from individual questionnaire
use pl.dta, clear

* Keep only identifiers and vars pertaining to individual health, health behaviors, and health insurance (+ sexual orientation) and individual risk taking behavior
keep pid hid cid hpmax pnr syear pld0298 ple* plh0188 plh0198 plh0201 plh0204_h plh0226 plh0244 plh0253 plh0254

* Merge this data from pl with generated files from the DIW: pgen (employment, occupation, education, family status), health (generated health indices), pwealth (personal wealth information) and pequiv.dta (generated file to be comparable with other surveys) using pid and syear as unique identifiers
merge 1:1 pid syear using pgen.dta
rename _merge merge_pgen
merge 1:1 pid syear using health.dta
rename _merge merge_health
merge 1:1 pid syear using pwealth.dta
rename _merge merge_wealth
merge 1:1 pid syear using pequiv.dta
rename _merge merge_pequiv


* Keep only data from 2004-2018 (since then, health behaviors are frequently asked)
keep if syear >=2004

****** 2. Variable Generation and further sample selection

*** Age
gen age = d11101
gen age2 = age^2

* Define Age Categories
gen long_age_cat = .
replace long_age_cat =1 if age >=25 & age <45
replace long_age_cat =2 if age >=45 & age <65
replace long_age_cat =3 if age >=65 & age <85
replace long_age_cat =4 if age >=85

egen age_cat10 = cut(age), at(25,35,45,55,65,75,85,99)

*** Education
* College education
gen college =.
replace college = 1 if d11108 == 3
replace college = 0 if d11108 == 1 | d11108 == 2

* Number of years in education (from 7 to 18)
gen yschool = d11109 if d11109 >0

*** Health (using SRHS and same classification as in De Nardi et al '23, French '05)
gen health_cat=.
replace health_cat = 0 if ple0008 == 4 | ple0008 == 5
replace health_cat = 1 if ple0008 >=1 & ple0008 <=3

*** Keep only observations with college, health, and age information
keep if college !=. & health_cat !=. & age >=25


*** Labor Income 
/* Includes the (imputed) gross labor income in the previous year generated for all SOEP respondents who are employed in a main job in each respective wave. Interviewees should consider overtime payments but no irregular one-time payments such as holiday or bonuses. Income details are in Euros and converted to 2016 values using a CPI index 
*/ 
* Create discount factor
gen discount = 100/y11101
gen HH_labinc = i11103 * discount if i11103 >=0
gen labinc = i11110 * discount if i11110 >=0

*** Labor Supply (use definition in pequiv file)
/* 
Strategy:
If the individual had positive wages in the previous year and worked at least 52 hours and had a labor income > 5400 EUR, then the individual was working. Otherwise, the individual was not working.

To aggregate to 2-year employment status, we count a person as working part-time if she worked exactly one of 2 consecutive years full-time. We count a person as working full-time if she worked 2 consecutive years.

As in the model, we normalize full-time working hours to 1, and part-time working hours to 0.5.  
*/
gen working = e11102 if age>=25 & age <65
replace working = 0 if labinc <= 5400 & age>=25 & age <65    
replace working = 0 if age >= 65 

tsset pid syear
gen working2y = 1 if (working == 1 & f.working == 1)
replace working2y = 0.5 if (working == 0 & f.working == 1) | (working == 1 & f.working == 0) 	
replace working2y = 0 if working == 0 & f.working == 0

* We add two consecutive years to get 2-year labor income (to correspond with model frequency)
tsset pid syear
gen labinc2y = labinc + f.labinc if working2y !=.		// only for working age population

* Log labor income (cond. on working) in working age population by education
gen loglabinc2y_workers = log(labinc2y) if working2y >0

* Hours worked 
gen hours = e11101

* Hours worked conditional on working 1 or 2 years among working age population
tsset pid syear
gen hours2y = hours + f.hours if working2y !=.

* Individual Wage of workers
gen wage = labinc/hours
replace wage = . if wage == 0

*** Individual Wealth
* Wealth information is asked from 2002 on every 5 years
* Includes property, real estate, financial assets, savings contracts, business assets, tangible assets, debts

* Average over five different imputation techniques to create net wealth
foreach i in "a" "b" "c" "d" "e" {
gen w0111`i'15 = w0111`i' * discount
}
gen net_wealth = (w0111a15 + w0111b15 + w0111c15 +w0111d15 +w0111e15) / 5  // average over the imputations of net wealth (the 2012 measure)

gen net_wealth_full = net_wealth
* As in Heathcote, Storesletten, Violante (2010, RED), we drop the bottom 0.5% of the wealth observations for moments construction
* Moreover, we group everyone with negative net worth into having zero wealth.
_pctile net_wealth, p(0.5)
scalar cutoff_wealth = r(r1)
replace net_wealth = . if net_wealth <cutoff_wealth
gen net_wealth0 = net_wealth
replace net_wealth0 = 0 if net_wealth <0

* Express wealth in terms of average 2-year labor income (48201 Euros) <- This is an old number but we keep it as a normalization throughout
gen net_wealth0_norm = net_wealth0 / 48201

*** Individual Pension Benefits
* Include old-age, disability, and widowhood social security pensions: payments of the German Pension Insurance (GRV), Miner's social Insurance (Knappschaft), Civil Servant Pension (Beamtenpension), War Victim Benefits (Kriegsopferversorgung), Farmer's Benefits and accident pension (GUV), pensions from abroad and pension for liberal professions (berufsständische Versorgungswerke)
* Pension benefits are recorded at the HH level, so we divide by the number of persons in the HH to get individual pension benefits
gen pension_ind = i11108/d11106 if age>65 & i11108 >0 

*** Further Variables
* Martial status (1: married, 2: single, 3: widowed, 4: divorced, 5: separated)
gen marital_st = d11104 if d11104 >0

* Type of Health Insurance (1: public, 2: private)
gen insurance = ple0097 if ple0097 >0

* Work Status (1: full time, 2: part time, 3: not working)
gen work_st = e11103

* Disability status (in percent)
gen disstat = ple0041
replace disstat = 0 if m11124 == 0

* Days not worked due to sickness
gen noworkdays = ple0046
replace noworkdays = 0 if ple0044_h == 1

* Gender
gen gender = d11102ll

* Birth Cohort
gen birthy = ple0010_h if ple0010_h >=0

* Primary Activity 
gen activity =e11104 if e11104>0

* Highest Education Degree
gen degree = pgisced97 if pgisced97>=0

* Employment Status
gen status = pgemplst
replace status = . if pgemplst <0


*** Generate Individual Health effort
do "$path\dofiles\construct_healtheffort.do"

*********************************************
*** Save data in working data folder
cd "$path\workdata"
save workdata04-18.dta, replace

*********************************************
*** Prepare and add consumption data (available in survey year 2010)
* Consumption expenditures are recorded at the HH level
do "$path\dofiles\create_work_consdata.do"

use workdata04-18.dta, replace

merge m:1 hid syear using work_consdata.dta 
keep if _merge !=2
foreach var in cons_exp_mthly cons_exp_annl health_exp_mthly{
	replace `var' = `var'*discount
}

*** Prepare and add data from time preference experiment in 2006 (Money now or later experiment with 2 list for each participant) to calculate time preferences
drop _merge
merge 1:1 pid syear using "$path\rawdata\timepref.dta"
drop if _merge==2
do "$path\dofiles\get_preftype.do"
drop split teber1 exper teber2 zwl1 zwl2 gzahl zzahl liste zeile gewinn monate pnot1 pnot2 pnot3 pnot4 pnot5 pnot6 pelz01 pelz02 pelz03 pelz05 pelz04 pelz06 pelz07 pelz08 z04 z05ka z05 z06 z07 _merge interest11 interest12 interest21 interest22 interest31 interest32 z_interest11 z_interest12 z_interest21 z_interest22 z_interest31 z_interest32

*** Keep only variables that are used in analysis
keep pid hid syear plh0253 d11107 pglfs m11126 pcs mcs ple0012 ple0011 ple0013 ple0014 ple0015 ple0016 ple0017 ple0018 ple0019 ple0020 ple0021 ple0022 ple0023 m11127 m11102 e11105 d11106 m11107 m11106 m11108 m11111 ple0010_h m11127 m11102 bmi m11125 e11105 age age2 long_age_cat age_cat10 college yschool health_cat discount HH_labinc labinc working working2y labinc2y loglabinc2y_workers hours hours2y wage net_wealth net_wealth_full net_wealth0 net_wealth0_norm pension_ind marital_st insurance work_st disstat noworkdays gender birthy activity degree status ncig_new ncig_new_sd sport sport_new_sd nut_new nut_new_sd InvNcig0 Lsport0 Nut0 Lsport0_resid InvNcig0_resid Nut0_resid eff eff_norm cons_exp_mthly cons_exp_annl health_exp_mthly avrz_interest

*** Save data in working data folder
cd "$path\workdata"
save workdata04-18.dta, replace

*********************************************
*** Prepare data to estimate survival rates
cd "$path\rawdata"

* Use individual data, keep only health info
use pl.dta, clear
keep pid syear ple0008

* Merge with data containing exit and death info
merge 1:1 pid syear using ppathl.dta, keepusing(todjahr)

* Merge with education and age information
drop _merge
merge 1:1 pid syear using pequiv.dta, keepusing(d11101 d11108)

keep pid syear d11101 d11108 ple0008 todjahr

cd "$path\workdata"
save survival_dat.dta, replace




