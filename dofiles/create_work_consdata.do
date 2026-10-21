**********************************************************
******* Prepare a consumption data set to use for estimation of Table A1 and Consumption Moments in Table A5
******* Detailed consumption expenditures are asked only in survey wave 2010
**********************************************************

* Change directory to where the raw data is stored
cd "$path/rawdata/"

* Use raw data from individual questionnaire
use hconsum.dta, clear

/* Construct monthly HH consumption expenditures using info on
- Food at home
- Food outside
- Clothing/Shoes
- Body Care
- Health
- Telecommunication
- Education/child care
- Culture
- Leisure activities
- Vacation
- Car repairs
- Transport
- Furntiture
- Home help
- Rent
- Electricity
- Heating
- Other housing costs
- life insurance
- other insurances
- mortgage payments
- other loans
- alimony
- other expenditures
*/

* First, average each consumption category over the 5 imputation versions a-e 
foreach var in rheat_ rrent_ relectr_ oheat_ oelectr_ outil_ consum1 consum2 consum3 consum4 consum5 consum6 consum7 consum8 consum9 consum10 consum11 consum12 consum13 consum14 consum15 consum16 mortgage_ housing_ loan_ alimony_ clhelp_ chcare_ {
	gen avr_`var' = 1/5*(`var'a + `var'b + `var'c + `var'd + `var'e)
}

* Since oheat oelectr outil reflect annual consumption expenses, divide them by 12
replace avr_oheat = avr_oheat/12
replace avr_oelectr = avr_oelectr/12
replace avr_outil = avr_outil/12

* Create total monthly consumption expenditures
egen cons_exp_mthly = rowtotal(avr_*)

* Subtract private pensions/insurance and loans + mortgage payments because these are already in net wealth
replace cons_exp_mthly = cons_exp_mthly - avr_consum11 - avr_loan

* Ge annual consumption and health-relatec monthly consumption
gen cons_exp_annl = 12* cons_exp_mthly
gen health_exp_mthly = avr_consum5

* Save data set
keep hid syear cid cons_exp_mthly cons_exp_annl health_exp_mthly

cd "$path\workdata"
save work_consdata.dta, replace