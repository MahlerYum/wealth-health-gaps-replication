**********************************************************
******* Construct fixed Time Preference / Discount Rate Types
******

*** Uses data from time preference experiment in 2006 (Money now or later experiment with 2 list for each participant) 
* Calculate the degree of present bias by calculating the interest rate respondents required for accepting delayed payment
gen interest11 = .
gen interest12 = .
gen interest21 = .
gen interest22 = .
gen interest31 = .
gen interest32 = .
replace interest11 = (200*1.025^zwl1)/200 if split == 1 & zwl1 >0 & zwl1 <=20
replace interest12 = (200*1.0125^zwl2)/200 if split == 1 & zwl2 >0 & zwl2 <=20
replace interest21 = (200*1.025^zwl1)/200 if split == 2 & zwl1 >0 & zwl1 <=20
replace interest22 = (200*1.002^zwl2)/200 if split == 2 & zwl2 >0 & zwl2 <=20
replace interest31 = (200*1.004^zwl1)/200 if split == 3 & zwl1 >0 & zwl1 <=20
replace interest32 = (200*1.004^zwl2)/200 if split == 3 & zwl2 >0 & zwl2 <=20

* Standardize implied interest rate in each list (because lists differed in duration etc.)
forvalues s=1/3{
	forvalues z=1/2{
		egen z_interest`s'`z'= std(interest`s'`z')
	}
}
* Take the average (implied interest rate) of the two lists for an individuals
egen avrz_interest = rowmean(z_interest11-z_interest32)

* Assume preference type (from experiment in 2006) is constant for an individual 
sort pid avrz_interest
quietly by pid: replace avrz_interest = avrz_interest[1] if avrz_interest ==. & _n >1 
sum avrz_interest


