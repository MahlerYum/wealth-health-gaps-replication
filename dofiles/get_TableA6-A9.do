**********************************************************
******* Estimate Health Transition Technology (Table A6)
******

cd "$path\workdata"
use workdata04-18.dta, clear

*** Get fixed health types
do "$path\dofiles\get_healthtype.do"

*** Estimate health technology using logit model and group fixed effects
mat tabA6 = J(14,4,.)
matrix colnames tabA6 = "End_Health_Estimate" "End_Health_SE" "Ex_Health_Estimate" "Ex_Health_SE"
matrix rownames tabA6 = "$\lambda_1$" "$\lambda_2$" "$h=1$" "Age_Group" "35" "45" "55" "65" "75" "$\eta=1$" "e=1$" "Constant" "Pseudo_R2" "Number_of_Obs"

sort pid syear

* One coefficient for age group 75+ (as in paper)
gen age_cat10_bunch = age_cat10
replace age_cat10_bunch = 75 if age_cat10_bunch == 85

logit f2.health_cat eff_norm l2.eff_norm health_cat i.age_cat10_bunch health_type college, vce(robust) 

* Predict probabilities to be healthy in data (needed for validation)
predict p2healthy_est, pr

* Save in matrix
forvalues kk=1/12{
	mat tabA6[`kk',1]=e(b)[1,`kk']
	mat tabA6[`kk',2]=sqrt(e(V)[`kk',`kk'])
}
mat tabA6[13,1]=e(r2_p)
mat tabA6[14,1]=e(N)


*** Estimate exogenous health model (i.e. without effort)
logit f2.health_cat health_cat i.age_cat10_bunch health_type college if e(sample) == 1, vce(robust) 

forvalues kk=3/12{
	mat tabA6[`kk',3]=e(b)[1,`kk'-2]
	mat tabA6[`kk',4]=sqrt(e(V)[`kk'-2,`kk'-2])
}
mat tabA6[13,3]=e(r2_p)
mat tabA6[14,3]=e(N)

forvalues cc=1/4{
	mat tabA6[4,`cc']=.
}

*** Export
outtable using "$path\tables\TableA6", mat(tabA6) replace format (%9.3f) caption("Logit Estimation of Probability of being Healthy in 2 years") nobox center

**********************************************************
******* Calculate Average Marginal Effects (Table A7)
*******

logit f2.health_cat eff_norm l2.eff_norm health_cat i.age_cat10_bunch health_type college, vce(robust) 

* Calculate average marginal effects conditioinal on covariates
gen la0_pr =.
gen la1_pr =.
gen la2_pr =.
forvalues a=25(10)75{
    forvalues h=0/1{
		forvalues d=0/1{
				forvalues c=0/1{
					qui margins, at(eff=0 l2.eff=0 health_cat = `h' age_cat10 = `a' health_type = `d' college = `c')
					replace la0_pr = r(b)[1,1] if health_cat == `h' & age_cat10 == `a' & health_type == `d' & college==`c'
					qui margins, dydx(eff l2.eff) at(health_cat = `h' age_cat10 = `a' health_type = `d' college = `c')
					replace la1_pr = r(b)[1,1] if health_cat == `h' & age_cat10 == `a' & health_type == `d' & college==`c'
					replace la2_pr = r(b)[1,2] if health_cat == `h' & age_cat10 == `a' & health_type == `d' & college==`c'
				}
		}
	}
}

mat tabA7 = J(12,12,.)

forvalues a=25(10)75{
    forvalues h=0/1{
		forvalues d=0/1{
				forvalues c=0/1{
					
					if `d'==0 {
					local rr = `a'/10 - 1.5
					}
					else {
					local rr = `a'/10 + 4.5
					}
					
					if `c'==0 {
						if `h'== 0{
							local cc = 1
						}
						else {
							local cc = 4
						}
						
					} 
					else  {
						if `h'== 0{
							local cc = 7
						}
						else {
							local cc = 10
						}
						
					}
					
					qui sum la0_pr if health_cat == `h' & age_cat10 == `a' & health_type == `d' & college == `c'			
					mat tabA7[`rr',`cc']=r(mean)
					
					qui sum la1_pr if health_cat == `h' &age_cat10 == `a' & health_type == `d' & college == `c'			
					mat tabA7[`rr',`cc'+1]=r(mean)
					
					qui sum la2_pr if health_cat == `h' & age_cat10 == `a' & health_type == `d' & college == `c'
					mat tabA7[`rr',`cc'+2]=r(mean)
				}
			}
		}
	}

matrix colnames tabA7 = "$\pi^0$" "$\lambda_1$" "$\lambda_2$" "$\pi^0$" "$\lambda_1$" "$\lambda_2$" "$\pi^0$" "$\lambda_1$" "$\lambda_2$" "$\pi^0$" "$\lambda_1$" "$\lambda_2$" 
matrix rownames tabA7 = "25-34" "35-44" "45-54" "55-64" "65-74" "75+" "25-34" "35-44" "45-54" "55-64" "65-74" "75+"        

outtable using "$path\tables\TableA7", mat(tabA7) replace format (%9.2f) caption("Average Marginal Effects from Health Technology Estimates") nobox center

**********************************************************
******* Summarize Disease Prevalance by Age Group, Education, and Health Status (Table A8)
*******
*** Prevalance of Diabetes, Cancer, Hypertension, and Heart Condition
mat tabA8=J(12,8,.)
 
	forvalues h=0/1{
		forvalues c=0/1{	
			
				if `h'==0 {
					local rr = 1
					}
				else {
					local rr = 2
				}
					
				if `c'==0 {
					local cc = 1	
				} 
				else  {
					local cc = 2						
				}
				
				forvalues a=25(10)75{ 
			qui sum m11107 if m11107>=0 & health_cat == `h' & age_cat10_bunch == `a' & college == `c'			
				mat tabA8[`rr',`cc']=r(mean)
			qui sum m11108 if m11108>=0 & health_cat == `h' & age_cat10_bunch == `a' & college == `c'			
				mat tabA8[`rr',`cc'+2]=r(mean)
			qui sum m11106 if m11106>=0 & health_cat == `h' & age_cat10_bunch == `a' & college == `c'			
				mat tabA8[`rr',`cc'+4]=r(mean)
			qui sum m11111 if m11111>=0 & health_cat == `h' & age_cat10_bunch == `a' & college == `c'			
				mat tabA8[`rr',`cc'+6]=r(mean)
				local rr = `rr'+2
						
		}
	}	
}

matrix colnames tabA8 = "NCl_Diabetes" "Cl_Diabetes" "NCl_Cancer" "Cl_Cancer" "NCl_Hypertension" "Cl_Hypertension" "NCl_Heart_Condition" "Cl_Heart_Condition" 
matrix rownames tabA8 = "25-34" "35-44" "45-54" "55-64" "65-74" "75+" "25-34" "35-44" "45-54" "55-64" "65-74" "75+"  

outtable using "$path\tables\TableA8", mat(tabA8) replace format (%9.3f) caption("Prevalence of Diseases in Population by Age Group and Health Status") nobox center

*** Table A9 Implied Probability of Disease by Past and Current Effort Tercile
gen prev_diab_good = .
gen prev_diab_bad =.
gen prev_canc_good =.
gen prev_canc_bad = .
gen prev_hyp_good = .
gen prev_hyp_bad =.
gen prev_heart_good = .
gen prev_heart_bad = .

forvalues a=25(10)75{
	    forvalues c=0/1{
			qui su m11107 if m11107>=0 & health_cat == 1 & age_cat10 == `a' & college == `c'
			replace prev_diab_good = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11107 if m11107>=0 & health_cat == 0 & age_cat10 == `a' & college == `c'
			replace prev_diab_bad = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11108 if m11108>=0 & health_cat == 1 & age_cat10 == `a' & college == `c'
			replace prev_canc_good = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11108 if m11108>=0 & health_cat == 0 & age_cat10 == `a' & college == `c'
			replace prev_canc_bad = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11106 if m11106>=0 & health_cat == 1 & age_cat10 == `a' & college == `c'
			replace prev_hyp_good = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11106 if m11106>=0 &health_cat == 0 & age_cat10 == `a' & college == `c'
			replace prev_hyp_bad = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11111 if m11111>=0 & health_cat == 1 & age_cat10 == `a' & college == `c'
			replace prev_heart_good = r(mean) if age_cat10 == `a' & college == `c'
			qui su m11111 if m11111>=0 & health_cat == 0 & age_cat10 == `a' & college == `c'
			replace prev_heart_bad = r(mean) if age_cat10 == `a' & college == `c'
		}
}


* Average precent deviation from status-mean of disease by status and effort tercile (using estimated health transition coefficients)
cap drop dev_* prob_*
foreach disease in diab canc hyp heart{
gen dev_`disease' =.
gen dev_`disease'_lf =.
gen dev_`disease'_bf =.
gen prob_`disease'=.
}

* Average implied probability of onset of specific disease using model-implied probability of being healthy in t+2
foreach disease in diab canc hyp heart{
	replace prob_`disease' = p2healthy_est * prev_`disease'_good + (1-p2healthy_est) * prev_`disease'_bad
}

* Construct effort terciles by age-group, current health and college
cap drop eff_terc
gen eff_terc = .
forvalues a=25(10)75{
    forvalues h=0/1{
	    forvalues c=0/1{
			_pctile eff_norm if health_cat == `h' & age_cat10 == `a' & college == `c', n(3)
			replace eff_terc = 1 if eff_norm <= r(r1) & health_cat == `h' & age_cat10 == `a' & college == `c'
			replace eff_terc = 2 if eff_norm >r(r1) & eff_norm <= r(r2) & health_cat == `h' & age_cat10 == `a' & college == `c'
			replace eff_terc = 3 if eff_norm >r(r2) & eff_norm !=. & health_cat == `h' & age_cat10 == `a' & college == `c'
		}
	}
}

* Calculate percent deviation of probability of disease onset within status for different effort terciles
tsset pid syear
forvalues a=25(10)65{
    forvalues h=0/1{
	    forvalues c=0/1{
			foreach disease in diab canc hyp heart{
				qui sum prob_`disease' if health_cat == `h' & age_cat10 == `a' & college == `c' 
				scalar mm_`disease' = r(mean)
				
				* Calculate percent deviations of average disease probability in each effort tercile from mean
				forvalues ee=1/3{
					qui sum prob_`disease' if health_cat == `h' & age_cat10 == `a' & college == `c' & eff_terc == `ee'	
					replace dev_`disease' = (r(mean) / mm_`disease' - 1)*100 if health_cat == `h' & age_cat10 == `a' & college == `c' & eff_terc == `ee'
					qui sum prob_`disease' if health_cat == `h' & age_cat10 == `a' & college == `c' & l2.eff_terc == `ee'	
					replace dev_`disease'_lf = (r(mean) / mm_`disease' - 1)*100 if l2.eff_terc == `ee' & health_cat == `h' & age_cat10 == `a' & college == `c'
					qui sum prob_`disease' if health_cat == `h' & age_cat10 == `a' & college == `c' & eff_terc == `ee' & l2.eff_terc == `ee' 
					replace dev_`disease'_bf = (r(mean) / mm_`disease' - 1)*100 if eff_terc == `ee' & l2.eff_terc == `ee' & health_cat == `h' & age_cat10 == `a' & college == `c'
				}	
			}
		}
	}
}

* Save and export
mat tabA9 = J(9,4,.)
local cc = 1
foreach disease in diab canc hyp heart{
	forvalues ee = 1/3{
		qui {
		sum dev_`disease' if eff_terc == `ee'
		mat tabA9[`ee',`cc']=r(mean)
		sum dev_`disease'_lf if eff_terc == `ee'
		mat tabA9[`ee'+3,`cc']=r(mean)
		sum dev_`disease'_bf if eff_terc == `ee'
		mat tabA9[`ee'+6,`cc']=r(mean)		
		}
	}
	local cc = `cc'+1
}

matrix colnames tabA9 = "Diabetes" "Cancer" "Hypertension" "Heart_Condition" 
matrix rownames tabA9 = "Low" "Middle" "High" "Low" "Middle" "High" "Low" "Middle" "High"  

outtable using "$path\tables\TableA9", mat(tabA9) replace format (%9.2f) caption("Implied Probability of Disease by Past and Current Effort Tercile") nobox center