**********************************************************
*******  Relationship between Labor Earnings or Wealth and Consumption Expenditures on health-related goods and services (Table A1)
******

cd "$path\workdata"
use workdata04-18.dta, clear

* Correlation of Labor Income with Average household expenditure on health-related goods and services conditional on other observables
gen avr_hh_hcons = health_exp_mthly*12

reg avr_hh_hcons health_cat age yschool labinc, vce(robust)
mat tabA1 = J(12,2,.)
forvalues k=1/4{
	mat tabA1[`k'*2-1,1]=e(b)[1,`k']
	mat tabA1[`k'*2,1]=sqrt(e(V)[`k',`k'])
}
mat tabA1[11,1]=e(N)
mat tabA1[12,1]=e(r2)

mat tabA1[7,1]=tabA1[7,1]*1000
mat tabA1[8,1]=tabA1[8,1]*1000

* Since consumption is asked in 2010, but wealth only in 2007 and 2012, use the average of those two wealth points 
tsset pid syear
gen wealth_avr = (f2.net_wealth + l3.net_wealth)/2
reg avr_hh_hcons health_cat age yschool wealth_avr, vce(robust)
forvalues k=1/4{
	mat tabA1[`k'*2-1,2]=e(b)[1,`k']
	mat tabA1[`k'*2,2]=sqrt(e(V)[`k',`k'])
}
mat tabA1[11,2]=e(N)
mat tabA1[12,2]=e(r2)

mat tabA1[9,2]=tabA1[7,2]*1000
mat tabA1[10,2]=tabA1[8,2]*1000

mat tabA1[7,2]=.
mat tabA1[8,2]=.

matrix rownames tabA1 = "Good_Health" "SE" "Age" "SE" "Years_of_Education" "SE" "Earnings" "SE" "Wealth" "SE" "N" "R^2"

outtable using "$path\tables\TableA1", mat(tabA1) replace format (%9.2f) caption("Effect of Earnings and Wealth on Spending on Health Goods") nobox center

