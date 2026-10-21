**********************************************************
******* Create Table A3
******

cd "$path\workdata"
use workdata04-18.dta, clear

*** For the details on the health effort components, see construct_healtheffort.do
mat effort_comp = J(3,4,.)
corr Lsport0 Nut0 InvNcig0

mat effort_comp=r(C)

pca Lsport0_resid Nut0_resid InvNcig0_resid
mat effort_comp = effort_comp,e(L)[1..3,1]

matrix colnames effort_comp = "Physical_Exercise" "Healthy_Nutrition" "Abstention_from_Smoking" "Loading"
matrix rownames effort_comp = "Physical_Exercise" "Healthy_Nutrition" "Abstention_from_Smoking" 

outtable using "$path\tables\TableA3", mat(effort_comp) replace format (%9.3f) caption("Health Effort Components and Weights") nobox center


