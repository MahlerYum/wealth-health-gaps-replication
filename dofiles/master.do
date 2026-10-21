******* 
*** Master Do-File for all Empirical Analysis done in Mahler & Yum (2024)
*** using data from SOEP-Core, v35, doi:10.5684/soep-core.v35 and Stata version 18.0

clear 

ssc install tabstatmat.pkg
ssc install outtable.pkg
ssc install inequal7.pkg

**********************************************************
******* O. Set Directories to wherever the Replication Folder is stored
**********************************************************

global path "...\Replication"

**********************************************************
******* 1. Create Working Data Set from Raw Data (only possible if raw data is available)
**********************************************************
do "$path\dofiles\create_workdata.do"


******* If working data (two data sets: workdata04-18.dta and survival_dat.dta) are constructed and stored in the workdata folder, comment line 20 and run from the top.


**********************************************************
******* 2. Get Tables and Figures (that require only Stata program)
**********************************************************
do "$path\dofiles\get_TableA5.do"
do "$path\dofiles\get_TableA3.do"
do "$path\dofiles\get_TableA2FiguresA1A2.do"
do "$path\dofiles\get_TableA14.do"
do "$path\dofiles\get_TableA6-A9.do"
do "$path\dofiles\get_Table1.do"
do "$path\dofiles\get_TableA4.do"
do "$path\dofiles\get_TableA1.do"

*********************************************************
******* 3. Export data used for figures and tables created in Matlab
**********************************************************
do "$path\dofiles\create_proddata.do"
do "$path\dofiles\export_data_figures.do"
do "$path\dofiles\export_data_tables.do"

*********************************************************
******* 4. Calculate conditional survival probabilities
**********************************************************
do "$path\dofiles\get_survivalrates.do"