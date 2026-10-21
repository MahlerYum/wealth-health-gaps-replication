! Module to declare parameters 
include 'mkl_vsl.fi'    
module mod_parameters   

use mkl_vsl

    implicit none

    integer, parameter :: estimation = 0               ! 1: data storage for standard error estimation
    integer, parameter :: recordsimstat = 0+1            ! 1: record simulation statistics for REPLICATION
    integer, parameter :: recordsimdata = 0            ! 1: record simulated data for counterfactuals
        
    ! Experiments    
    integer, parameter :: exp_prod = 0!+1      ! 1: Productivity at healthy levels
    integer, parameter :: exp_lab = 0!+1        ! 1: Labor disutility at healthy levels
    
    integer, parameter :: equalsurv = 0!+1      ! set to 1 for equalizing survival rates by health at GOOD health rates
    integer, parameter :: exp_no_consp = 0!+1 ! 1: no consumption utility penalty from bad health for everyone (includes quality of life)
    
    integer, parameter :: exp_hesav = 0     ! 1: savings decision averaged over health types
    integer, parameter :: exp_no_consp_ret = 0  ! 1: no consumption utility penalty from bad health only for old (includes quality of life)
    
    integer, parameter :: exp_wealthtax_pension_ret = 0!+1     ! 1: Tax savings for retired at tau_a, give uniform pension scale_tt * tt                                   
            integer, parameter :: tau_a = 1.0d0                     ! Wealth Tax Rate used above
            integer, parameter :: scale_tt = 13.426d0                   ! Scale for uniform pensions used above
    
    ! Baseline Distribution
    integer, parameter :: use_BM_health = 0!+1      ! 1: use baseline distribution over health state in simulation
    
    integer, parameter :: use_BM_assets = 0!+1      ! 1: use baseline distribution over savings choice in simulation
    integer, parameter :: use_BM_work = 0!+1        ! 1: use baseline distribution over labor supply choice in simulation
    
    integer, parameter :: use_BM_survival = 0    ! 1: use baseline distribution over survival state in simulation (needs to switched on if other BM states are used)
    integer, parameter :: use_BM_effort = 0      ! 1: use baseline distribution over effort choice in simulation
        
    
    ! ***********************************
    ! Fixed parameters 
    real(8), parameter :: sigma0 = 2.0d0            ! EIS for consumption
    real(8), parameter :: gamma = 1.0d0             ! Frisch Elasticity of Labor Supply
    real(8), parameter :: rho = 0.975d0             ! Persistence of productivity process     
    
    real(8), parameter :: r = 1.04d0**2.0d0 - 1.0d0         ! interest rate  
    
    real(8), parameter :: taul = 0.128d0              ! labor income tax progressivity parameter
    real(8), parameter :: lambda = 1.0d0 - 0.321d0    ! labor tax scale
    real(8), parameter :: mincon0 = 0.10d0            ! minimum consumption relative to average earnings
    real(8), parameter :: tt0 = 0.115d0               ! Sickness benefits: 70% replacement rate and maximum eligible periods 0.7 * 120 / (365*2)    
            
    ! Health technology estimates (NO EFFORT)
    real(8), parameter :: const_healthtr = 0.0129956d0
    real(8), parameter, dimension(6) :: age_const = (/0.0d0,-0.3007491d0,-0.6550076d0,-0.871416d0,-1.07407d0,-1.527431d0/)
    real(8), parameter, dimension(2) :: eff_param = (/0.0d0,0.0d0/)             ! Coefficient on current health and past EFFORT
    real(8), parameter :: eff_sq = 0d0
    real(8), parameter :: healthy_dummy = 2.340282d0
    real(8), parameter :: htype_dummy = 0.6543543d0
    real(8), parameter :: college_dummy = 0.3379621d0    
   
    
    ! ***************************************************************
    ! Parameters estimated internally
    ! 2026 Update: newly estimated with the corrected implementation
    ! ***************************************************************
    																																		
 	! Parameters for disutility of work    ! 
    ! healthy
    real(8), parameter :: nuh_10 = 2.66211197060956d0
    real(8), parameter :: nuh_20 = 1.7336221854522d0 
    real(8), parameter :: nuh_30 = 1.27921296006633d0
    real(8), parameter :: nuh_40 = 1.62184939862195d0
    
    ! unhealthy
    real(8), parameter :: nuu_10 = 2.80199514201671d0
    real(8), parameter :: nuu_20 = 1.80908438446131d0
    real(8), parameter :: nuu_30 = 1.38582547555823d0
    real(8), parameter :: nuu_40 = 1.9993696366332d0
        
    real(8), parameter :: nuadd0 = 0.829665748476618d0        
    
    real(8), parameter :: sigx0 = 0.0286816760967762d0        
    
    real(8), parameter :: penre0 = 0.366071123435765d0        
    
    real(8), parameter :: betamean0 = 0.958243726766528d0     
    real(8), parameter :: betadev0 = 0.0325187822742555d0     
    
    real(8), parameter :: bb0 = 12.6724977403706d0            
    
    real(8), parameter :: conp0 = 0.97784050879236d0          
    
    ! Wage profile for HS + healthy
    real(8), parameter :: ytHS_s0 = 0.0605668965222253d0           
    real(8), parameter :: ytHS_sq0 = -0.0024432648508462d0  
    
    ! Wage profile for CL + healthy
    real(8), parameter :: ytCL_s0 = 0.0840941030142781d0  
    real(8), parameter :: ytCL_sq0 = -0.00283222546167337d0 
    
    ! wage penalty: depends on education and age
    real(8), parameter :: wagep_HS0 = 0.202245800050694d0      ! 
    real(8), parameter :: wagep_CL0 = 0.129319898922546d0      !     
       
    
    ! Initial yt(1) for HS relative to CL
    real(8), parameter :: y1_HS0 = 1.24763388335393d0
    real(8), parameter :: y1_CL0 = 1.6103572846017d0																					
					
    

    
    ! **********************************************************************************************
    ! Technical settings
    real(8), parameter :: tol = 1d-4               
    real(8), parameter :: deggrida = 2.0d0
    integer, parameter :: diagnum = 10              ! display result in every diagnum interval (frequency)    
    integer, parameter :: nw = 3                       ! number of work choices    
    real(8), parameter :: fulltime = 1.0d0      
    
    ! Grids
    integer, parameter :: nt = 38                   ! 38 age periods (corresponding to ages 25-99)
    integer, parameter :: tr = 21                   ! First period of retirement (corresponding to ages 65)    
    integer, parameter :: tj = 11                   ! End of young age group
    
    integer, parameter :: na = 50                   ! Grid size of the asset grid
    integer, parameter :: nx = 5                    ! Grid size of the productivity grid    
    integer, parameter :: nh = 2                    ! Grid size of the health status grid
    
    integer, parameter :: ra = na-2
    
    real(8), parameter :: amax0 = 30.0d0              ! Highest Asset Position
    
    integer, parameter :: nz = nx            ! same length as permanent productivity
    integer, parameter :: ns = 200000       ! set size of simulated agents (Calculation of Gini coefficient gets to costly with higher numbers)
    integer, parameter :: nc = 20000         ! compressed dimension
    
    ! Fixed permanent heterogeneity (types); education x (beta x prod type x health type)
    integer, parameter :: ne = 2                     ! Education;            index ie
    integer, parameter :: nb = 2                     ! Discount factor;      index ib
    integer, parameter :: ntheta = 2                 ! Productivity type;    index itheta
    integer, parameter :: neta = 2                   ! Health type;          index ieta
    
    integer, parameter :: ntype = nb*ntheta*neta     ! # type excluding education  index itype
    
    real(8), parameter, dimension(ntheta) :: theta_val = (/exp(-0.2898d0),exp(0.2898d0)/)     

    real(8), parameter,dimension(2) :: winit =(/43978d0,48201d0/)    ! Among age 25-26 (healthy), wages adjusted to observed hours => winit(2) is normalized to be 1 in the model
    
    real(8), parameter :: avrgearn_dat = 57706.57d0             ! average 2-year earnings (including zeros) among working age population 
    

    !
    
    !***************************************
    ! Calibration
    
    ! Calibration option
    integer, parameter :: calind = 0         ! 0: No Calibration, 1: Calibration
    integer, parameter :: efweight = 1        ! 1: Using moment s.e. for weights, 0: Log difference with equal weights
    
    integer, parameter :: parm_est = 0        ! Index of parameter of efficient estimation for moment selection; 0 means no exclusions
    
    integer, parameter :: nparam = 42-19      ! # of internally calibrated parameters
    integer, parameter :: nmom = 64-28        ! # of targeted moments
    real(8), parameter :: ntarget = 82      ! number of (potential) target statistics: dimension of Datamoments text file
    integer, parameter :: nm = nmom
    
    integer, parameter :: calstep = 1  ! 
    
    integer, parameter :: deriv=1           !Parameter for the direction of numerical derivative: Right direction
    real(8), parameter :: gsnum = (3.0d0-sqrt(5.0d0))/2.0d0     ! Golden number
    real(8), parameter :: ddelta = tol*1d-3 ! delta for the secant method.
    
    ! ***************************************************************************************     

        
    real(8), parameter :: tdatamm(nmom)=(/0.6508581d0,0.7660204d0,0.8232445d0,0.6193264d0, &           ! 1-4. 2-year employment share among HEALTHY (4 moments) 
        0.5055072d0,0.5830671d0,0.6008949d0,0.4091998d0, &                                                  ! 5-8. 2-year employment share among UNHEALTHY (4 moments)       
        0.0619297d0,0.516081d0,1.165899d0,1.651459d0, 1.567324d0, 1.006182d0, &       ! 9-14. Median wealth age group 1-6, relative to norm=winit(2)
        1.237489d0, &          ! 15. Education gradients in employment            
        8.49264942390098d0, &                         ! 16. VSL multiple
        0.7456731d0, &                        ! 17.  Gini wealth
        1.163207d0, &                        ! 18. c(healthy)/c(unhealthy) for number of household = 1 or 2
                                                            ! 19-34. vector of earnings profiles in 1,000 EUR (4x4);
        35.39329d0, 49.37886d0, 55.95501d0, 42.21932d0,	&   ! HS-Healthy
        24.94774d0, 33.16593d0, 36.69067d0, 25.31111d0, &   ! HS-Unhealthy
        59.48338d0, 89.53806d0, 107.9282d0, 98.27698d0, &   ! CL-Healthy
        50.38816d0, 66.25301d0, 78.31755d0, 63.1325d0, &    ! CL-Unhealthy
        0.5952184d0, &                        ! 35. Variance of log earnings       
        0.4770515d0 /)                       ! 36. Pension replacement rate
        
        
        
        
    real(8), parameter :: tdatase(nmom)=(/0.0022079d0,0.001673d0,0.0015903d0,0.0024375d0, &            ! 1-4. 2-year employment share among HEALTHY (4 moments) 
        0.0078668d0,0.0054486d0,0.0045718d0,0.0045788d0, &                                                  ! 5-8. 2-year employment share among UNHEALTHY (4 moments)                
        0.0031501d0,0.0146831d0,0.023547d0,0.037393d0, 0.042682d0, 0.0473329d0, &      ! 9-14. Median wealth age group 1-6, relative to norm=winit(2)                   
        0.0029621d0, &          ! 15. Education gradients in employment        
        0.594830904063775d0, &                         ! 16. VSL multiple 
        0.0035907d0, &                        ! 17.  Gini wealth
        0.0221391d0, &                        ! 18. c(healthy)/c(unhealthy) for number of household = 1 or 2
                                                                ! 19-34. vector of earnings profiles in 1,000 EUR (4x4);
        0.1955369d0, 0.2318309d0, 0.2660378d0, 0.3528976d0,	&   ! HS-Healthy
        0.5630693d0, 0.5187444d0, 0.4988166d0, 0.4986972d0, &   ! HS-Unhealthy
        0.4875483d0, 0.631705d0, 0.7607303d0, 1.108492d0, &     ! CL-Healthy
        1.848723d0, 1.65571d0, 1.688008d0, 1.78551d0, &         ! CL-Unhealthy
        0.0023382d0, &                        ! 35. Variance of log earnings       
        0.0015815d0 /)                       ! 36. Pension replacement rate       
        
        
    
    !***************************************
    ! Parameters realted to Nelder-Mead method
    real(8), parameter :: rho_nm=1, chi_nm=2, gamma_nm=0.5, sigma_nm=0.5        
    Character(12), parameter :: calfileplace='.\cal_files\'
    
    ! Random number generation
    integer, parameter :: method = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
    Integer, parameter :: brng = vsl_brng_sfmt19937    
    Integer, parameter :: rnseed=68161

    
    
end module mod_parameters
