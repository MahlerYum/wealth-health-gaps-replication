! Module to declare parameters 
include 'mkl_vsl.fi'    
module mod_parameters      

use mkl_vsl

    implicit none
    
    ! Parameters to estimate the model
    integer, parameter :: calind = 0          ! 0: No Calibration, 1: Calibration
    integer, parameter :: efweight = 1        ! 1: Using moment s.e. for weights, 0: Log difference with equal weights   
    integer, parameter :: parm_est = 0!        ! Index of parameter of efficient estimation for moment selection; 0 means no exclusions
    
    ! Parameters to switch on for experiments
    integer, parameter :: run_experiments = 1         ! 1: run all experiments for REPLICATION
    integer, parameter :: recordsimstat = 1          ! 1: record simulation statistics for REPLICATION
    
    ! Number of Experiment Sepcifications  and parameters governing experiments   
    integer, parameter :: nspecs = 16   
    
    integer, parameter :: exp_hesav = 0         ! 1: savings decision averaged over health types                                               
    integer, parameter :: exp_nowe_eff = 0      ! 1: Effort decision depends on average wealth and not actual individual wealth               
    
    ! 2026 Update: tau_a, scale_tt now real (integer declaration truncated scale_tt); scale_tt recomputed for re-estimated baseline (Table 4)
    real(8), parameter :: tau_a = 1.0d0                     ! Wealth Tax Rate used in Table 4 experiment
    real(8), parameter :: scale_tt = 14.914d0               ! Scale for uniform pensions used above

        
    
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
    
    ! Health technology estimates
    real(8), parameter :: const_healthtr = -0.906d0
    real(8), parameter, dimension(6) :: age_const = (/0.0d0,-0.289d0,-0.644d0,-0.881d0,-1.138d0,-1.586d0/)
    real(8), parameter, dimension(2) :: eff_param = (/0.693d0,0.734d0/)             ! Coefficient on current health and past EFFORT
    real(8), parameter :: eff_sq = 0d0
    real(8), parameter :: healthy_dummy = 2.311d0
    real(8), parameter :: htype_dummy = 0.632d0
    real(8), parameter :: college_dummy = 0.238d0
    
    
    ! ***************************************************************
    ! Parameters estimated internally
    ! 2026 Update: newly estimated with the corrected implementation
    ! ***************************************************************
    															
    ! Parameters for disutility of work    ! 
    ! healthy
    real(8), parameter :: nuh_10 = 2.69264906924743d0
    real(8), parameter :: nuh_20 = 1.73406945202999d0
    real(8), parameter :: nuh_30 = 1.30154786743289d0
    real(8), parameter :: nuh_40 = 1.76162848181229d0
    
    ! unhealthy
    real(8), parameter :: nuu_10 = 2.15838263070013d0
    real(8), parameter :: nuu_20 = 1.8575631558528d0
    real(8), parameter :: nuu_30 = 1.37414382944875d0
    real(8), parameter :: nuu_40 = 2.24998813045837d0
        
    real(8), parameter :: nuadd0 = 0.781909921869279d0         
    
    ! direct utility cost of effort
    ! HS-Healthy
    real(8), parameter :: xiHSh_10 = 0.323894011038064d0
    real(8), parameter :: xiHSh_20 = 1.1585447580704d0
    real(8), parameter :: xiHSh_30 = 2.01288238078007d0
    real(8), parameter :: xiHSh_40 = 2.4599232261209d0!
    
    
    ! HS-Unhealthy    
    real(8), parameter :: xiHSu_10 = 0.990641349451605d0
    real(8), parameter :: xiHSu_20 = 2.59084419824129d0
    real(8), parameter :: xiHSu_30 = 3.09226997127447d0
    real(8), parameter :: xiHSu_40 = 1.36098692014578d0
    																							
    
    ! CL-Healthy
    real(8), parameter :: xiCLh_10 = 0.172850459668474d0
    real(8), parameter :: xiCLh_20 = 0.68902057472667d0
    real(8), parameter :: xiCLh_30 = 1.54649052382899d0
    real(8), parameter :: xiCLh_40 = 2.47912783620706d0

    
    ! CL-Unhealthy    
    real(8), parameter :: xiCLu_10 = 0.667053672622597d0
    real(8), parameter :: xiCLu_20 = 1.87384591090647d0
    real(8), parameter :: xiCLu_30 = 3.00363658743826d0
    real(8), parameter :: xiCLu_40 = 2.10374118600681d0
    
    
    ! effort habit adjustment cost max
    real(8), parameter :: chi_10 = 0.000132144222460612d0          
    real(8), parameter :: chi_20 = 0.122525346794748d0              
    
    real(8), parameter :: sigx0 = 0.0291215886457768d0               
    
    real(8), parameter :: penre0 = 0.363223483122939d0           
    
    real(8), parameter :: betamean0 = 0.951374014899581d0 
    real(8), parameter :: betadev0 = 0.0347943968416017d0
    
    real(8), parameter :: bb0 = 12.8916624219592d0          
    
    real(8), parameter :: conp0 = 0.875420726634848d0       
    real(8), parameter :: psi0 = 1.04267952238969d0         
    
    ! Wage profile for HS + healthy
    real(8), parameter :: ytHS_s0 = 0.061276683941943d0           
    real(8), parameter :: ytHS_sq0 = -0.0025219478076732d0  
    
    ! Wage profile for CL + healthy
    real(8), parameter :: ytCL_s0 = 0.0871398363805312d0  
    real(8), parameter :: ytCL_sq0 = -0.00294739291224355d0 
    
    ! wage penalty: depends on education and age
    real(8), parameter :: wagep_HS0 = 0.202001661757094d0      ! 
    real(8), parameter :: wagep_CL0 = 0.153845746298704d0      !     
       
    
    ! Initial yt(1) for HS relative to CL
    real(8), parameter :: y1_HS0 = 1.25878322704252d0
    real(8), parameter :: y1_CL0 = 1.60305899945338d0
    
    
    
    ! **********************************************************************************************
    ! Technical settings
    real(8), parameter :: tol = 1d-4               
    real(8), parameter :: deggrida = 2.0d0
    integer, parameter :: diagnum = 10                 ! display result in every diagnum interval (frequency)           
    integer, parameter :: puregrid = 2                 ! 1: use pure grid method for effort choice; 2: use golden section search over effort    
    integer, parameter :: nfnow = 40                   ! Number of the effort choice (not used if puregrid = 2)    
    integer, parameter :: nw = 3                       ! number of work choices
    real(8), parameter :: fulltime = 1.0d0          
    integer, parameter :: nxi = 1                      ! number of effort disutility shocks (removed)
    
    ! Grids
    integer, parameter :: nt = 38                   ! 38 age periods (corresponding to ages 25-99)
    integer, parameter :: tr = 21                   ! First period of retirement (corresponding to ages 65)    
    integer, parameter :: tj = 11                   ! End of young age group
    
    integer, parameter :: na = 50                   ! Grid size of the asset grid
    integer, parameter :: nx = 5                    ! Grid size of the productivity grid
    integer, parameter :: nf = 30                   ! Grid size of the effort grid    
    integer, parameter :: nh = 2                    ! Grid size of the health status grid
    integer, parameter :: nchi = 200                 ! length of adjustment shock grid
    
    integer, parameter :: ra = na-2
    integer, parameter :: rf = nf-2
    
    real(8), parameter :: amax0 = 30.0d0              ! Highest Asset Position
    
    real(8), parameter :: fmin = 0.0d0                ! Minimal possible health effort
    real(8), parameter :: fmax = 1.0d0                ! Maximum possible health effort
    
    real(8), parameter :: chimin = 0.0d0
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
     
    
    !***************************************
    ! Calibration
    integer, parameter :: nparam = 42      ! # of internally calibrated parameters
    integer, parameter :: nmom = 64        ! # of targeted moments
    real(8), parameter :: ntarget = 82      ! number of (potential) target statistics: dimension of Datamoments text file
    integer, parameter :: nm = nmom
    
    integer, parameter :: calstep = 1  ! 
        
    integer, parameter :: deriv=1           !Parameter for the direction of numerical derivative: Right direction
    real(8), parameter :: gsnum = (3.0d0-sqrt(5.0d0))/2.0d0     ! Golden number
    real(8), parameter :: ddelta = tol*1d-3 ! delta for the secant method.
    
    ! ***************************************************************************************     

        
    real(8), parameter :: tdatamm(nmom)=(/0.6508581d0,0.7660204d0,0.8232445d0,0.6193264d0, &           ! 1-4. 2-year employment share among HEALTHY (4 moments) 
        0.5055072d0,0.5830671d0,0.6008949d0,0.4091998d0, &                                                  ! 5-8. 2-year employment share among UNHEALTHY (4 moments)       
    
        0.6777659d0,0.6769325d0,0.6802505d0,0.6992036d0,0.7301746d0,0.7237555d0, &  ! 9-14. Average Health Effort Profile of HS & HEALTHY by 10-year age groups (6 moments)
        0.6426084d0,0.6227545d0,0.627258d0,0.6552106d0,0.6968261d0,0.6921402d0, &  ! 15-20. Average Health Effort Profile of HS & UNHEALTHY by 10-year age groups (6 moments)
        0.7790819d0,0.7702285d0,0.7660254d0,0.7634262d0,0.779154d0,0.7724553d0, &  ! 21-26. Average Health Effort Profile of CL & HEALTHY by 10-year age groups (6 moments)
        0.7517721d0,0.7435739d0,0.736526d0,0.7381558d0,0.750504d0,0.734436d0, &  ! 27-32. Average Health Effort Profile of HS & UNHEALTHY by 10-year age groups (6 moments)    
        0.0619297d0,0.516081d0,1.165899d0,1.651459d0, 1.567324d0, 1.006182d0, &       ! 33-38. Median wealth age group 1-6, relative to norm=winit(2)
        1.237489d0, &          ! 39. Education gradients in employment        
        0.2672905d0,0.3283083d0,0.4041793d0, &        ! 40-42. non-adjuster shares by long age groups (with val_inaction = 0.0)           
        8.49264942390098d0, &                         ! 43. VSL multiple 
        0.1610319d0, &                        ! 44. Standard Deviation of Effort        
        0.7456731d0, &                        ! 45.  Gini wealth
        1.163207d0, &                        ! 46. c(healthy)/c(unhealthy) for number of household = 1 or 2
                                                            ! 47-62. vector of earnings profiles in 1,000 EUR (4x4);
        35.39329d0, 49.37886d0, 55.95501d0, 42.21932d0,	&   ! HS-Healthy
        24.94774d0, 33.16593d0, 36.69067d0, 25.31111d0, &   ! HS-Unhealthy
        59.48338d0, 89.53806d0, 107.9282d0, 98.27698d0, &   ! CL-Healthy
        50.38816d0, 66.25301d0, 78.31755d0, 63.1325d0, &    ! CL-Unhealthy
        0.5952184d0, &                        ! 63. Variance of log earnings       
        0.4770515d0 /)                       ! 64. Pension replacement rate          
        
        
        
        
    real(8), parameter :: tdatase(nmom)=(/0.0022079d0,0.001673d0,0.0015903d0,0.0024375d0, &            ! 1-4. 2-year employment share among HEALTHY (4 moments) 
        0.0078668d0,0.0054486d0,0.0045718d0,0.0045788d0, &                                                  ! 5-8. 2-year employment share among UNHEALTHY (4 moments)        
        0.0019615d0,0.0016137d0,0.0016517d0,0.0018318d0,0.0016836d0,0.0022494d0, &  ! 9-14. Average Health Effort Profile of HS & HEALTHY by 10-year age groups (6 moments)
        0.0066662d0,0.0047753d0,0.0035851d0,0.0031197d0,0.0027306d0,0.0025937d0, &  ! 15-20. Average Health Effort Profile of HS & UNHEALTHY by 10-year age groups (6 moments)
        0.0024741d0,0.0019636d0,0.0019423d0,0.0022411d0,0.0024561d0,0.0037815d0, &  ! 21-26. Average Health Effort Profile of CL & HEALTHY by 10-year age groups (6 moments)
        0.0107689d0,0.0082543d0,0.0063126d0,0.0051546d0,0.0050761d0,0.0057938d0, &  ! 27-32. Average Health Effort Profile of HS & UNHEALTHY by 10-year age groups (6 moments)         
        0.0031501d0,0.0146831d0,0.023547d0,0.037393d0, 0.042682d0, 0.0473329d0, &      ! 33-38. Median wealth age group 1-6, relative to norm=winit(2)                   
        0.0029621d0, &          ! 39. Education gradients in employment        
        0.0037247d0,0.0030799d0,0.0039969d0, &        ! 40-42. non-adjuster shares by long age groups (with val_inaction = 0.0)           
        0.594830904063775d0, &                         ! 43. VSL multiple
        0.0004399d0, &                        ! 44. Standard Deviation of Effort        
        0.0035907d0, &                        ! 45.  Gini wealth
        0.0221391d0, &                        ! 46. c(healthy)/c(unhealthy) for number of household = 1 or 2
                                                                ! 47-62. vector of earnings profiles in 1,000 EUR (4x4);
        0.1955369d0, 0.2318309d0, 0.2660378d0, 0.3528976d0,	&   ! HS-Healthy
        0.5630693d0, 0.5187444d0, 0.4988166d0, 0.4986972d0, &   ! HS-Unhealthy
        0.4875483d0, 0.631705d0, 0.7607303d0, 1.108492d0, &     ! CL-Healthy
        1.848723d0, 1.65571d0, 1.688008d0, 1.78551d0, &         ! CL-Unhealthy
        0.0023382d0, &                        ! 63. Variance of log earnings       
        0.0015815d0 /)                       ! 64. Pension replacement rate   
        
        
    
    !***************************************
    ! Parameters realted to Nelder-Mead method
    real(8), parameter :: rho_nm=1, chi_nm=2, gamma_nm=0.5, sigma_nm=0.5        
    Character(12), parameter :: calfileplace='.\cal_files\'
    
    ! Random number generation
    integer, parameter :: method = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
    Integer, parameter :: brng = vsl_brng_sfmt19937      
    Integer, parameter :: rnseed=68161

    
    
end module mod_parameters
