subroutine sub_main(x_param,y,evalnum)
    
    use mod_parameters
    use mod_utils
    use mod_solve
    use mod_spline
    use sort_all
    use mod_simul
    
    implicit none
    
    real(8), intent(in) :: x_param(nparam) 
    real(8), intent(out) :: y(nmom)
    integer, intent(inout) :: evalnum
    
    ! Define variables
    real(8) V_adjtemp, V_nadjtemp, V_w_temp, V_n_temp, V_n_nadj_temp, V_w_nadj_temp, trash, objval, trash1, trash2, trash3, valtemp, mucontemp
    real(8) agrid(na), xgrid(nx), zgrid(nz), pij(nx,nx), pz(nz,nz), bgrid(nb), zeta(nt,ne,nh)
    real(8) yt(tr-1,ne,nh), wageg(8,ne), ggrid(8), cfg(4,7,ne), yg(8,ne), phimat(tr-1,ne,nh)    
    real(8) pi0(6,nh,ne,neta), l1(6,nh,ne,neta), l2(6,nh,ne,neta), amin, amax, avrgearn, hours, Vvec_nadj(nw) 
    real(8) beta, anow, ap, ahigh, pbnow, spnow, pinow, phinow, awhigh, anhigh, ynow, ynow_net, eval(nh), aptemp, Vvaltemp, Mvaltemp
    real(8) tnugrid(4),  xiHShval(4), xiHSuval(4), xiCLhval(4), xiCLuval(4), nuhval(4), nuuval(4)    
    integer t, i, j, ix, jx, ia, ja, ik, jk, iz, jz, ib, jb, ih, jh, ie, ind, m, g, ind2, ind3, itheta, ieta, itype, indmat(3,ntype), maxind, iter
    real(8), allocatable :: Vmat(:,:,:,:,:,:), Vmatnext(:,:,:,:,:), V_nadj(:,:,:,:,:,:,:), ap_nadj(:,:,:,:,:,:,:)        
    
    integer :: st, ct1, ct2, ct0, ctf, cr, st1
    real(8) eptime
    
    ! Simulation
    integer :: is, count_active
    integer, allocatable :: apos(:,:), zpos(:,:), educ(:), thetapos(:), etapos(:), dpos(:), health(:,:), survival(:,:)
    real(8), allocatable :: assets(:,:), earn(:,:), savings(:), consum(:,:), bequests(:,:), work(:,:)    
    real(8) WtempS, NtempS, hy_contr, le_contr, fh_contr
    real(8) a, b, b1, b2, cumzss(nz), onedraw(1), cumpz(nz,nz)
    
    ! Calculation of Moments
    real(8) num, numa, numh, numwe, numf, numwo, numwo_part, numinc, norm, effcorr(3), weight_mat(nmom,nparam)
    real(8) avwealth2y(nt,nh), avcons2y(nt,nh), avwork2y(tr-1,nh), avearn2y(tr-1,nh),avrhealth2y(nt),avrhealthedu2y(nt,ne), aveff2yall(nt)
    real(8) wequart2y(nt,nh,3), wequart10y_av(7,3), hequ(7,4), earn_contr, avrg_lsdisutil, avrg_fdisutil,vsl2y
    real(8) avreffHS, avreffCL, avrempHS, avrempCL, num1, num2, num3, num4, shareW(5), negW, w2y_h, w2y_u, sum1
    integer :: denom1, denom2, denom3, denom4,  denom, denomall, denom2564, denomwo, inda, indinc, s, e, ihw, t1, t2
    
    real(8) numcons, avrcons25, le25, svl, totalwealth, totaly, w2y, penavg, aa, bbb, Sginiwe, Sginiy, Sginiy_HS, Sginiy_CL, we_contr    
    real(8), allocatable :: sagrid(:),sygrid(:),sadistr(:),sydistr(:),lsdisutil(:,:),fdisutil(:,:), Swevec(:), Sincvec(:), Swedist(:), Sincdist(:), wequart(:,:)    
    
    real(8) Swedistcomp(nc), Sweveccomp(nc), Sincdistcomp(nc), Sincveccomp(nc), Swedistcomp2(nc), Sweveccomp2(nc)
    real(8) long_avr_earn(2,nh), long_avr_work(2,nh), long_avr_work_part(2,nh), long_avr_cons(3,nh)    
    
    real(8) avrgearn_sim, avrcons, avrriskdeath, vsl, avreffAll, hwgrid(nw), avrg_consutil, avrg_marg_consutil, avgcons, avghours, avgf, avgkappa, avgphi, avgxi    
    real(8) wequart10y(7,nh,3), long_avr_earn10(4,nh), long_avr_work10(4,nh), long_avr_work_part10(4,nh), long_avr_cons10(7,nh), long_avr_health10(7,ne), &
        long_avr_earn10_edu(4,ne), long_avr_work10_edu(4,ne), long_eff_noadj(3), long_eff_pos(3), long_eff_neg(3), &
        eff_pos, eff_neg
    real(8) long_avr_earn10_health_edu(4,2,2), long_avr_work10_health_edu(4,2,2)
    real(8) long_avr_wealth10_ter(7,2,3), avr_cons_work(nh), avr_wealth10_ter4(3), wequart10y_edu_health(7,nh,ne,3), wequart10y_edu(7,ne,3)
    real(8) wequart20y_edu(3,ne,3)
    
    type(VSL_STREAM_STATE) :: stream

    ! Calibration
    !
    real(8) phi1_i, phi1_s, phi1_sq, phi2_i, phi2_s, phi2_sq, xi1_i, xi1_s, xi1_sq, xi2_i, xi2_s, xi2_sq, mincon, y1_HS, y1_CL
    real(8) wagep_HS, wagep_CL, conp, betamean, betadev, psi, penre, beta1, beta2, nu1_i, nu1_s, nu1_sq, nu2_i,nu2_s,nu2_sq, xiadd, xiadd_1, xiadd_2, xiadd_3, xiadd_4, xiadd_5, xiadd_6, nuadd, tt, sigma
    real(8) xiHSh_1, xiHSh_2, xiHSh_3, xiHSh_4, xiHSu_1, xiHSu_2, xiHSu_3, xiHSu_4, xiCLh_1, xiCLh_2, xiCLh_3, xiCLh_4, xiCLu_1, xiCLu_2, xiCLu_3, xiCLu_4, chi_1, chi_2, bbval, chi_3, sdxi, sdztemp
    real(8) ytHS_s, ytHS_sq, ytCL_s, ytCL_sq, sigx, nuh_1, nuh_2, nuh_3, nuh_4, nuu_1, nuu_2, nuu_3, nuu_4
    character(len=100) :: filename
    
    ! spline
    real(8), allocatable :: EV(:,:), Tnumat(:,:), cfnuh(:,:), cfnuu(:,:)
    real(8), allocatable :: agedeath(:), pensionvec(:), marg_consutil(:,:), phisim(:,:)
    
	allocate(EV(na,nh), Tnumat(2,2), cfnuh(4,2), cfnuu(4,2))           
    allocate(agedeath(ns), pensionvec(ns), marg_consutil(ns,nt),phisim(ns,tr-1))
    
    ! Define the parameters that are evaluated
    
    nuh_1 = x_param(1)
    nuh_2 = x_param(2)    
    nuh_3 = x_param(3)
    nuh_4 = x_param(4)
    nuu_1 = x_param(5)
    nuu_2 = x_param(6)
    nuu_3 = x_param(7)
    nuu_4 = x_param(8)
    
    nuadd = x_param(9)
         
    sigx = x_param(10)              
    penre = x_param(11)
    
    betamean = x_param(12)
    betadev = x_param(13)
        
    bbval = x_param(14)
    conp = x_param(15)    
    
    ytHS_s = x_param(16)
    ytHS_sq = x_param(17)
    ytCL_s = x_param(18)
    ytCL_sq = x_param(19)
    
    wagep_HS = x_param(20)
    wagep_CL = x_param(21)
    
    y1_HS = x_param(22)
    y1_CL = x_param(23)
    
    sigma = sigma0
           
        
    sdxi= 0.0d0
    
    if (exp_no_consp == 1) then
        conp = 1.0d0
    endif
            
    !**************************** Preliminaries  *******************************
    
    call system_clock(ct0, cr)
    
    ! Get survival probabilites
    ! non-college
    open(12, file="surv_HS.txt")
    do i=1,nt
        read(12,*) zeta(i,1,:)
    enddo
    close(12)
    open(12, file="surv_CL.txt")
    do i=1,nt
        read(12,*) zeta(i,2,:)
    enddo
    close(12)
    
    
    ! Load Weight Matrix
    open(18, file="Weights.txt")
    do i = 1,nmom
        read(18,*) weight_mat(i,:)
    enddo
    close(18)
    
    
    ! Hours grid: 0 to full-time with number of positive choices being nw-1
    hwgrid = linspace(0.0d0,fulltime,nw)
    
    ! Grids for labor productivity
    ! Approximate permanent component using rouwenhorst method
    ! nz = nx 
    call rouwenhorst(rho,sqrt(sigx),nx,xgrid,pij)     
    zgrid = exp(xgrid)
    pz = pij
       
    if (exp_prod == 1) then        
        wagep_HS = 0.0d0
        wagep_CL = 0.0d0    
    end if
    
    ! winit(2) is normalized to be 1 in the model
    do t=1,tr-1        
        ! HS + healthy       
        yt(t,1,1) = y1_HS*exp( ytHS_s*(dble(t)-1.0d0) + ytHS_sq*(dble(t)-1.0d0)**2.0d0 )
        ! HS + unhealthy        
        yt(t,1,2) = yt(t,1,1)*(1.0d0-wagep_HS)
        
        ! CL + healthy                
        yt(t,2,1) = y1_CL*exp( ytCL_s*(dble(t)-1.0d0) + ytCL_sq*(dble(t)-1.0d0)**2.0d0 )
        ! CL + unhealthy
        yt(t,2,2) = yt(t,2,1)*(1.0d0-wagep_CL)
        
        if (t==tr-1) then
            write(*,"('  wagep (HS) =', F10.5)") wagep_HS
	        write(*,"('  wagep (CL) =', F10.5)") wagep_CL
        end if
               
    end do
    
    ! adjustment
    ! 2026 Update: normalize wages by E[exp(theta)]E[exp(z)]; sigx is a variance
    sdztemp = (sigx/(1.0d0-rho**2.0d0))**0.5d0
    yt = yt/(exp((log(theta_val(2))**2.0d0)/2.0d0)*exp((sdztemp**2.0d0)/2.0d0))
    
    nuhval = (/nuh_1, nuh_2, nuh_3, nuh_4/)
    nuuval = (/nuu_1, nuu_2, nuu_3, nuu_4/)
    
    if (exp_lab == 1) then
        nuuval = nuhval        
    end if
     
    tnugrid = (/1.0d0,8.0d0,13.0d0,20.0d0/) 
            
    ! setup splines    
    Tnumat = spbas(4-2,tnugrid)
    Tnumat = invmat(Tnumat)
    
    ! fit spline    
    cfnuh = spfit(Tnumat,nuhval,4-2,tnugrid)
    cfnuu = spfit(Tnumat,nuuval,4-2,tnugrid)
    
    do t=1,tr-1
        do ie=1,ne
            do ih=1,nh                
                if (ie==1) then
                    ! HS
                    if (ih==1) then
                        ! Healthy
                        call speva(cfnuh,dble(t),4-2,tnugrid,valtemp,trash,trash)
                        
                        phimat(t,ie,ih) = max(0.0d0,valtemp) *exp(nuadd)
                    else
                        ! Unhealthy
                        call speva(cfnuu,dble(t),4-2,tnugrid,valtemp,trash,trash)                
                        
                        phimat(t,ie,ih) = max(0.0d0,valtemp) *exp(nuadd)
                    end if
                else
                    ! CL
                    if (ih==1) then
                        ! Healthy
                        call speva(cfnuh,dble(t),4-2,tnugrid,valtemp,trash,trash)                
                        phimat(t,ie,ih) = max(0.0d0,valtemp)
                    else
                        ! Unhealthy
                        call speva(cfnuu,dble(t),4-2,tnugrid,valtemp,trash,trash)                
                        phimat(t,ie,ih) = max(0.0d0,valtemp)
                    end if
                end if
            end do
        end do
    end do
    
    
    write(*, '(x)')
    write(*, '(a)') ' Deterministic wage profiles'
    write(*, '(a)') ' Age           H       U'
    write(*, '(a)') '           HS        CL       HS       CL'
    do t=1,tr-1
        write(*,'(i2,a,4f10.4)') t, ' | ', yt(t,1,1), yt(t,2,1), yt(t,1,2), yt(t,2,2)         
    end do
    
    avrgearn = avrgearn_dat/winit(2) 
    tt = tt0*avrgearn     
    mincon = mincon0*avrgearn

    ! Define asset grid
    amin = 0.0d0
    agrid = loggrid(amin,amax0,na,deggrida)
    
    
    ! Define preference type grid
    bgrid(1) = betamean-betadev
    bgrid(2) = betamean+betadev
    
    write(*,"('  bgrid =', F10.5, F10.5)") bgrid(1), bgrid(2)
       
   
    ! adjusted to fit our order of states (nb,ntheta,neta)
    indmat = dimreduce3_new(nb,ntheta,neta)
    
    ! Allocate Space        
    allocate(Vmat(na,nz,ne,ntype,nh,nt), Vmatnext(na,nz,ne,ntype,nh))  
    allocate(V_nadj(na,nz,ne,ntype,nh,nw,nt))
    allocate(ap_nadj(na,nz,ne,ntype,nh,nw,nt))
    
    call system_clock(ct1, cr)
    
    ! ----------------------------------------------------------------------------------------- !
    ! Backward induction
    do t = nt,1,-1 
    
        if (t<nt) then
            Vmatnext = Vmat(:,:,:,:,:,t+1)                    
        else
            Vmatnext = 0.0d0
        end if
        
        if (mod(t,diagnum)==0 ) then
            write(*,'(a,x,i2)')'Solving value functions and decision rules at age',t
        end if        
        
        if (t >= tr) then       
            ! Retirement age group only face uncertainty regarding health
            if (t >= 26) then
                j = 6               
            else
                j = 5
            endif
            
            !$omp parallel do default(shared) private(iz,ie,ib,itheta,ieta,itype,maxind,jh,ih,EV,ia,beta,anow,pbnow,phinow,spnow,ahigh,pinow,aptemp,mucontemp)                
            
            do itype = 1,ntype
                do ie = 1,ne
                    do iz = 1,nz
                        do jh = 1,nh                                                                                                
                            EV(:,jh) = Vmatnext(:,iz,ie,itype,jh)          ! No need to take expectations over productivity realizations
                            
                        end do
                        
                        do ih =1,nh    
                            do ia = 1,na
                                    
                                if (ih==1) then                                            
                                    mucontemp = 1.0d0
                                else     
                                    if (exp_no_consp_ret == 1) then
                                        mucontemp = 1.0d0
                                    else
                                        mucontemp = conp
                                    endif
                                            
                                end if
                                    
                                ib = indmat(1,itype)
                                itheta = indmat(2,itype)
                                ieta = indmat(3,itype)
                                        
                                beta = bgrid(ib)
                                anow = agrid(ia)
                                if (exp_wealthtax_pension_ret == 1) then
                                    pbnow = scale_tt*tt             
                                else
                                    pbnow = pension(ie,iz,ih,yt(tr-1,:,:)*theta_val(itheta),zgrid,penre)                               
                                end if
                                        
                                phinow = 0.0d0  
                                spnow = zeta(t,ie,ih)            
                                       
                                if (exp_wealthtax_pension_ret == 1) then
                                    ahigh = (1.0d0+r)*anow - tau_a*max(anow,0d0) + pbnow
                                else
                                    ahigh = (1.0d0+r)*anow + pbnow 
                                end if                                    
                                
                                ! ************************
                                ! (1) Value of NOT ADJUST                                        
                                pinow = tr2yp(j,ie,ih,ieta)                                        
                                                                                
                                call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV,beta,mucontemp,bbval,sigma,spnow,pinow,phinow,psi,ih,t,ie,hwgrid(1),mincon)   ! not working                                    
                                V_nadj(ia,iz,ie,itype,ih,:,t) = -objneg_w(aptemp,ahigh,agrid,EV,beta,mucontemp,bbval,sigma,spnow,pinow,phinow,psi,ih,t,ie,hwgrid(1),mincon)  ! not working
                                ap_nadj(ia,iz,ie,itype,ih,:,t) = aptemp
                                    
                                if (ahigh - aptemp .le. 0d0) then
                                    write(*,'(a)') ' Negative Consumption: Ret, no-adjust'
                                endif
                                        
                                    
                            end do   
                        end do   
                                                
                    end do   
                end do   
            end do               
            !$omp end parallel do
            
            do ih = 1,nh
                do itype = 1,ntype
                    do ie = 1,ne
                        do iz = 1,nz
                            do ia = 1,na
                                    
                                Vmat(ia,iz,ie,itype,ih,t) = V_nadj(ia,iz,ie,itype,ih,1,t)                                    
                                    
                            end do
                        end do
                    end do
                end do
            end do
            
        else
            ! =============================     
            ! Decisions in Work Life
            ! =============================     
            
            if (t>=16) then
                j = 4
            elseif (t>=11) then
                j = 3
            elseif (t>=6) then
                j = 2
            else
                j = 1
            end if
            
            !$omp parallel do default(shared) private(iz,ie,ib,itheta,ieta,itype,maxind,jh,ih,EV,ia,beta,anow,pbnow,phinow,spnow,pinow,aptemp,ynow,ynow_net,anhigh,awhigh,mucontemp,hours,ihw)
            do itype = 1,ntype
                do ie = 1,ne
                    do iz = 1,nz
                        do jh = 1,nh
                            if (t == tr-1) then                 ! No expectations over productivity prior to retirement 
                                EV(:,jh) = Vmatnext(:,iz,ie,itype,jh)                                                  
                                
                            else                                ! Compute conditional expected value functions on grid points
               
                                EV(:,jh) = 0.0d0                
                                do ja = 1,na
                                    do jz = 1,nz
                                        EV(ja,jh) = EV(ja,jh) + pz(iz,jz)*Vmatnext(ja,jz,ie,itype,jh)
                                    end do
                                end do
                                
                            end if 
                                
                        end do
                            
                        do ih=1,nh
                            do ia=1,na    
                                if (ih==1) then
                                    mucontemp = 1.0d0
                                else                                            
                                    mucontemp = conp
                                end if
                                        
                                ib = indmat(1,itype)
                                itheta = indmat(2,itype)
                                ieta = indmat(3,itype)
                                    
                                beta = bgrid(ib)
                                anow = agrid(ia)
                                spnow = zeta(t,ie,ih)                              
                                phinow = phimat(t,ie,ih)                                                              
                                    
                                do ihw = 1,nw
                                    hours = hwgrid(ihw)
                                    ynow = yt(t,ie,ih)*zgrid(iz)*theta_val(itheta)*hours                
                                            
                                    ynow_net = ynow - tax(ynow,avrgearn,taul,lambda)       
                                            
                                    ! highest possible savings if working            
                                    if ((hours==0.0d0).and.(ih==2)) then
                                        awhigh = (1.0d0+r)*anow + ynow_net + tt                                            
                                    else
                                        awhigh = (1.0d0+r)*anow + ynow_net 
                                    end if    
                                                                                  
                                    ! ************************ 
                                    ! (1) Values of NOT ADJUST                                            
                                    pinow = tr2yp(j,ie,ih,ieta)                                            
                                            
                                    call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV,beta,mucontemp,bbval,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)  
                                    V_nadj(ia,iz,ie,itype,ih,ihw,t) = -objneg_w(aptemp,awhigh,agrid,EV,beta,mucontemp,bbval,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)   
                                            
                                    ap_nadj(ia,iz,ie,itype,ih,ihw,t) = aptemp
                                        
                                end do
                                    
                            end do
                        end do                          
                    end do   
                end do 
            end do 
            !$omp end parallel do
        
            do ih = 1,nh
                do itype = 1,ntype
                    do ie = 1,ne
                        do iz = 1,nz
                            do ia = 1,na
                                Vvec_nadj = 0.0d0
                                            
                                do ihw = 1,nw                                                
                                    Vvec_nadj(ihw) = V_nadj(ia,iz,ie,itype,ih,ihw,t)                                             
                                end do
                                    
                                Vmat(ia,iz,ie,itype,ih,t) = maxval(Vvec_nadj)    
                                    
                            end do
                        end do
                    end do
                end do
            end do
                    
        end if   ! j
    
    end do   ! t
                     

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    
    write(*,*) '    Values done'
    write(*,"('  Elasped time = ', F10.5)") eptime
   
    ! ************************************* SIMULATION ***************************************************
    call system_clock(ct1, cr)    
    
    allocate(educ(ns), dpos(ns), thetapos(ns), etapos(ns), health(ns,nt), work(ns,tr-1), survival(ns,nt))
    allocate(assets(ns,nt), earn(ns,tr-1), savings(ns), consum(ns,nt), bequests(ns,nt), zpos(ns,nt))  

    call simulation(V_nadj, ap_nadj, hwgrid, agrid, zeta, yt, zgrid, pz, penre, &
        assets, earn, consum, bequests, agedeath, work, educ, dpos, thetapos, etapos, zpos, health, survival, pensionvec, indmat, tt, mincon, avrgearn)  

    if (recordsimdata==1) then
        ! Save simulated data for counterfactual exercises
    
        open(30, file="healthBM.txt")
        open(31, file="assetsBM.txt")
        open(34, file="workBM.txt")
        open(38, file="survivalBM.txt")
        
        do t=1,nt
            do is=1,ns
                if (t<tr) then
                    write(34, *) work(is,t)         
                end if
                write(30, *) health(is,t)
                write(31, *) assets(is,t)
                write(38, *) survival(is,t)
                
            end do
        end do  
    
        close(30)
        close(31)
        close(34)
        close(38)
    
    end if
    
    write(*,*) ' Calculating Model Statistics'
    ! ********************************* Calculate Model Statistics ******************************

    ! Average Health by education and period
    denom = 0
    do t=1,nt
            denom = sum(survival(:,t))
            avrhealth2y(t) =  dble(sum(survival(:,t), mask=health(:,t)==2)) / denom 
        do ie=1,ne
            denom = sum(survival(:,t), mask=educ(:)==ie)
            avrhealthedu2y(t,ie) = dble(sum(survival(:,t), mask=health(:,t)==2 .AND. educ(:)==ie)) / denom
        enddo
    enddo
    
    denom = 0.0d0
    sum1 = 0.0d0
    do t=1,nt
        sum1 = sum1 + sum(survival(:,t), mask=health(:,t)==2 .AND. educ==1)
        denom = denom + sum(survival(:,t), mask=educ==1)
    enddo
    write(*, '(x)')
    write(*, '(a,f10.4)') ' Average Share Unhealthy (HS) ', sum1/denom
    
    denom = 0.0d0
    sum1 = 0.0d0
    do t=1,nt
        sum1 = sum1 + sum(survival(:,t), mask=health(:,t)==2 .AND. educ==2)
        denom = denom + sum(survival(:,t), mask=educ==2)
    enddo
    write(*, '(x)')
    write(*, '(a,f10.4)') ' Average Share Unhealthy (CL) ', sum1/denom
    
    !write(*, '(x)')
    !write(*, '(a)') ' Simulation Shares in Bad Health'
    !write(*, '(a)') ' Age   Overall  HS  CL'
    !do t=1,nt
    !    write(*, '(x,i2,a,3f10.4)') t, '     |  ', avrhealth2y(t), avrhealthedu2y(t,1), avrhealthedu2y(t,2)
    !end do
    
    ! Average labor income (including zeros) to benchmark model values
    avrgearn_sim = sum(earn, mask=survival(:,1:tr-1)==1) / sum(survival(:,1:tr-1))
    
    norm = winit(2) ! 30.07.2022 winit(2) is normalized to be 1 in the model
    
    
    ! Average of wealth, consumption, earnings, work by health and 2-year period (only survivors)
    do ih = 1,nh
        do t=1,nt
    
            denom = sum(survival(:,t), mask=health(:,t)==ih)
            avwealth2y(t,ih) = sum(assets(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom
            avcons2y(t,ih) = sum(consum(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom
            
            if (t<tr) then
                denomwo = sum(work(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih)   
                avwork2y(t,ih) = dble(denomwo) / denom            
                avearn2y(t,ih) = sum(earn(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom      ! avr earn across health type including zeros 31.08.2021
            end if
        
        end do
    end do

    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages by period of HEALHTY (in thousand euros)'    
    write(*, '(a)') ' Age Group  Avr.Wealth  Avr.Cons  Empl.Rate  Avr.Earnings (all)'
    do t=1,tr-1
        write(*, '(x,i2,a,6f10.4)') t, '     |  ', avwealth2y(t,1)*norm/1000, avcons2y(t,1)*norm/1000, avwork2y(t,1), avearn2y(t,1)*norm/1000
    end do
    do t=tr,nt
        write(*, '(x,i2,a,3f10.4)') t, '     |  ', avwealth2y(t,1)*norm/1000, avcons2y(t,1)*norm/1000
    end do
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages by period of UNHEALHTY (in thousand euros)'    
    write(*, '(a)') ' Age Group  Avr.Wealth  Avr.Cons  Empl.Share  Avr.Earnings (all)'
    do t=1,tr-1        
        write(*, '(x,i2,a,6f10.4)') t, '     |  ', avwealth2y(t,2)*norm/1000, avcons2y(t,2)*norm/1000, avwork2y(t,2), avearn2y(t,2)*norm/1000
    end do
    do t=tr,nt
        write(*, '(x,i2,a,3f10.4)') t, '     |  ', avwealth2y(t,2)*norm/1000, avcons2y(t,2)*norm/1000
    end do
    
    ! Calculate share of people alive in each period
    write(*,'(x)')
    write(*,'(a)') ' Share of people alive by Age'
    do t=1,nt
        write(*, '(x,i2,a,f10.4)') t, '     |  ', sum(survival(:,t)) / dble(ns)
    enddo
    
    
    ! 15.09.2021 Average Effort and Employment by Education    
    num3 = 0.0d0
    num4 = 0.0d0
    denom3 = 0
    denom4 = 0
    do is = 1,ns
        if (educ(is) == 1) then
            
            do t=1,(tr-1)
                if (survival(is,t) == 1) then
                num3 = num3 + work(is,t)
                denom3 = denom3 + 1
                endif
            enddo
            
        else
            
            do t=1,(tr-1)
                if (survival(is,t) == 1) then
                num4 = num4 + work(is,t)
                denom4 = denom4 + 1
                endif
            enddo
            
        end if
    enddo
    
    avrempHS = num3 / dble(denom3)
    avrempCL = num4 / dble(denom4)
    
    write(*, '(x)')
    write(*, '(x,a)') '                     Model  |   Data'

    write(*, '(x,a,2f10.4)') ' Avrg. Emp among HS ', avrempHS, 0.6537d0
    write(*, '(x,a,2f10.4)') ' Avrg. Emp among CL', avrempCL, 0.8089d0
    
    
    ! -------------                   
    
    do j= 1,7
        if (j==1) then
            s = 1
        elseif (j==2) then
            s = 6
        elseif (j==3) then
            s = 11
        elseif (j==4) then
            s = 16
        elseif (j==5) then
            s = 21
        elseif (j==6) then
            s = 26
        else
            s = 31
        end if
    
        if (j<7) then
            e = s+4
        else
            e = nt
        endif
        
        do ih = 1,nh
            if (j < 5) then
               long_avr_earn10(j,ih) = sum(earn(:,s:e), mask=survival(:,s:e)==1 .AND. health(:,s:e)==ih) / dble(sum(survival(:,s:e), mask=health(:,s:e)==ih))     ! Average earnings across everyone from that health group
               long_avr_work10(j,ih) = sum(work(:,s:e), mask=health(:,s:e)==ih .AND. survival(:,s:e)==1) / dble(sum(survival(:,s:e), mask=health(:,s:e)==ih))     
               long_avr_work_part10(j,ih) = sum(survival(:,s:e), mask=(health(:,s:e)==ih) .AND. (work(:,s:e)<=0.5)) / dble(sum(survival(:,s:e), mask=(health(:,s:e)==ih) .and. (work(:,s:e)>0.0d0)))
            end if   
            long_avr_cons10(j,ih) = sum(consum(:,s:e), mask=survival(:,s:e)==1 .AND. health(:,s:e)==ih) / sum(survival(:,s:e), mask=health(:,s:e)==ih)  
                        
           ! 13.07.2022 Calculate average wealth by wealth terciles and health
            denom = sum(survival(:,s:e), mask=health(:,s:e)==ih)
            if(allocated(sagrid))deallocate(sagrid)
            if(allocated(sadistr))deallocate(sadistr)
            allocate(sagrid(denom),sadistr(denom))
            sadistr = 1.0d0
            inda=1
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. health(is,t) == ih) then
                        sagrid(inda)=assets(is,t)
                        inda=inda+1
                    endif 
                enddo
            enddo
            call compressdist(sagrid,denom,nc,Swedistcomp,Sweveccomp)   
            call quart(Swedistcomp,Sweveccomp,nc,wequart10y(j,ih,:))          
            ! Get age and health specific wealth terciles
            t1 = ceiling(denom/3d0)
            t2 = ceiling(denom*(2d0/3d0))            
            
            ! Get average wealth by terciles
            long_avr_wealth10_ter(j,ih,1) = sum(sagrid(1:t1-1)) / dble(denom/3)
            long_avr_wealth10_ter(j,ih,2) = sum(sagrid(t1:t2-1)) / dble(denom/3)
            long_avr_wealth10_ter(j,ih,3) = sum(sagrid(t2:denom)) / dble(denom/3)
            
        end do
        
        ! Calculate unconditional average wealth by wealth terciles for age group 4
        !if (j==4) then
        if (j==4 .and. calind==0) then
            
            denom = sum(survival(:,s:e))
            if(allocated(sagrid))deallocate(sagrid)
            if(allocated(sadistr))deallocate(sadistr)
            allocate(sagrid(denom),sadistr(denom))
            sadistr = 1.0d0
            inda=1
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1) then
                        sagrid(inda)=assets(is,t)
                        inda=inda+1
                    endif 
                enddo
            enddo
            call quicksort(sagrid)
            ! Get age specific wealth terciles
            t1 = ceiling(denom/3d0)
            t2 = ceiling(denom*(2d0/3d0))            
            
            ! Get average wealth by terciles
            avr_wealth10_ter4(1) = sum(sagrid(1:t1-1)) / dble(denom/3)
            avr_wealth10_ter4(2) = sum(sagrid(t1:t2-1)) / dble(denom/3)
            avr_wealth10_ter4(3) = sum(sagrid(t2:denom)) / dble(denom/3)
           
            
        endif
        
        
        ! Calculate health (and effort) by wealth quartiles
        denom = sum(survival(:,s:e))
        if(allocated(sagrid))deallocate(sagrid)
        if(allocated(sadistr))deallocate(sadistr)
        allocate(sagrid(denom),sadistr(denom))
        sadistr = 1.0d0
        inda=1
        do t=s,e
            do is=1,ns
                if (survival(is,t) == 1) then
                    sagrid(inda)=assets(is,t)
                    inda=inda+1
                endif 
            enddo
        enddo
        call compressdist(sagrid,denom,nc,Swedistcomp,Sweveccomp)    
        call quart(Swedistcomp,Sweveccomp,nc,wequart10y_av(j,:))
        
        do t=s,e
            do is=1,ns
                if (survival(is,t) == 1) then
                    if (assets(is,t) <= wequart10y_av(j,1)) then
                        wequart(is,t) = 1
                    elseif ((assets(is,t) > wequart10y_av(j,1)) .and. (assets(is,t) <= wequart10y_av(j,2))) then
                        wequart(is,t) = 2
                    elseif ((assets(is,t) > wequart10y_av(j,2)) .and. (assets(is,t) <= wequart10y_av(j,3))) then
                        wequart(is,t) = 3
                    else
                        wequart(is,t) = 4
                    endif
                endif
            enddo
        enddo
        
        do i = 1,4
            hequ(j,i) = sum(survival(:,s:e), mask=(health(:,s:e)==2) .AND. (wequart(:,s:e)==i)) / dble(sum(survival(:,s:e), mask=(wequart(:,s:e)==i)))            
        enddo
        
        ! --- Moments by Education
        do ie=1,ne
            num = 0
            denom = 0
            numinc = 0
            numwo = 0
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. educ(is) == ie) then
                       denom = denom +1
                        if(health(is,t) == 2) then
                            num = num + 1
                        endif
                        if (j<5) then
                            numinc = numinc + earn(is,t)
                            numwo = numwo + work(is,t)
                        endif
                        !numf = numf + effort(is,t) 
                    endif
                enddo
            enddo
            
            long_avr_health10(j,ie) = dble(num) / dble(denom)
            if (j<5) then 
                long_avr_earn10_edu(j,ie) = dble(numinc) / dble(denom)
                long_avr_work10_edu(j,ie) = dble(numwo) / dble(denom)
            end if   
            
            ! Wealth Quartiles by Education
            if(allocated(sagrid))deallocate(sagrid)
            if(allocated(sadistr))deallocate(sadistr)
            allocate(sagrid(denom),sadistr(denom))
            sadistr = 1.0d0
            inda=1
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. educ(is) == ie) then
                        sagrid(inda)=assets(is,t)
                        inda=inda+1
                    endif 
                enddo
            enddo
            call compressdist(sagrid,denom,nc,Swedistcomp,Sweveccomp)     
            call quart(Swedistcomp,Sweveccomp,nc,wequart10y_edu(j,ie,:))            

        end do
        
        
        !---- Moments by Education and Health
        do ie=1,ne
            do ih=1,nh
            num = 0
            denom = 0
            numinc = 0
            numwo = 0
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. educ(is) == ie .and. health(is,t) == ih) then
                       denom = denom +1
                        if (j<5) then
                        numinc = numinc + earn(is,t)
                        numwo = numwo + work(is,t)
                        endif
                    endif
                enddo
            enddo
            
            if (j<5) then 
                long_avr_earn10_health_edu(j,ih,ie) = dble(numinc) / dble(denom)
                long_avr_work10_health_edu(j,ih,ie) = dble(numwo) / dble(denom)
            end if   

            ! Calculate Wealth Quartiles by Health and Education
            if(allocated(sagrid))deallocate(sagrid)
            if(allocated(sadistr))deallocate(sadistr)
            allocate(sagrid(denom),sadistr(denom))
            sadistr = 1.0d0
            inda=1
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. health(is,t) == ih .and. educ(is) == ie) then
                        sagrid(inda)=assets(is,t)
                        inda=inda+1
                    endif 
                enddo
            enddo
            call compressdist(sagrid,denom,nc,Swedistcomp2,Sweveccomp2)     
            call quart(Swedistcomp2,Sweveccomp2,nc,wequart10y_edu_health(j,ih,ie,:))           
  
        end do
    enddo
        
    end do      ! j
    
    ! Effort by Wealth Quartiles in Long Age Groups
    do j = 1,3
        if (j==1) then
            s = 1
        elseif (j==2) then
            s = 11
        else
            s = 21
        endif
        e = s+9
        
    ! Wealth Quartiles by Education
        do ie = 1,ne
            denom = 0
            do t = s,e
                 do is = 1,ns
                    if (survival(is,t) == 1 .and. educ(is) == ie) then
                        denom = denom + 1
                    endif
                enddo
            enddo
            
            if(allocated(sagrid))deallocate(sagrid)
            if(allocated(sadistr))deallocate(sadistr)
            allocate(sagrid(denom),sadistr(denom))
            sadistr = 1.0d0
            inda=1
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. educ(is) == ie) then
                        sagrid(inda)=assets(is,t)
                        inda=inda+1
                    endif 
                enddo
            enddo
            call compressdist(sagrid,denom,nc,Swedistcomp,Sweveccomp)     
            call quart(Swedistcomp,Sweveccomp,nc,wequart20y_edu(j,ie,:))           
            
            ! Assign groups and commpute average effort
            wequart = 0
            do t=s,e
            do is=1,ns
                if (survival(is,t) == 1 .and. educ(is) == ie) then
                    if (assets(is,t) <= wequart20y_edu(j,ie,1)) then
                        wequart(is,t) = 1
                    elseif ((assets(is,t) > wequart20y_edu(j,ie,1)) .and. (assets(is,t) <= wequart20y_edu(j,ie,2))) then
                        wequart(is,t) = 2
                    elseif ((assets(is,t) > wequart20y_edu(j,ie,2)) .and. (assets(is,t) <= wequart20y_edu(j,ie,3))) then
                        wequart(is,t) = 3
                    else
                        wequart(is,t) = 4
                    endif
                endif
            enddo
        enddo
                
        enddo ! ie
    enddo   ! j

    write(*, '(x)')
    write(*, '(a)') ' Model Simulation Averages (4 age groups groups) of HEALHTY and UNHEALTHY by EDUCATION'
    write(*, '(a)') ' Age Group               HS                           CL'
    write(*, '(a)') ' Age Group  Avr.Earn(H) Avr.Earn(U)  Avr.Earn(H) Avr.Earn(U)'
    do j = 1,4
        write(*, '(x,i2,a,4f10.4)') j, '     |  ', long_avr_earn10_health_edu(j,1,1)*norm/1000, long_avr_earn10_health_edu(j,2,1)*norm/1000, long_avr_earn10_health_edu(j,1,2)*norm/1000, long_avr_earn10_health_edu(j,2,2)*norm/1000           
    enddo
    
    write(*, '(x)')
    write(*, '(a)') ' Model Simulation Employment (4 age groups groups) of HEALHTY and UNHEALTHY by EDUCATION'
    write(*, '(a)') ' Age Group               HS                           CL'
    write(*, '(a)') ' Age Group  Avr.Empl(H) Avr.Empl(U)  Avr.Empl(H) Avr.Empl(U)'
    do j = 1,4
        write(*, '(x,i2,a,4f10.4)') j, '     |  ', long_avr_work10_health_edu(j,1,1), long_avr_work10_health_edu(j,2,1), long_avr_work10_health_edu(j,1,2), long_avr_work10_health_edu(j,2,2)           
    enddo
        
        
    write(*, '(x)')
    write(*, '(a)') ' Model Simulation Averages (4 age groups groups) of HEALHTY and UNHEALTHY'
    write(*, '(a)') ' Age Group  EmplShare(H) EmplShare(U) WorkPart(H) WorkPart(U) Avr.Earn(H) Avr.Earn(U)'
    do j = 1,4
        write(*, '(x,i2,a,6f10.4)') j, '     |  ', long_avr_work10(j,1), long_avr_work10(j,2), long_avr_work_part10(j,1), long_avr_work_part10(j,2), long_avr_earn10(j,1)*norm/1000, long_avr_earn10(j,2)*norm/1000        
    enddo
    
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages (7 age groups groups) of HEALHTY and UNHEALTHY'
    write(*, '(a)') ' Age Group Bad Health(HS) Bad Health (CL) '
    do j = 1,7
        write(*, '(x,i2,a,4f10.4)') j, '     |  ', long_avr_health10(j,1), long_avr_health10(j,2)
    enddo
    
    
    write(*, '(x)')
    write(*, '(a)') ' Comparison of Average Wealth by Wealth Terciles for 6 age groups by Health (in thousand euros)'
    write(*, '(a)') '                   Healthy         |       Unhealthy         |       Gap '       
    write(*, '(a)') ' Age          Bottom  Middle  Top  |   Bottom  Middle  Top   |   Bottom  Middle  Top '
    do j=1,7
        write(*, '(x,i2,a,3f8.2,a,3f8.2,a,3f8.2)') j, '     |  ', long_avr_wealth10_ter(j,1,1)*norm/1000, long_avr_wealth10_ter(j,1,2)*norm/1000, long_avr_wealth10_ter(j,1,3)*norm/1000, ' | ', long_avr_wealth10_ter(j,2,1)*norm/1000, long_avr_wealth10_ter(j,2,2)*norm/1000, long_avr_wealth10_ter(j,2,3)*norm/1000, ' | ', (long_avr_wealth10_ter(j,1,1)-long_avr_wealth10_ter(j,2,1))*norm/1000, (long_avr_wealth10_ter(j,1,2)-long_avr_wealth10_ter(j,2,2))*norm/1000, (long_avr_wealth10_ter(j,1,3)-long_avr_wealth10_ter(j,2,3))*norm/1000        
    end do
    
    
    write(*, '(x)')
    write(*, '(a)') ' Unconditional Average Wealth by Wealth Terciles at age 55-64 (in thousand euros)'
    write(*, '(a)') '         Bottom   Middle   Top  '
    write(*, '(a,3f10.4)') '   Model: ', avr_wealth10_ter4(1)*norm/1000, avr_wealth10_ter4(2)*norm/1000, avr_wealth10_ter4(3)*norm/1000
    
    
    write(*, '(x)')
    write(*, '(a)') ' Comparison of Wealth Quartiles each period by Health (in thousand euros)'
    write(*, '(a)') '                   Healthy         |       Unhealthy         |       Gap '        
    write(*, '(a)') ' Age          25th  50th  75th     |   25th  50th  75th      |     25th  50th  75th   '
    do t=1,7
        write(*, '(x,i2,a,3f8.2,a,3f8.2,a,3f8.2)') t, '     |  ', wequart10y(t,1,1)*norm/1000, wequart10y(t,1,2)*norm/1000, wequart10y(t,1,3)*norm/1000, ' | ', wequart10y(t,2,1)*norm/1000, wequart10y(t,2,2)*norm/1000, wequart10y(t,2,3)*norm/1000, ' | ', (wequart10y(t,1,1)-wequart10y(t,2,1))*norm/1000, (wequart10y(t,1,2)-wequart10y(t,2,2))*norm/1000, (wequart10y(t,1,3)-wequart10y(t,2,3))*norm/1000
    end do
    
    write(*, '(a)') ' Unconditional Wealth Quantiles  (in thousand euros)'
    write(*, '(x)')
    write(*, '(a)') ' Age          25th  50th  75th      '
    do t=1,7
    write(*, '(x,i2,a,3f8.2)') t, '     |  ', wequart10y_av(t,1)*norm/1000, wequart10y_av(t,2)*norm/1000, wequart10y_av(t,3)*norm/1000
    enddo
    
    write(*, '(x)')
    write(*, '(a)') ' Wealth Quartiles each period by Education (in thousand euros)'
    write(*, '(a)') '                   HS         |       CL         |       Gap '        
    write(*, '(a)') ' Age          25th  50th  75th     |   25th  50th  75th      |     25th  50th  75th   '
    do t=1,7
        write(*, '(x,i2,a,3f8.2,a,3f8.2,a,3f8.2)') t, '     |  ', wequart10y_edu(t,1,1)*norm/1000, wequart10y_edu(t,1,2)*norm/1000, wequart10y_edu(t,1,3)*norm/1000, ' | ', wequart10y_edu(t,2,1)*norm/1000, wequart10y_edu(t,2,2)*norm/1000, wequart10y_edu(t,2,3)*norm/1000, ' | ', (wequart10y_edu(t,2,1)-wequart10y_edu(t,1,1))*norm/1000, (wequart10y_edu(t,2,2)-wequart10y_edu(t,1,2))*norm/1000, (wequart10y_edu(t,2,3)-wequart10y_edu(t,1,3))*norm/1000
    end do
    
    
    write(*, '(x)')
    write(*, '(a)') ' Wealth Quartiles each period by Health and Education (in thousand euros)'
    write(*, '(a)') '                   Healthy - HS       |       Unhealthy - HS        |       Healthy - CL      |      Unhealthy - CL '        
    write(*, '(a)') ' Age          25th  50th  75th     |   25th  50th  75th      |     25th  50th  75th      |     25th  50th  75th  '
    do t=1,6
        write(*, '(x,i2,a,3f8.2,a,3f8.2,a,3f8.2,a,3f8.2)') t, '     |  ', wequart10y_edu_health(t,1,1,1)*norm/1000, wequart10y_edu_health(t,1,1,2)*norm/1000, wequart10y_edu_health(t,1,1,3)*norm/1000, ' | ', &
            wequart10y_edu_health(t,2,1,1)*norm/1000, wequart10y_edu_health(t,2,1,2)*norm/1000, wequart10y_edu_health(t,2,1,3)*norm/1000, ' | ', &
            wequart10y_edu_health(t,1,2,1)*norm/1000, wequart10y_edu_health(t,1,2,2)*norm/1000, wequart10y_edu_health(t,1,2,3)*norm/1000, ' | ', &
            wequart10y_edu_health(t,2,2,1)*norm/1000, wequart10y_edu_health(t,2,2,2)*norm/1000, wequart10y_edu_health(t,2,2,3)*norm/1000
    end do
    
    write(*, '(x)')
    write(*, '(a)') ' Statistics by Wealth Quartiles - Model'
    write(*, '(a)') '                   Share in Bad Health           |'
    write(*, '(a)') ' Age           1st     2nd     3rd     4th       |'
    do j=1,7
        write(*, '(x,i2,a,4f10.3,a,4f10.3)') j, '     |  ',  hequ(j,1), hequ(j,2), hequ(j,3), hequ(j,4), ' | '
    end do
       
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages (7 age groups groups) of HS and CL'
    write(*, '(a)') ' Age Group  EmplShare(HS) EmplShare(CL) Avr.Earn(HS) Avr.Earn(CL)'
    do j = 1,5   
        write(*, '(x,i2,a,6f10.4)') j, '     |  ', long_avr_work10_edu(j,1), long_avr_work10_edu(j,2), long_avr_earn10_edu(j,1)*norm/1000, long_avr_earn10_edu(j,2)*norm/1000
    end do
    
    ! Calculate average based on longer age groups
    do ih = 1,nh
        long_avr_earn(1,ih) = sum(earn(:,1:tj-1), mask=survival(:,1:tj-1)==1 .AND. health(:,1:tj-1)==ih) / sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih)     ! Average earnings across everyone from that health group
        long_avr_work(1,ih) = sum(work(:,1:tj-1), mask=health(:,1:tj-1)==ih .AND. survival(:,1:tj-1)==1) / dble(sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih))     ! This is now average hours worked can be 0, 0.5, or 1 (for 2 years of full time work)
        long_avr_work_part(1,ih) = sum(survival(:,1:tj-1), mask=(health(:,1:tj-1)==ih) .AND. (work(:,1:tj-1)==0.5)) / dble(sum(survival(:,1:tj-1), mask=(health(:,1:tj-1)==ih) .and. (work(:,1:tj-1)>0.0d0)))
        long_avr_cons(1,ih) = sum(consum(:,1:tj-1), mask=survival(:,1:tj-1)==1 .AND. health(:,1:tj-1)==ih) / sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih)  
!        
        long_avr_earn(2,ih) = sum(earn(:,tj:tr-1), mask=survival(:,tj:tr-1)==1 .AND. health(:,tj:tr-1)==ih) / sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih)     ! Average earnings across everyone from that health group
        long_avr_work(2,ih) = sum(work(:,tj:tr-1), mask=health(:,tj:tr-1)==ih .AND. survival(:,tj:tr-1)==1) / dble(sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih))
        long_avr_work_part(2,ih) = sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih .AND. work(:,tj:tr-1)==0.5) / dble(sum(survival(:,tj:tr-1), mask=(health(:,tj:tr-1)==ih) .and. (work(:,tj:tr-1)>0.0d0)))
        long_avr_cons(2,ih) = sum(consum(:,tj:tr-1), mask=survival(:,tj:tr-1)==1 .AND. health(:,tj:tr-1)==ih) / sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih)  
!    
        ! Average Consumption here only until age 84
        long_avr_cons(3,ih) = sum(consum(:,tr:nt), mask=survival(:,tr:nt)==1 .AND. health(:,tr:nt)==ih) / sum(survival(:,tr:nt), mask=health(:,tr:nt)==ih)      
        
    enddo
    
    do ih = 1,2
        avr_cons_work(ih) =  sum(consum(:,1:tr-1), mask=survival(:,1:tr-1)==1 .AND. health(:,1:tr-1)==ih) / sum(survival(:,1:tr-1), mask=health(:,1:tr-1)==ih)      
    enddo
    
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages (3 age groups groups) of HEALHTY and UNHEALTHY'
    write(*, '(a)') ' Age Group EmplShare(H) EmplShare(U) WorkPart(H) WorkPart(U) Avr.Earn(H) Avr.Earn(U)'
    write(*, '(x,i2,a,8f10.4)') 1, '     |  ', long_avr_work(1,1), long_avr_work(1,2), long_avr_work_part(1,1), long_avr_work_part(1,2), long_avr_earn(1,1)*norm/1000, long_avr_earn(1,2)*norm/1000
    write(*, '(x,i2,a,8f10.4)') 2, '     |  ', long_avr_work(2,1), long_avr_work(2,2), long_avr_work_part(2,1), long_avr_work_part(2,2), long_avr_earn(2,1)*norm/1000, long_avr_earn(2,2)*norm/1000
    
   write(*, '(x)')
   write(*, '(a)') ' Age Group Avr.Cons(H) Avr.Cons(U) '
   write(*, '(x,i2,a,8f10.4)') 1, '     |  ', long_avr_cons(1,1)*norm/1000, long_avr_cons(1,2)*norm/1000
   write(*, '(x,i2,a,8f10.4)') 2, '     |  ', long_avr_cons(2,1)*norm/1000, long_avr_cons(2,2)*norm/1000
   
   write(*, '(x)')
   write(*, '(a)') ' Age Group Earn Ratio(H/U) Cons.Ratio(H/U) '
   write(*, '(x,i2,a,8f10.4)') 1, '     |  ', (long_avr_earn(1,1)*norm/1000)/(long_avr_earn(1,2)*norm/1000), (long_avr_cons(1,1)*norm/1000)/(long_avr_cons(1,2)*norm/1000)
   write(*, '(x,i2,a,8f10.4)') 2, '     |  ', (long_avr_earn(2,1)*norm/1000)/(long_avr_earn(2,2)*norm/1000), (long_avr_cons(2,1)*norm/1000)/(long_avr_cons(2,2)*norm/1000)

   write(*, '(x)')
   write(*, '(a)') ' Consumption Ratio for all workers (H/U)'
   write(*, '(x,8f10.4)') (avr_cons_work(1)*norm/1000) / (avr_cons_work(2)*norm/1000) 
   
! -------------     
   
    
    ! ----------------
    ! Get simulated wealth and earnings Gini coefficients    
    
    denomall = sum(survival(:,1:nt))              ! number of all alive people 
    denom2564 = sum(survival(:,1:tr-1))           ! number of alive people aged 25-64
    denomwo = sum(survival(:,1:tr-1), mask=work(:,1:tr-1)>0.0d0)         
    
    
    allocate(Swevec(denomall),Swedist(denomall),Sincvec(denomwo),Sincdist(denomwo))
    Swedist = 1.0d0
    Sincdist = 1.0d0

    inda = 0
    indinc = 0
    do t=1,nt
        do is=1,ns            
            if (survival(is,t) == 1) then
                inda=inda+1
                Swevec(inda)=assets(is,t)                
                if (t<tr) then
                    if (work(is,t) > 0.0d0) then
                        indinc=indinc+1
                        Sincvec(indinc)=log(earn(is,t)*norm)        ! Log Earings of Workers
                        if (earn(is,t) == 0.0d0) then
                            print*, ' Error: Zero earn recorded'
                            print*, earn(is,t)
                            print*, work(is,t)
                        end if                        
                    end if
                end if
            end if 
        end do
    end do
    
    call compressdist(Swevec,denomall,nc,Swedistcomp,Sweveccomp)    
    
    call ginicoef(Swedistcomp,Sweveccomp,nc,5,Sginiwe,shareW,negW)

    write(*,"('  Simulated Wealth Gini (all age)  = ', F6.3)") Sginiwe

    Sginiy = varsim(Sincvec)
    write(*,"('  Variance of log Earnings (of workers) = ', F6.3)") Sginiy
    
    ! Calculate variance of log earnings by education
    indinc = 0
    do ie = 1,ne
        do t=1,(tr-1)
            do is=1,ns            
                if (survival(is,t) == 1 .and. educ(is) == ie) then
                    if (work(is,t) > 0.0d0) then
                        indinc=indinc+1
                        Sincvec(indinc)=log(earn(is,t)*norm)        ! Log Earnings of Workers
                        if (earn(is,t) == 0.0d0) then
                            print*, ' Error: Zero Earnings recorded'
                        end if                        
                    end if
                end if
            end do
        end do
        if (ie==1) then            
            Sginiy_HS = varsim(Sincvec)
        else
            Sginiy_CL = varsim(Sincvec)
        end if
        
        write(*,'(a,i2,f10.4)') '  Variance of log Earnings (of workers) with education type:', ie, varsim(Sincvec)
    enddo
    
  
    
    ! ---------------------
    ! Calculate Total Wealth to Earnings (age<65) Ratio 

    totalwealth = sum(assets, mask=survival(:,1:tr-1) == 1) 
    totaly = sum(earn, mask=survival(:,1:tr-1) == 1)
    w2y = totalwealth/totaly 
     
    write(*,"('  Avg.Wealth = ', f6.3,' or', f8.0,' EUR   Avg.Earnings (all ages,inc. 0) =', f6.3,' or', f8.0,' EUR')") &
        totalwealth/dble(denom2564), (totalwealth/dble(denom2564))*norm, avrgearn_sim, avrgearn_sim*norm    
    write(*,"('  Wealth(avg)/Earnings =', F7.4)") w2y

    ! Pension
    bbb = sum(pensionvec)/sum(survival(:,tr))
    
    penavg = bbb/avrgearn_sim
    
    write(*,'(x)')
    write(*,"('  Average Pensions = ', F8.3, '  Agg. Pension Replacement rate (%) = ', F7.2)") bbb, penavg*100.d0
    write(*,'(a,f10.3)') ' Average wealth in retirement years ' , (sum(assets(:,tr:nt), mask=survival(:,tr:nt)==1) / sum(survival(:,tr:nt)) ) *norm/1000
    
    ! ---------------------
    ! Calculate model-implied life expectancy at age 25
    le25 = 2.0d0*sum(agedeath-1) / dble(ns)    ! 2 year period
    
    write(*,'(a,f10.4)') ' Average Life Expectancy at age 25: ', le25

    ! Calculate Value of a statistical life year using the formula in Glover et al.
    
    ! 1a. Get model-implied 2-year average consumption, hours, efforts for a representative agent
    avgcons = sum(consum, mask=survival==1) / denomall
    avghours = sum(work, mask=survival(:,1:tr-1)==1) / denomall   ! 2026 Update: mask restricted to working ages to conform with work(ns,tr-1) 

    ! 1b. Get average kappa, phi, xi for a representative agent
    do is=1,ns
        do t=1,nt
            if (survival(is,t) == 1 .and. consum(is,t) > 0) then
                if (health(is,t)==1) then
                    mucontemp = 1.0d0
                else
                    mucontemp = conp
                end if                
                marg_consutil(is,t) = mucontemp
                                
                if (t<tr) then                    
                    phisim(is,t) = phimat(t,educ(is),health(is,t))  
                end if                
            end if
        end do
    end do
    
    avgkappa = sum(marg_consutil, mask=survival == 1 .and. consum >0) / sum(survival, mask=survival == 1 .and. consum >0)
    avgphi = sum(phisim, mask=survival(:,1:tr-1)==1) / sum(survival(:,1:tr-1))   ! 2026: mask restricted to working ages to conform with phisim(ns,tr-1)
    
    ! 2a. Obtain LHS of vsly equation 
    if (sigma == 1d0) then
        vsl2y = avgkappa*(log(avgcons)) + bbval - avgphi*lsutil(avghours,gamma)
    else
        vsl2y = avgkappa*((avgcons)**(1.0d0-sigma))/(1.0d0-sigma) + avgkappa*bbval - avgphi*lsutil(avghours,gamma) 
    endif    
    ! 2b. then divide it by MU of RHS    
    vsl2y = vsl2y / ( avgkappa*(avgcons**(-sigma)) )
    
    write(*, '(x)') 
    write(*,'(a,f18.4)') ' Average Value of a 2 Statistical Life Years: ', vsl2y
    write(*,'(a,f15.4,a,f10.4)') ' This corresponds to ', vsl2y/(avgcons), ' times average 2-year consumption. Target: ', 8.49d0

    
    ! ************************************************************************************************
    ! Target moments
    
    ! 1-4. 2-year employment share among HEALTHY (4 moments) 
    do g=1,4
        y(g) = long_avr_work10(g,1)      
    end do
    
    ! 5-8. 2-year employment share among UNHEALTHY (4 moments)        
    do g=5,8
        j = g-4
        y(g) = long_avr_work10(j,2)
    end do
     
    ! 9-14. Median wealth age group 1-6, relative to norm=winit(2)
    y(9) = wequart10y_av(1,2)
    y(10) = wequart10y_av(2,2)    
    y(11) = wequart10y_av(3,2)   
    y(12) = wequart10y_av(4,2)
    y(13) = wequart10y_av(5,2)    
    y(14) = wequart10y_av(6,2)   
        
    ! 15. Education gradients in employment       
    y(15) = avrempCL/avrempHS
    
    ! 16. VSL multiple
    y(16) = vsl2y/(avgcons)   

    ! 17. Wealth Gini
    y(17) = Sginiwe
    
    ! 18. c(healthy)/c(unhealthy) for number of household = 1 or 2    
    y(18) = avr_cons_work(1)/avr_cons_work(2)       
    
    ! 19-34. vector of earnings profiles in 1,000 EUR (4x4)
    do g=19,22
        ! HS-Healthy
        j = g-18
        y(g) = long_avr_earn10_health_edu(j,1,1)*norm/1000.0d0
    end do
    
    do g=23,26
        ! HS-Unhealthy
        j = g-22
        y(g) = long_avr_earn10_health_edu(j,2,1)*norm/1000.0d0
    end do
    
    do g=27,30
        ! CL-Healthy
        j = g-26
        y(g) = long_avr_earn10_health_edu(j,1,2)*norm/1000.0d0
    end do
    
    do g=31,34
        ! CL-Unhealthy
        j = g-30
        y(g) = long_avr_earn10_health_edu(j,2,2)*norm/1000.0d0
    end do
        
    ! 35. Variance of log earnings
    y(35) = Sginiy
        
    ! 36. Pension replacement rate
    y(36) = penavg
       
    
    objval = 0.0d0
    do g = 1,nm
        if (efweight==1) then
            ! Using efficient weights
            if (parm_est==0) then
                objval = objval +(tdatase(g)**(-2.0d0))*((y(g)-tdatamm(g))**2.0d0)
            else
                objval = objval + weight_mat(g,parm_est)*(tdatase(g)**(-2.0d0))*((y(g)-tdatamm(g))**2.0d0)
            end if            
        else
            objval = objval + log(y(g)/tdatamm(g))**2.0d0
        end if
        
    end do    
    
    
    
    call system_clock(ctf, cr)
    eptime = (dble(ctf-ct0)/cr)/60.0d0
    
    open(10, file="./equil_summary.txt")
    
    !if (estimation==1) then
        
        write(10,*) objval
        ! parameters    
	    do g=1,nparam
            write(10,*) x_param(g)
        end do
        
        do g=1,nmom
            write(10,*) y(g)
        end do
    !end if
    
    close(10)

    
    write(*, '(x)')
    write(*,"('===  Elapsed time for one evaluation (min) = ', F6.2,'    Obj value = ', F14.5,'    Eval num = ', I5,' ===')") eptime, objval, evalnum
    
    evalnum = evalnum + 1
     
    
    
    if (recordsimstat==1) then
     
    ! ********************************************
    ! Replication: model outputs
    
    ! ********************************************
    ! FIGURE 9: RIGHT PANEL ONLY
    ! ********************************************
    open(19, file="Figure9model2.txt")
    
    ! median wealth among HEALTHY
    do g=2,6
        write(19,*) wequart10y(g,1,2)*norm/1000
    end do
    ! median wealth among UNHEALTHY
    do g=2,6
        write(19,*) wequart10y(g,2,2)*norm/1000
    end do
    
    close(19)
    
    end if

    end subroutine sub_main



