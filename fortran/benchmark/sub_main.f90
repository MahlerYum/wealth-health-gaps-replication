subroutine sub_main(x_param,y,evalnum,exper_run)
    
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
    integer, intent(in) :: exper_run
    
    ! Define variables
    real(8) V_adjtemp, V_nadjtemp, V_w_temp, V_n_temp, V_n_nadj_temp, V_w_nadj_temp, chitemp, trash, objval, trash1, trash2, trash3, effa, effb, effc, effd, fc, fd, effz, diff, valtemp, mucontemp
    real(8) agrid(na), xgrid(nx), zgrid(nz), pij(nx,nx), pz(nz,nz), fgrid(nf), fnowgrid(nfnow), bgrid(nb), zeta(nt,ne,nh), chimax, chigrid(nchi,nt) 
    real(8) yt(tr-1,ne,nh), wageg(8,ne), ggrid(8), cfg(4,7,ne), yg(8,ne), ximat(nt,ne,nh), phimat(tr-1,ne,nh), fmeanvec(nt)    
    real(8) pi0(6,nh,ne,neta), l1(6,nh,ne,neta), l2(6,nh,ne,neta), fpast, amin, amax, avrgearn, hours, Vvec_adj(nw), Vvec_nadj(nw)  
    real(8) beta, anow, ap, ahigh, pbnow, spnow, pinow, xinow, phinow, awhigh, anhigh, ynow, ynow_net, fnow, eval(nh), aptemp, Vvaltemp, Mvaltemp, aptempf(nfnow), vtempf(nfnow)
    real(8) prxi(nxi), prtemp2(nxi,nxi), xigrid(nxi), txigrid(4), tnugrid(4),  xiHShval(4), xiHSuval(4), xiCLhval(4), xiCLuval(4), nuhval(4), nuuval(4)    
    integer t, i, j, indf, ix, jx, ia, ja, ik, jk, iz, jz, ib, jb, ih, jh, ie, ichi, bf, ind, m, g, ind2, ind3, ixi, itheta, ieta, itype, indmat(3,ntype), maxind, iter     
    real(8), allocatable :: Vmat(:,:,:,:,:,:,:), Vmatnext(:,:,:,:,:,:), V_adj(:,:,:,:,:,:,:,:,:), V_nadj(:,:,:,:,:,:,:,:,:), ap_adj(:,:,:,:,:,:,:,:,:), ap_nadj(:,:,:,:,:,:,:,:,:), f_adj(:,:,:,:,:,:,:,:,:)   
            
    integer :: st, ct1, ct2, ct0, ctf, cr, st1
    real(8) eptime    
    
    ! Simulation
    integer :: is, count_active
    integer, allocatable :: apos(:,:), zpos(:,:), educ(:), thetapos(:), etapos(:), dpos(:), health(:,:), survival(:,:), fpos(:,:), chipos(:,:), adj(:,:), spike10(:,:), spike20(:,:), xipos(:,:)
    real(8), allocatable :: assets(:,:), effort(:,:), earn(:,:), savings(:), consum(:,:), bequests(:,:), effhabit(:,:), work(:,:)    
    real(8) cumchi(nchi), WtempS, NtempS, hy_contr, le_contr, fh_contr
    real(8) a, b, b1, b2, cumzss(nz), onedraw(1), cumpz(nz,nz)
    
    ! Calculation of Moments
    real(8) num, numa, numh, numwe, numf, numwo, numwo_part, numinc, norm, effcorr(3), weight_mat(nmom,nparam)
    real(8) avwealth2y(nt,nh), avcons2y(nt,nh), aveff2y(nt,nh), avwork2y(tr-1,nh), avearn2y(tr-1,nh),avrhealth2y(nt),avrhealthedu2y(nt,ne), aveff2yall(nt)
    real(8) wequart2y(nt,nh,3), wequart10y_av(7,3), hequ(7,4), efqu(7,4), earn_contr, avrg_lsdisutil, avrg_fdisutil,vsl2y, asqu(3,4,ne)
    real(8) avreffHS, avreffCL, avrempHS, avrempCL, num1, num2, num3, num4, shareW(5), negW, w2y_h, w2y_u, sum1, stdeff, avrdeltaeff
    integer :: denom1, denom2, denom3, denom4, denom, denomall, denom2564, denomwo, inda, indinc, s, e, ihw, t1, t2
    
    real(8) numcons, avrcons25, le25, svl, totalwealth, totaly, w2y, penavg, aa, bbb, Sginiwe, Sginiy, Sginiy_HS, Sginiy_CL, we_contr
    real(8), allocatable :: sagrid(:),sygrid(:),sadistr(:),sydistr(:),lsdisutil(:,:),fdisutil(:,:), Swevec(:), Sincvec(:), Swedist(:), Sincdist(:), wequart(:,:), effcomp(:), effcomp2(:)
    
    real(8) Swedistcomp(nc), Sweveccomp(nc), Sincdistcomp(nc), Sincveccomp(nc), Swedistcomp2(nc), Sweveccomp2(nc)
    real(8) long_avr_earn(2,nh), long_avr_work(2,nh), long_avr_work_part(2,nh), long_avr_eff(3,nh), long_avr_cons(3,nh)
    
    real(8) avrgearn_sim, avrcons, avrriskdeath, vsl, avreffAll, hwgrid(nw), avrg_consutil, avrg_marg_consutil, avgcons, avghours, avgf, avgkappa, avgphi, avgxi    
    real(8) wequart10y(7,nh,3), long_avr_earn10(4,nh), long_avr_work10(4,nh), long_avr_work_part10(4,nh), long_avr_eff10(7,nh), long_avr_cons10(7,nh), long_avr_health10(7,ne), &
        long_avr_eff10_edu(7,ne), long_avr_earn10_edu(4,ne), long_avr_work10_edu(4,ne), long_eff_noadj(3), long_eff_pos(3), long_eff_neg(3), pos_spike10, neg_spike10, pos_spike20, neg_spike20, &
        long_pos_spike10(3), long_neg_spike10(3), long_pos_spike20(3), long_neg_spike20(3), eff_pos, eff_neg
    real(8) long_avr_eff10_health_edu(7,2,2), long_avr_earn10_health_edu(4,2,2), long_avr_work10_health_edu(4,2,2)
    real(8) long_avr_wealth10_ter(7,2,3), avr_cons_work(nh), avr_wealth10_ter4(3), wequart10y_edu_health(7,nh,ne,3), wequart10y_edu(7,ne,3)
    real(8) efqu_educ(3,4,ne), wequart20y_edu(3,ne,3), ef_contr
    
    type(VSL_STREAM_STATE) :: stream

    ! Calibration
    !
    real(8) phi1_i, phi1_s, phi1_sq, phi2_i, phi2_s, phi2_sq, xi1_i, xi1_s, xi1_sq, xi2_i, xi2_s, xi2_sq, mincon, haddft, y1_HS, y1_CL
    real(8) wagep_HS, wagep_CL, conp, betamean, betadev, psi, penre, beta1, beta2, nu1_i, nu1_s, nu1_sq, nu2_i,nu2_s,nu2_sq, xiadd, xiadd_1, xiadd_2, xiadd_3, xiadd_4, xiadd_5, xiadd_6, nuadd, tt, sigma
    real(8) xiHSh_1, xiHSh_2, xiHSh_3, xiHSh_4, xiHSu_1, xiHSu_2, xiHSu_3, xiHSu_4, xiCLh_1, xiCLh_2, xiCLh_3, xiCLh_4, xiCLu_1, xiCLu_2, xiCLu_3, xiCLu_4, chi_1, chi_2, bbval, chi_3, sdxi, sdztemp
    real(8) ytHS_s, ytHS_sq, ytCL_s, ytCL_sq, sigx, nuh_1, nuh_2, nuh_3, nuh_4, nuu_1, nuu_2, nuu_3, nuu_4
    
    ! 2026 Update: work array for continuation value at new effort
    real(8) :: EV_adj(na,nh)   
    logical :: equal_effort_now, equal_effort_prev
    
    character(len=100) :: filename
    
    ! spline
    real(8), allocatable :: EV(:,:), cff(:,:), Tximat(:,:), Tnumat(:,:), cfxiHSh(:,:), cfxiCLh(:,:), cfnuh(:,:), cfxiHSu(:,:), cfxiCLu(:,:), cfnuu(:,:)
    real(8), allocatable :: agedeath(:), pensionvec(:), marg_consutil(:,:), xisim(:,:), phisim(:,:), effcomptemp(:,:)
    
    integer :: specs(nspecs),recordsimdata,write_fmeanvec,exp_prod,exp_lab,equalsurv,exp_no_consp,exp_wealthtax_pension_ret,use_BM_health 
    integer :: use_BM_assets,use_BM_work,use_BM_survival,use_BM_effort,alleffort_young,alleffort_mid,alleffort_retire,use_CF_health 
    
	allocate(EV(na,nh), Tximat(2,2), Tnumat(2,2), cfxiHSh(4,2), cfxiCLh(4,2), cfnuh(4,2), cfxiHSu(4,2), cfxiCLu(4,2), cfnuu(4,2))           
    allocate(agedeath(ns), pensionvec(ns), marg_consutil(ns,nt), xisim(ns,nt), phisim(ns,tr-1), effcomptemp(ns*nt,4))
    

    
    ! Call the experiment specifications
    ! 2026 Update: exp_specs called for every run (incl. baseline) so all experiment flags are initialised
    call exp_specs(specs,exper_run) 
    if (exper_run>0) then
        write(*,'(a,i2)') 'Experiments Run = ', exper_run
    end if
    
    recordsimdata = specs(1)
    write_fmeanvec = specs(2)
    exp_prod= specs(3)
    exp_lab= specs(4)
    equalsurv= specs(5)
    exp_no_consp = specs(6)
    exp_wealthtax_pension_ret = specs(7)
    use_BM_health = specs(8)
    use_BM_assets = specs(9)
    use_BM_work = specs(10)
    use_BM_survival= specs(11)
    use_BM_effort= specs(12)
    alleffort_young = specs(13)
    alleffort_mid= specs(14)
    alleffort_retire = specs(15)
    use_CF_health = specs(16)
   
    
    
    !do i=1,nspecs
    !    write(*,'(a,i2)') 'Experiments Specification = ', specs(i)
    !enddo
    
    
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
        
    xiHSh_1 = x_param(10)
    xiHSh_2 = x_param(11)
    xiHSh_3 = x_param(12)
    xiHSh_4 = x_param(13)
    
    xiHSu_1 = x_param(14)
    xiHSu_2 = x_param(15)    
    xiHSu_3 = x_param(16)
    xiHSu_4 = x_param(17)
    
    xiCLh_1 = x_param(18)
    xiCLh_2 = x_param(19)
    xiCLh_3 = x_param(20)
    xiCLh_4 = x_param(21)
        
    xiCLu_1 = x_param(22)
    xiCLu_2 = x_param(23)
    xiCLu_3 = x_param(24)
    xiCLu_4 = x_param(25)  
        
    chi_1 = x_param(26)
    chi_2 = x_param(27)    
    
    sigx = x_param(28)              
    penre = x_param(29)
    
    betamean = x_param(30)
    betadev = x_param(31)
        
    bbval = x_param(32)
    conp = x_param(33)    
    psi = x_param(34)
    
    ytHS_s = x_param(35)
    ytHS_sq = x_param(36)
    ytCL_s = x_param(37)
    ytCL_sq = x_param(38)
    
    wagep_HS = x_param(39)
    wagep_CL = x_param(40)
    
    y1_HS = x_param(41)
    y1_CL = x_param(42)
    
    sigma = sigma0
           
    haddft = 0.0d0
        
    sdxi= 0.0d0
    chi_3 = 0.0d0    
    
    
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
    
    ! Equalizing survival prob
    if (equalsurv==1) then        
        zeta(:,:,2) = zeta(:,:,1)        
    end if
    if (equalsurv==2) then
        zeta(:,:,1) = zeta(:,:,2)        
    end if
    
    ! Load Benchmark Effort Averages
    open(17, file="fmeanvec.txt")
    do i = 1,nt
        read(17,*) fmeanvec(i)
    enddo
    close(17)
    
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
    
    xiHShval = (/xiHSh_1, xiHSh_2, xiHSh_3, xiHSh_4 /)
    xiHSuval = (/xiHSu_1, xiHSu_2, xiHSu_3, xiHSu_4 /)
    xiCLhval = (/xiCLh_1, xiCLh_2, xiCLh_3, xiCLh_4 /)
    xiCLuval = (/xiCLu_1, xiCLu_2, xiCLu_3, xiCLu_4 /)
    nuhval = (/nuh_1, nuh_2, nuh_3, nuh_4/)
    nuuval = (/nuu_1, nuu_2, nuu_3, nuu_4/)
    
    if (exp_lab == 1) then
        nuuval = nuhval        
    end if
        
    txigrid = (/1.0d0,12.0d0,20.0d0,31.0d0/)       
    tnugrid = (/1.0d0,8.0d0,13.0d0,20.0d0/) 
        
    ! setup splines
    Tximat = spbas(4-2,txigrid)    
    Tximat = invmat(Tximat)
    Tnumat = spbas(4-2,tnugrid)
    Tnumat = invmat(Tnumat)
    
    ! fit spline
    cfxiHSh = spfit(Tximat,xiHShval,4-2,txigrid)
    cfxiHSu = spfit(Tximat,xiHSuval,4-2,txigrid)
    cfxiCLh = spfit(Tximat,xiCLhval,4-2,txigrid)
    cfxiCLu = spfit(Tximat,xiCLuval,4-2,txigrid)
    cfnuh = spfit(Tnumat,nuhval,4-2,tnugrid)
    cfnuu = spfit(Tnumat,nuuval,4-2,tnugrid)
    
    do t=1,nt
        if (t<=30) then
        do ie=1,ne
            do ih=1,nh                
                if (ie==1) then
                    ! HS
                    if (ih==1) then
                        ! Healthy
                        call speva(cfxiHSh,dble(t),4-2,txigrid,valtemp,trash,trash)                        
                        ximat(t,ie,ih) = max(0.0d0,valtemp)
                    else
                        ! Unhealthy
                        call speva(cfxiHSu,dble(t),4-2,txigrid,valtemp,trash,trash)                                     
                        ximat(t,ie,ih) = max(0.0d0,valtemp)
                    end if
                else
                    ! CL
                    if (ih==1) then
                        ! Healthy
                        call speva(cfxiCLh,dble(t),4-2,txigrid,valtemp,trash,trash)
                        ximat(t,ie,ih) = max(0.0d0,valtemp)
                    else
                        ! Unhealthy
                        call speva(cfxiCLu,dble(t),4-2,txigrid,valtemp,trash,trash)
                        ximat(t,ie,ih) = max(0.0d0,valtemp)
                    end if
                end if
            end do
        end do
        else
            ! use the last value
            do ie=1,ne
            do ih=1,nh                
                if (ie==1) then
                    ! HS
                    if (ih==1) then
                        ! Healthy                        
                        ximat(t,ie,ih) = xiHShval(4)
                    else
                        ! Unhealthy                        
                        ximat(t,ie,ih) = xiHSuval(4)
                    end if
                else
                    ! CL
                    if (ih==1) then
                        ! Healthy                        
                        ximat(t,ie,ih) = xiCLhval(4)
                    else
                        ! Unhealthy                        
                        ximat(t,ie,ih) = xiCLuval(4)
                    end if
                end if
            end do
            end do
        end if        
            
    end do
     
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
     
    
    ! Define effort grid
    fgrid = linspace(fmin,fmax,nf)     
    ! separately for choice
    fnowgrid = linspace(fmin,fmax,nfnow)      ! 
   
    if (alleffort_young==1) then
        print*, '*** Counterfactual: Effort choices are the same for the young ***'        
    end if
    if (alleffort_mid==1) then
        print*, '*** Counterfactual: Effort choices are the same for the mid ***'        
    end if
    if (alleffort_retire==1) then
        print*, '*** Counterfactual: Effort choices are the same for the retired ***'        
    end if
    
    
    ! Define preference type grid
    bgrid(1) = betamean-betadev
    bgrid(2) = betamean+betadev
    
    if (nxi>2) then
        call tauchen(nxi,0.0d0,0.0d0,sdxi,2.0d0,xigrid,prtemp2)
        prxi = prtemp2(1,:)
        xigrid = exp(xigrid)    
    else
        xigrid(1) = 1.0d0
        prxi = 1.0d0
    end if
    
    write(*,"('  bgrid =', F10.5, F10.5)") bgrid(1), bgrid(2)
    
    ! Define stochastic effort adjustment costs grid (uniformly distributed)      
    do t = 1,nt                    
        chimax = max(0.0d0, chi_1*exp(chi_2*(dble(t)-1.0d0) + chi_3*(dble(t)-1.0d0)**2.0d0)  )
         chigrid(:,t) = linspace(chimin,chimax,nchi)                
    end do
    
    ! adjusted to fit our order of states (nb,ntheta,neta)
    indmat = dimreduce3_new(nb,ntheta,neta)
   
    ! Allocate Space        
    allocate(Vmat(na,nz,ne,ntype,nh,nf,nt), Vmatnext(na,nz,ne,ntype,nh,nf))  
    allocate(V_adj(na,nz,ne,ntype,nh,nf,nxi,nw,nt), V_nadj(na,nz,ne,ntype,nh,nf,nxi,nw,nt))
    allocate(ap_adj(na,nz,ne,ntype,nh,nf,nxi,nw,nt), ap_nadj(na,nz,ne,ntype,nh,nf,nxi,nw,nt), f_adj(na,nz,ne,ntype,nh,nf,nxi,nw,nt))
    
    call system_clock(ct1, cr)
    
    ! ----------------------------------------------------------------------------------------- !
    ! Backward induction    
    do t = nt,1,-1 
    
        if (t<nt) then
            Vmatnext = Vmat(:,:,:,:,:,:,t+1)                    
        else
            Vmatnext = 0.0d0
        end if
        
        if (mod(t,diagnum)==0 ) then
            write(*,'(a,x,i2)')'Solving value functions and decision rules at age',t
        end if       
        
        ! 2026 Update: flag whether efforts are equalised in t and t-1 (Fig. 11, A.5, Table A.13)
        equal_effort_now = .false.

        if ((alleffort_young  == 1 .and. t <= 10) .or. &
            (alleffort_mid    == 1 .and. t >= 11 .and. t < tr) .or. &
            (alleffort_retire == 1 .and. t >= tr)) then
            equal_effort_now = .true.
        end if


        equal_effort_prev = .false.

        if (t > 1) then
            if ((alleffort_young  == 1 .and. t-1 <= 10) .or. &
                (alleffort_mid    == 1 .and. t-1 >= 11 .and. t-1 < tr) .or. &
                (alleffort_retire == 1 .and. t-1 >= tr)) then
                equal_effort_prev = .true.
            end if
        end if
        
        
        if (t >= tr) then       
            ! Retirement age group only face uncertainty regarding health
            if (t >= 26) then
                j = 6               
            else
                j = 5
            endif
            
            !$omp parallel do default(shared) private(iz,ie,ib,itheta,ieta,itype,bf,indf,aptempf,vtempf,maxind,jh,ih,EV,EV_adj,ia,beta,anow,pbnow,xinow,phinow,spnow,ahigh,fnow,fpast,pinow,aptemp,mucontemp,hours,ixi,effa,effb,effc,effd,fc,fd,effz,diff,iter)                
            do bf = 1,nf
                do itype = 1,ntype
                    do ie = 1,ne
                        do iz = 1,nz
                            do jh = 1,nh                                                                
                                
                                EV(:,jh) = Vmatnext(:,iz,ie,itype,jh,bf)          ! No need to take expectations over productivity realizations
                                                                
                            end do
                            
                            do ixi=1,nxi
                                do ih =1,nh    
                                    do ia = 1,na
                                    
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
                                        if (exp_wealthtax_pension_ret == 1) then
                                            pbnow = scale_tt*tt             
                                        else
                                            pbnow = pension(ie,iz,ih,yt(tr-1,:,:)*theta_val(itheta),zgrid,penre)                           
                                        end if
                                        
                                        xinow = xigrid(ixi)*ximat(t,ie,ih)     
                                        phinow = 0.0d0 
                                        spnow = zeta(t,ie,ih)             
                                        fpast = fgrid(bf)             
                                    
                                        if (exp_wealthtax_pension_ret == 1) then
                                            ahigh = (1.0d0+r)*anow - tau_a*max(anow,0d0) + pbnow
                                        else
                                            ahigh = (1.0d0+r)*anow + pbnow            
                                        end if                                          
                                
                                        ! ************************
                                        ! (1) Value of NOT ADJUST                                        
                                        ! 2026 Update: inherited effort set to mean only if t-1 also equalised; continuation value at imposed effort 
                                        fpast = fgrid(bf)

                                        if (equal_effort_now) then

                                            ! Current effort is imposed at the benchmark age-specific mean
                                            fnow = fmeanvec(t)

                                            ! Replace inherited effort only if previous period was
                                            ! also part of the equal-effort experiment
                                            if (equal_effort_prev) then
                                                fpast = fmeanvec(t-1)
                                            end if

                                            pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)

                                            ! Current imposed effort becomes next period's habit
                                            ! No productivity uncertainty in retirement
                                            call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)

                                            call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)

                                            V_nadj(ia,iz,ie,itype,ih,bf,ixi,:,t) = -objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)

                                        else

                                            fnow = fgrid(bf)

                                            pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)

                                            call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)

                                            V_nadj(ia,iz,ie,itype,ih,bf,ixi,:,t) = -objneg_w(aptemp,ahigh,agrid,EV,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)

                                        end if
                                        
                                        ap_nadj(ia,iz,ie,itype,ih,bf,ixi,:,t) = aptemp
                                    
                                        if (ahigh - aptemp .le. 0d0) then
                                            write(*,'(a)') ' Negative Consumption: Ret, no-adjust'
                                        endif
                                        
                                    
                                        ! ************************
                                        ! (2) Value of ADJUST (w\o adjustment costs)                                        
                                        ! 2026 Update: no adjustment when efforts are equalized 
                                        if (equal_effort_now) then
                                            V_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = -1d+20                                            
                                            ap_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = 0.0d0
                                            f_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = 0.0d0
                                        else
                                            
                                            if (puregrid==1) then
                                                ! find effort using pure grid method
                                                aptempf = -1d+20
                                                vtempf = -1d+20
                                                do indf = 1,nfnow
                                                    fnow = fnowgrid(indf)                                                    
                                                    pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)                                                                                                        
                                                    ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                    call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)
                                                    call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)   ! not working                                    
                                                    vtempf(indf) = -objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)  ! not working
                                                    aptempf(indf) = aptemp
                                                end do
                                                maxind = maxloc(vtempf,1)
                                                
                                                V_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = vtempf(maxind) 
                                                ap_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = aptempf(maxind)
                                                f_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = fnowgrid(maxind)
                                                
                                                if (ahigh - aptempf(maxind) .le. 0d0) then
                                                    write(*,'(a)') ' Negative Consumption: Ret, adjust'
                                                end if
                                                
                                            elseif (puregrid==2) then
                                                ! golden section search over effort
                                                effa = fmin
                                                effb = fmax
                                                
                                                effc = effa + gsnum*(effb-effa)
                                                effd = effa + (1.0d0-gsnum)*(effb-effa)
                                                
                                                fnow = effc
                                                pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)    
                                                ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)
                                                call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)   ! not working                                    
                                                fc = objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)
                                                
                                                fnow = effd
                                                pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)                                                     
                                                ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)
                                                call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)   ! not working                                    
                                                fd = objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)
                                                
                                                diff = 1d+4
                                                iter = 0
                                                
                                                do while (diff>tol) 
                                                    if (fc>=fd) then
                                                        effz = effc + (1.0d0-gsnum)*(effb-effc)
                                                        effa = effc
                                                        effc = effd
                                                        fc = fd
                                                        effd = effz
                                                        fnow = effd
                                                        pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)                                                             
                                                        ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                        call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)
                                                        call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)   ! not working                                    
                                                        fd = objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)
                                                    else
                                                        effz = effa + gsnum*(effd-effa)
                                                        effb = effd
                                                        effd = effc
                                                        fd = fc
                                                        effc = effz
                                                        fnow = effc
                                                        pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)     
                                                        ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                        call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)
                                                        call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)   ! not working                                    
                                                        fc = objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)
                                                    end if
                                                    
                                                    diff = effd-effc
                                                    iter = iter+1
                                                end do
                                                
                                                fnow = (effc+effd)/2.0
                                                pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)     
                                                ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,.false.)
                                                call gss(amin,ahigh,aptemp,tol**2,0,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)
                                                V_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = - objneg_w(aptemp,ahigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hwgrid(1),mincon,haddft)
                                                
                                                ap_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = aptemp
                                                f_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) = fnow
                                                
                                            end if
                                                                                        
                                        end if                                    
                                                   
                                    end do   
                                end do   
                            end do                            
                        end do   
                    end do   
                end do   
            end do   
            !$omp end parallel do
            
            ! Compare values of adjusting with non-adjusting while accounting for fixed cost draw            
            do bf = 1,nf
                do ih = 1,nh
                    do itype = 1,ntype
                        do ie = 1,ne
                            do iz = 1,nz
                                do ia = 1,na
                                    fpast = fgrid(bf)
                                    Vvaltemp = 0.0d0
                                    
                                    do ixi = 1,nxi
                                        
                                        V_nadjtemp = V_nadj(ia,iz,ie,itype,ih,bf,ixi,1,t)
                                    
                                        do ichi=1,nchi
                                            
                                            V_adjtemp = V_adj(ia,iz,ie,itype,ih,bf,ixi,1,t) - chigrid(ichi,t)
                                            
                                            Mvaltemp = max(V_adjtemp, V_nadjtemp)
                                          
                                            Vvaltemp = Vvaltemp + prxi(ixi)*(1.0d0/dble(nchi))*Mvaltemp
                                        
                                        end do
                                        
                                    end do
                                    
                                    Vmat(ia,iz,ie,itype,ih,bf,t) = Vvaltemp                                    
                                    
                                end do
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
            
            !$omp parallel do default(shared) private(iz,ie,ib,itheta,ieta,itype,bf,indf,aptempf,vtempf,maxind,jh,ih,EV,EV_adj,ia,beta,anow,pbnow,xinow,phinow,spnow,fnow,fpast,pinow,aptemp,ynow,ynow_net,anhigh,awhigh,mucontemp,hours,ihw,ixi,effa,effb,effc,effd,fc,fd,effz,diff,iter)         
            do bf = 1,nf
                do itype = 1,ntype
                    do ie = 1,ne
                        do iz = 1,nz
                            do jh = 1,nh
                                
                                if (t == tr-1) then                 ! No expectations over productivity prior to retirement 
                                    EV(:,jh) = Vmatnext(:,iz,ie,itype,jh,bf)                                                 
                                    
                                else                                ! Compute conditional expected value functions on grid points
               
                                    EV(:,jh) = 0.0d0                
                                    do ja = 1,na
                                        do jz = 1,nz
                                            EV(ja,jh) = EV(ja,jh) + pz(iz,jz)*Vmatnext(ja,jz,ie,itype,jh,bf)
                                        end do
                                    end do
                                    
                                end if 
                                
                            end do
                            
                            do ixi=1,nxi
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
                                        xinow = xigrid(ixi)*ximat(t,ie,ih) 
                                        phinow = phimat(t,ie,ih)                                                     
                                        fpast = fgrid(bf)
                                    
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
                                            ! 2026 Update: inherited effort set to mean only if t-1 also equalised; continuation value at imposed effort 
                                            fpast = fgrid(bf)

                                            if (equal_effort_now) then

                                                ! Current effort is imposed at the benchmark age-specific mean
                                                fnow = fmeanvec(t)

                                                ! Only overwrite lagged effort if previous period
                                                ! was ALSO part of the equal-effort experiment
                                                if (equal_effort_prev) then
                                                    fpast = fmeanvec(t-1)
                                                end if

                                                pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)

                                                ! Current imposed effort becomes next period's habit
                                                call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)

                                                call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)

                                                V_nadj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = -objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)

                                            else

                                                fnow = fgrid(bf)

                                                pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)

                                                call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)

                                                V_nadj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = -objneg_w(aptemp,awhigh,agrid,EV,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)

                                            end if
                                            
                                            
                                            ap_nadj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = aptemp
                                        
                                            ! ************************ 
                                            ! (2) Values of ADJUST
                                            ! 2026 Update: no adjustment when efforts are equalized 
                                            if (equal_effort_now) then                                                
                                                V_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = -1d+20
                                                ap_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = 0.0d0
                                                f_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = 0.0d0    

                                            else
                                                
                                                if (puregrid==1) then
                                                    ! find effort using pure grid method
                                                    aptempf = -1d+20
                                                    vtempf = -1d+20
                                                    
                                                    do indf = 1,nfnow
                                                        fnow = fnowgrid(indf)
                                                        pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)
                                                        ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                        call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)
                                                        call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)  
                                                        vtempf(indf) = -objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) 
                                                        aptempf(indf) = aptemp
                                                    end do
                                                    maxind = maxloc(vtempf,1)
                                                        
                                                    V_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = vtempf(maxind)
                                                    ap_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = aptempf(maxind)
                                                    f_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = fnowgrid(maxind)
                                                    
                                                elseif (puregrid==2) then
                                                    ! golden section search over effort
                                                    effa = fmin
                                                    effb = fmax
                                                    
                                                    effc = effa + gsnum*(effb-effa)
                                                    effd = effa + (1.0d0-gsnum)*(effb-effa)
                                                    
                                                    fnow = effc
                                                    pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)                                                    
                                                    ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                    call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)
                                                    call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)  
                                                    fc = objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) 
                                                    
                                                    fnow = effd
                                                    pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)
                                                    ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                    call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)
                                                    call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)  
                                                    fd = objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) 
                                                    
                                                    diff = 1d+4
                                                    iter = 0
                                                    
                                                    do while (diff>tol)
                                                        if (fc>=fd) then
                                                            effz = effc + (1.0d0-gsnum)*(effb-effc)
                                                            effa = effc
                                                            effc = effd
                                                            fc = fd
                                                            effd = effz
                                                            fnow = effd
                                                            pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)  
                                                            ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                            call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)
                                                            call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)  
                                                            fd = objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) 
                                                        else
                                                            effz = effa + gsnum*(effd-effa)
                                                            effb = effd
                                                            effd = effc
                                                            fd = fc
                                                            effc = effz
                                                            fnow = effc
                                                            pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)    
                                                            ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                            call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)
                                                            call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)  
                                                            fc = objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) 
                                                        end if
                                                        
                                                        diff = effd-effc
                                                        iter = iter+1
                                                    end do
                                                    
                                                    fnow = (effc+effd)/2.0
                                                    pinow = tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)                                                       
                                                    ! 2026 Update: continuation value evaluated at candidate effort fnow
                                                    call calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,t<tr-1)
                                                    call gss(amin,awhigh,aptemp,tol**2,1,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)
                                                    
                                                    V_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = - objneg_w(aptemp,awhigh,agrid,EV_adj,beta,mucontemp,bbval,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) 
                                                    ap_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = aptemp
                                                    f_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) = fnow
                                                end if
                                                                                                                                
                                            end if
                                        end do
                                    
                                    end do
                                end do 
                            end do                            
                        end do   
                    end do 
                end do 
            end do  
            !$omp end parallel do
         
            do bf = 1,nf
                do ih = 1,nh
                    do itype = 1,ntype
                        do ie = 1,ne
                            do iz = 1,nz
                                do ia = 1,na
                                    fpast = fgrid(bf)
                                    
                                    Vvaltemp = 0.0d0
                                    
                                    do ixi=1,nxi
                                        do ichi = 1,nchi                                        
                                            Vvec_nadj = 0.0d0
                                            Vvec_adj = 0.0d0
                                            
                                            do ihw = 1,nw
                                                
                                                Vvec_nadj(ihw) = V_nadj(ia,iz,ie,itype,ih,bf,ixi,ihw,t)                                                                                        
                                                Vvec_adj(ihw) = V_adj(ia,iz,ie,itype,ih,bf,ixi,ihw,t) - chigrid(ichi,t)
                                            
                                            end do
                                            
                                                                                    
                                            V_adjtemp = maxval(Vvec_adj)                                            
                                            V_nadjtemp = maxval(Vvec_nadj)
                                        
                                            Mvaltemp = max(V_adjtemp,V_nadjtemp)
                                            
                                            Vvaltemp = Vvaltemp + prxi(ixi)*(1.0d0/dble(nchi))*Mvaltemp
                                        
                                        end do
                                    end do                                    
                                    
                                    Vmat(ia,iz,ie,itype,ih,bf,t) = Vvaltemp    
                                    
                                end do
                            end do
                        end do
                    end do
                end do
            end do
                    
        end if  
    
    end do   
    
    
    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    
    write(*,*) '    Values done'
    write(*,"('  Elasped time = ', F10.5)") eptime
  
    ! ************************************* SIMULATION ***************************************************
    call system_clock(ct1, cr)    
    
    allocate(educ(ns), dpos(ns), thetapos(ns), etapos(ns), health(ns,nt), work(ns,tr-1), survival(ns,nt), chipos(ns,nt), adj(ns,nt), spike10(ns,nt), spike20(ns,nt), xipos(ns,nt)) 
    allocate(assets(ns,nt), effort(ns,nt), earn(ns,tr-1), savings(ns), consum(ns,nt), bequests(ns,nt), zpos(ns,nt), effhabit(ns,nt)) 

    call simulation(V_adj, V_nadj, ap_adj, ap_nadj, f_adj, prxi, hwgrid, agrid, EV, pi0, l1, l2, fgrid, zeta, chigrid, yt, zgrid, pz, penre, &
        assets, effort, earn, consum, bequests,  effhabit, agedeath, work, educ, dpos, thetapos, etapos, zpos, health, survival, pensionvec, indmat, tt, mincon, adj, spike10, spike20, fmeanvec, avrgearn, xipos, specs)  


    if (recordsimdata==1) then
        ! Save simulated data for counterfactual exercises
            
        open(30, file="healthBM.txt")
        open(31, file="assetsBM.txt")
        open(33, file="effortBM.txt")    
        open(34, file="workBM.txt")
        open(38, file="survivalBM.txt")
    
        do t=1,nt
            do is=1,ns
                if (t<tr) then
                    write(34, *) work(is,t) 
                end if
                write(30, *) health(is,t)
                write(31, *) assets(is,t)
                write(33, *) effort(is,t)
                write(38, *) survival(is,t)
            end do
        end do  
    
        close(30)
        close(31)
        close(33)
        close(34)
        close(38)
         
    end if
    
    write(*,*) ' Calculating Model Statistics'
    ! ********************************* Calculate Model Statistics ******************************

    ! Inequality Decomposition
    call inequ_decomp(assets, effort, earn, consum, effhabit, agedeath, educ, dpos, thetapos, etapos, health, survival, zpos, earn_contr, hy_contr, le_contr, fh_contr, we_contr)
    
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
    
    norm = winit(2) ! winit(2) is normalized to be 1 in the model
    
    
    ! Average of wealth, consumption, earnings, work, effort by health and 2-year period (only survivors)
    do ih = 1,nh
        do t=1,nt
            denom = sum(survival(:,t), mask=health(:,t)==ih)
            avwealth2y(t,ih) = sum(assets(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom
            avcons2y(t,ih) = sum(consum(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom
            aveff2y(t,ih) = sum(effort(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom
            aveff2yall(t) = sum(effort(:,t), mask=survival(:,t)==1) / sum(survival(:,t))
            if (t<tr) then
                denomwo = sum(work(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih)
                avwork2y(t,ih) = dble(denomwo) / denom            
                avearn2y(t,ih) = sum(earn(:,t), mask=survival(:,t)==1 .AND. health(:,t)==ih) / denom      ! avr earn across health type including zeros 31.08.2021
            end if
        
        end do
    end do

    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages by period of HEALHTY (in thousand euros)'
    write(*, '(a)') ' Age Group  Avr.Wealth  Avr.Cons  Avr.Effort  Empl.Rate  Avr.Earnings (all)'
    do t=1,tr-1
        write(*, '(x,i2,a,6f10.4)') t, '     |  ', avwealth2y(t,1)*norm/1000, avcons2y(t,1)*norm/1000, aveff2y(t,1), avwork2y(t,1), avearn2y(t,1)*norm/1000
    end do
    do t=tr,nt
        write(*, '(x,i2,a,3f10.4)') t, '     |  ', avwealth2y(t,1)*norm/1000, avcons2y(t,1)*norm/1000, aveff2y(t,1)
    end do
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages by period of UNHEALHTY (in thousand euros)'
    write(*, '(a)') ' Age Group  Avr.Wealth  Avr.Cons Avr.Effort  Empl.Share  Avr.Earnings (all)'
    do t=1,tr-1
        write(*, '(x,i2,a,6f10.4)') t, '     |  ', avwealth2y(t,2)*norm/1000, avcons2y(t,2)*norm/1000, aveff2y(t,2), avwork2y(t,2), avearn2y(t,2)*norm/1000
    end do
    do t=tr,nt
        write(*, '(x,i2,a,3f10.4)') t, '     |  ', avwealth2y(t,2)*norm/1000, avcons2y(t,2)*norm/1000, aveff2y(t,2)
    end do
    
    ! Calculate share of people alive in each period
    write(*,'(x)')
    write(*,'(a)') ' Share of people alive by Age'
    do t=1,nt
        write(*, '(x,i2,a,f10.4)') t, '     |  ', sum(survival(:,t)) / dble(ns)
    enddo
    
    
    if (write_fmeanvec == 1) then
    open(10, file="fmeanvec.txt")
    do i = 1,nt
        write(10,*) aveff2yall(i)
    end do
    close(10)
    endif
    
    
    ! 15.09.2021 Average Effort and Employment by Education
    num1 = 0.0d0
    num2 = 0.0d0
    num3 = 0.0d0
    num4 = 0.0d0
    denom1 = 0
    denom2 = 0
    denom3 = 0
    denom4 = 0
    do is = 1,ns
        if (educ(is) == 1) then
            do t=1,nt
                if (survival(is,t) == 1) then
                num1 = num1 + effort(is,t)
                denom1 = denom1 + 1
                endif
            enddo
            
            do t=1,(tr-1)
                if (survival(is,t) == 1) then
                num3 = num3 + work(is,t)
                denom3 = denom3 + 1
                endif
            enddo
            
        else
            do t=1,nt
                if (survival(is,t) == 1) then
                num2 = num2 + effort(is,t)
                denom2 = denom2 + 1
                endif
            enddo
            
            do t=1,(tr-1)
                if (survival(is,t) == 1) then
                num4 = num4 + work(is,t)
                denom4 = denom4 + 1
                endif
            enddo
            
        end if
    enddo
    
    avreffHS = num1 / dble(denom1)
    avreffCL = num2 / dble(denom2)
    avrempHS = num3 / dble(denom3)
    avrempCL = num4 / dble(denom4)
    avreffAll = sum(effort, mask=survival==1) / dble(sum(survival))
    
    write(*, '(x)')
    write(*, '(x,a)') '                     Model  |   Data'
    write(*, '(x,a,2f10.4)') ' Avrg. Effort among HS ', avreffHS, 0.6868d0
    write(*, '(x,a,2f10.4)') ' Avrg. Effort among CL', avreffCL, 0.7661d0

    write(*, '(x,a,2f10.4)') ' Avrg. Emp among HS ', avrempHS, 0.6537d0
    write(*, '(x,a,2f10.4)') ' Avrg. Emp among CL', avrempCL, 0.8089d0

    write(*, '(x)')
    write(*, '(x,a,2f10.4)') ' Avrg. Effort All ', avreffAll, 0.7091d0
    
    
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
            endif   
            long_avr_eff10(j,ih) = sum(effort(:,s:e), mask=survival(:,s:e)==1 .AND. health(:,s:e)==ih) / sum(survival(:,s:e), mask=health(:,s:e)==ih)
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
        if (j==4) then
            
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
            efqu(j,i) = sum(effort(:,s:e), mask=wequart(:,s:e)==i) / dble(sum(survival(:,s:e), mask=(wequart(:,s:e)==i)))
        enddo
        
        ! --- Moments by Education
        do ie=1,ne
            num = 0
            denom = 0
            numinc = 0
            numwo = 0
            numf = 0
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
                        numf = numf + effort(is,t) 
                    endif
                enddo
            enddo
            
            long_avr_health10(j,ie) = dble(num) / dble(denom)
            long_avr_eff10_edu(j,ie) = dble(numf) / dble(denom)
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
            numf = 0
            do t=s,e
                do is=1,ns
                    if (survival(is,t) == 1 .and. educ(is) == ie .and. health(is,t) == ih) then
                       denom = denom +1
                        if (j<5) then
                        numinc = numinc + earn(is,t)
                        numwo = numwo + work(is,t)
                        endif
                        numf = numf + effort(is,t) 
                    endif
                enddo
            enddo
            
            long_avr_eff10_health_edu(j,ih,ie) = dble(numf) / dble(denom)
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
        
        do i = 1,4
            efqu_educ(j,i,ie) = sum(effort(:,s:e), mask=wequart(:,s:e)==i) / dble(sum(survival(:,s:e), mask=(wequart(:,s:e)==i)))
            asqu(j,i,ie) = ( sum(assets(:,s:e), mask=wequart(:,s:e)==i) / dble(sum(survival(:,s:e), mask=(wequart(:,s:e)==i)))) *norm
        enddo
        
        enddo ! ie
    enddo   ! j
    
    
    write(*, '(x)')
    write(*, '(a)') ' Average Health Effort in Model by Wealth Quantile over long age groups and education'
    write(*, '(a)') ' Age Group               HS                               CL'
    write(*, '(a)') ' Age Group    1st      2nd       3rd     4th      1st      2nd       3rd     4th '
    do j = 1,3
        write(*, '(x,i2,a,8f10.4)') j, '     |  ',  efqu_educ(j,1,1), efqu_educ(j,2,1), efqu_educ(j,3,1), efqu_educ(j,4,1), efqu_educ(j,1,2), efqu_educ(j,2,2), efqu_educ(j,3,2), efqu_educ(j,4,2)     
    enddo
    
    write(*, '(x)')
    write(*, '(a)') ' Average Assets in Model by Wealth Quantile over long age groups and education'
    write(*, '(a)') ' Age Group               HS                               CL'
    write(*, '(a)') ' Age Group    1st      2nd       3rd     4th      1st      2nd       3rd     4th '
    do j = 1,3
        write(*, '(x,i2,a,12f10.1)') j, '     |  ',  asqu(j,1,1), asqu(j,2,1), asqu(j,3,1), asqu(j,4,1), asqu(j,1,2), asqu(j,2,2), asqu(j,3,2), asqu(j,4,2)     
    enddo
    
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
    write(*, '(a)') ' Age Group Avr.Eff(H) Avr.Eff(U) Cons(H) Cons(U) '
     do j = 1,7
        write(*, '(x,i2,a,4f10.4)') j, '     |  ', long_avr_eff10(j,1), long_avr_eff10(j,2), long_avr_cons10(j,1)*norm/1000, long_avr_cons10(j,2)*norm/1000
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
    write(*, '(a)') '                   Share in Bad Health           |          Average Effort'
    write(*, '(a)') ' Age           1st     2nd     3rd     4th       |       1st      2nd     3rd     4th'
    do j=1,7
        write(*, '(x,i2,a,4f10.3,a,4f10.3)') j, '     |  ',  hequ(j,1), hequ(j,2), hequ(j,3), hequ(j,4), ' | ', efqu(j,1), efqu(j,2), efqu(j,3), efqu(j,4)
    end do
            
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages (7 age groups groups) of HS and CL'
    write(*, '(a)') ' Age Group  Effort(HS) Effort(CL) EmplShare(HS) EmplShare(CL) Avr.Earn(HS) Avr.Earn(CL)'
     do j = 1,7
         if (j<5) then
            write(*, '(x,i2,a,6f10.4)') j, '     |  ', long_avr_eff10_edu(j,1), long_avr_eff10_edu(j,2), long_avr_work10_edu(j,1), long_avr_work10_edu(j,2), long_avr_earn10_edu(j,1)*norm/1000, long_avr_earn10_edu(j,2)*norm/1000
         else
            write(*, '(x,i2,a,2f10.4)') j, '     |  ', long_avr_eff10_edu(j,1), long_avr_eff10_edu(j,2)
         endif
     enddo
         
    
      
    ! Calculate average based on longer age groups
    do ih = 1,nh
        long_avr_earn(1,ih) = sum(earn(:,1:tj-1), mask=survival(:,1:tj-1)==1 .AND. health(:,1:tj-1)==ih) / sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih)     ! Average earnings across everyone from that health group
        long_avr_work(1,ih) = sum(work(:,1:tj-1), mask=health(:,1:tj-1)==ih .AND. survival(:,1:tj-1)==1) / dble(sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih))     ! This is now average hours worked can be 0, 0.5, or 1 (for 2 years of full time work)
        long_avr_work_part(1,ih) = sum(survival(:,1:tj-1), mask=(health(:,1:tj-1)==ih) .AND. (work(:,1:tj-1)==0.5)) / dble(sum(survival(:,1:tj-1), mask=(health(:,1:tj-1)==ih) .and. (work(:,1:tj-1)>0.0d0)))
        long_avr_eff(1,ih) = sum(effort(:,1:tj-1), mask=survival(:,1:tj-1)==1 .AND. health(:,1:tj-1)==ih) / sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih)
        long_avr_cons(1,ih) = sum(consum(:,1:tj-1), mask=survival(:,1:tj-1)==1 .AND. health(:,1:tj-1)==ih) / sum(survival(:,1:tj-1), mask=health(:,1:tj-1)==ih)  
!        
        long_avr_earn(2,ih) = sum(earn(:,tj:tr-1), mask=survival(:,tj:tr-1)==1 .AND. health(:,tj:tr-1)==ih) / sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih)     ! Average earnings across everyone from that health group
        long_avr_work(2,ih) = sum(work(:,tj:tr-1), mask=health(:,tj:tr-1)==ih .AND. survival(:,tj:tr-1)==1) / dble(sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih))
        long_avr_work_part(2,ih) = sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih .AND. work(:,tj:tr-1)==0.5) / dble(sum(survival(:,tj:tr-1), mask=(health(:,tj:tr-1)==ih) .and. (work(:,tj:tr-1)>0.0d0)))
        long_avr_eff(2,ih) = sum(effort(:,tj:tr-1), mask=survival(:,tj:tr-1)==1 .AND. health(:,tj:tr-1)==ih) / sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih)
        long_avr_cons(2,ih) = sum(consum(:,tj:tr-1), mask=survival(:,tj:tr-1)==1 .AND. health(:,tj:tr-1)==ih) / sum(survival(:,tj:tr-1), mask=health(:,tj:tr-1)==ih)  
!    
        ! Average Consumption here only until age 84
        long_avr_eff(3,ih) = sum(effort(:,tr:nt-8), mask=survival(:,tr:nt-8)==1 .AND. health(:,tr:nt-8)==ih) / sum(survival(:,tr:nt-8), mask=health(:,tr:nt-8)==ih)
        long_avr_cons(3,ih) = sum(consum(:,tr:nt), mask=survival(:,tr:nt)==1 .AND. health(:,tr:nt)==ih) / sum(survival(:,tr:nt), mask=health(:,tr:nt)==ih)      
        
    enddo
    
    ! Effort adjustment at the individual level          
    long_eff_noadj(1) = dble(count(adj(:,1:tj-1)==0.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_eff_noadj(2) = dble(count(adj(:,tj:tr-1)==0.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_eff_noadj(3) = dble(count(adj(:,tr:nt-8)==0.and.survival(:,tr:nt-8)==1))/ dble(sum(survival(:,tr:nt-8)))        
        
    long_eff_pos(1) = dble(count(adj(:,1:tj-1)==2.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_eff_pos(2) = dble(count(adj(:,tj:tr-1)==2.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_eff_pos(3) = dble(count(adj(:,tr:nt-8)==2.and.survival(:,tr:nt-8)==1)) / dble(sum(survival(:,tr:nt-8)))        
        
    long_eff_neg(1) = dble(count(adj(:,1:tj-1)==1.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_eff_neg(2) = dble(count(adj(:,tj:tr-1)==1.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_eff_neg(3) = dble(count(adj(:,tr:nt-8)==1.and.survival(:,tr:nt-8)==1)) / dble(sum(survival(:,tr:nt-8)))
    
    eff_pos = dble(count(adj(:,1:nt-8)==2.and.survival(:,1:nt-8)==1)) / dble(sum(survival(:,1:nt-8)))
    eff_neg = dble(count(adj(:,1:nt-8)==1.and.survival(:,1:nt-8)==1)) / dble(sum(survival(:,1:nt-8)))
    
    pos_spike10 = dble(count(spike10(:,1:nt-8)==2.and.survival(:,1:nt-8)==1)) / dble(sum(survival(:,1:nt-8)))
    neg_spike10 = dble(count(spike10(:,1:nt-8)==1.and.survival(:,1:nt-8)==1)) / dble(sum(survival(:,1:nt-8)))
    pos_spike20 = dble(count(spike20(:,1:nt-8)==2.and.survival(:,1:nt-8)==1)) / dble(sum(survival(:,1:nt-8)))
    neg_spike20 = dble(count(spike20(:,1:nt-8)==1.and.survival(:,1:nt-8)==1)) / dble(sum(survival(:,1:nt-8)))
    
    long_pos_spike10(1) = dble(count(spike10(:,1:tj-1)==2.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_pos_spike10(2) = dble(count(spike10(:,tj:tr-1)==2.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_pos_spike10(3) = dble(count(spike10(:,tr:nt-8)==2.and.survival(:,tr:nt-8)==1)) / dble(sum(survival(:,tr:nt-8)))
    
    long_neg_spike10(1) = dble(count(spike10(:,1:tj-1)==1.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_neg_spike10(2) = dble(count(spike10(:,tj:tr-1)==1.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_neg_spike10(3) = dble(count(spike10(:,tr:nt-8)==1.and.survival(:,tr:nt-8)==1)) / dble(sum(survival(:,tr:nt-8)))
    
    long_pos_spike20(1) = dble(count(spike20(:,1:tj-1)==2.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_pos_spike20(2) = dble(count(spike20(:,tj:tr-1)==2.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_pos_spike20(3) = dble(count(spike20(:,tr:nt-8)==2.and.survival(:,tr:nt-8)==1)) / dble(sum(survival(:,tr:nt-8)))
    
    long_neg_spike20(1) = dble(count(spike20(:,1:tj-1)==1.and.survival(:,1:tj-1)==1)) / dble(sum(survival(:,1:tj-1)))
    long_neg_spike20(2) = dble(count(spike20(:,tj:tr-1)==1.and.survival(:,tj:tr-1)==1)) / dble(sum(survival(:,tj:tr-1)))
    long_neg_spike20(3) = dble(count(spike20(:,tr:nt-8)==1.and.survival(:,tr:nt-8)==1)) / dble(sum(survival(:,tr:nt-8)))              ! 2026 Update: 65-84 group restricted to ages 65-84 (Table 3) 
    
    do ih = 1,2
        avr_cons_work(ih) =  sum(consum(:,1:tr-1), mask=survival(:,1:tr-1)==1 .AND. health(:,1:tr-1)==ih) / sum(survival(:,1:tr-1), mask=health(:,1:tr-1)==ih)      
    enddo
    
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Averages (3 age groups groups) of HEALHTY and UNHEALTHY'
    write(*, '(a)') ' Age Group Avr.Eff(H) Avr.Eff(U) EmplShare(H) EmplShare(U) WorkPart(H) WorkPart(U) Avr.Earn(H) Avr.Earn(U)'
    write(*, '(x,i2,a,8f10.4)') 1, '     |  ', long_avr_eff(1,1), long_avr_eff(1,2), long_avr_work(1,1), long_avr_work(1,2), long_avr_work_part(1,1), long_avr_work_part(1,2), long_avr_earn(1,1)*norm/1000, long_avr_earn(1,2)*norm/1000
    write(*, '(x,i2,a,8f10.4)') 2, '     |  ', long_avr_eff(2,1), long_avr_eff(2,2), long_avr_work(2,1), long_avr_work(2,2), long_avr_work_part(2,1), long_avr_work_part(2,2), long_avr_earn(2,1)*norm/1000, long_avr_earn(2,2)*norm/1000
    write(*, '(x,i2,a,8f10.4)') 3, '     |  ', long_avr_eff(3,1), long_avr_eff(3,2)

    write(*, '(x)')
    write(*, '(a)') ' Effort Adjustments at the Individual Level'
    write(*, '(a)') ' Age Group    Non-Adjusters   Pos Adjusters   Neg Adjusters   Pos Spikes10   Neg Spikes10   Pos Spikes20   Neg Spikes20'
    write(*, '(x,i2,a,8f10.4)') 1, '     |  ', long_eff_noadj(1), long_eff_pos(1), long_eff_neg(1), long_pos_spike10(1), long_neg_spike10(1), long_pos_spike20(1), long_neg_spike20(1)
    write(*, '(x,i2,a,8f10.4)') 2, '     |  ', long_eff_noadj(2), long_eff_pos(2), long_eff_neg(2), long_pos_spike10(2), long_neg_spike10(2), long_pos_spike20(2), long_neg_spike20(2)
    write(*, '(x,i2,a,8f10.4)') 3, '     |  ', long_eff_noadj(3), long_eff_pos(3), long_eff_neg(3), long_pos_spike10(3), long_neg_spike10(3), long_pos_spike20(3), long_neg_spike20(3)
        
    write(*,"('  10% Positive spike share = ', F6.3,'  Negative spike shares =', F6.3)") pos_spike10, neg_spike10
    write(*,"('  20% Positive spike share = ', F6.3,'  Negative spike shares =', F6.3)") pos_spike20, neg_spike20
    write(*,"('  Uncond Positive Change share = ', F6.3,'  Uncond. Negative Change shares =', F6.3)") eff_pos, eff_neg
    
    
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
   
   write(*, '(x)')
    write(*, '(a)') ' Model Simulation Average Effort(6 age groups groups) of HEALHTY and UNHEALTHY by EDUCATION'
    write(*, '(a)') ' Age Group               HS                           CL'
    write(*, '(a)') ' Age Group  Avr.Eff(H) Avr.Eff(U)  Avr.Eff(H) Avr.Eff(U)'
    do j = 1,7
        write(*, '(x,i2,a,4f10.4)') j, '     |  ', long_avr_eff10_health_edu(j,1,1), long_avr_eff10_health_edu(j,2,1), long_avr_eff10_health_edu(j,1,2), long_avr_eff10_health_edu(j,2,2)         
    enddo
                        
    ! Effort statistics
    count_active = 1    
    do is=1,ns
        do t=1,nt-8
            
            if (survival(is,t) == 1 .and. survival(is,t+1) == 1) then
                effcomptemp(count_active,1) = effort(is,t)
                effcomptemp(count_active,2) = effort(is,t+1)
                count_active = count_active + 1
            end if   
                    
        end do        
    end do
    
    allocate(effcomp(count_active),effcomp2(count_active))
    
    effcomp = effcomptemp(1:count_active,1)    
    effcomp2 = effcomptemp(1:count_active,2)    
    
    ! Autocorrelation of Effort
    effcorr(1) = calc_r(effcomp,effcomp2)
    
    ! Standard Deviation of Effort
    stdeff = sqrt((sum(effcomp**2)-sum(effcomp)**2/size(effcomp))/(size(effcomp)-1))
    
    ! Average Change of Effort
    avrdeltaeff = sum(effcomp - effcomp2) / size(effcomp)
    
    deallocate(effcomp,effcomp2)
    
    
    write(*, '(x)')
    write(*, '(a)') ' Simulation Autcorrelations of Effort'
    write(*, '(a)') ' Lag       Corr    '
    
    write(*, '(x,i2,a,f10.4)') 1, '     |  ', effcorr(1)
    
    write(*, '(x)')
    write(*, '(a,f10.4)') ' Simulation Std.Dev of Effort: ', stdeff
    
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
    avgf = sum(effort, mask=survival==1) / denomall

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
                
                xisim(is,t) = xigrid(xipos(is,t))*ximat(t,educ(is),health(is,t)) 
                if (t<tr) then                    
                    phisim(is,t) = phimat(t,educ(is),health(is,t))  
                end if                
            end if
        end do
    end do
    
    avgkappa = sum(marg_consutil, mask=survival == 1 .and. consum >0) / sum(survival, mask=survival == 1 .and. consum >0)
    avgxi = sum(xisim, mask=survival==1) / sum(survival, mask=survival == 1 .and. consum >0)
    avgphi = sum(phisim, mask=survival(:,1:tr-1)==1) / sum(survival(:,1:tr-1))   ! 2026 Update: mask restricted to working ages to conform with phisim(ns,tr-1)
    
    ! 2a. Obtain LHS of vsly equation 
    if (sigma == 1d0) then
        vsl2y = avgkappa*(log(avgcons)) + bbval - avgphi*lsutil(avghours,gamma) - avgxi*fcost(avgf,psi)
    else
        vsl2y = avgkappa*((avgcons)**(1.0d0-sigma))/(1.0d0-sigma) + avgkappa*bbval - avgphi*lsutil(avghours,gamma) - avgxi*fcost(avgf,psi)
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
        
    ! 9-14. Average Health Effort Profile of HS & HEALTHY by 10-year age groups (6 moments)
    do g=9,14
        j = g-8
        y(g) = long_avr_eff10_health_edu(j,1,1)
    end do
    
    ! 15-20. Average Health Effort Profile of HS & UNHEALTHY by 10-year age groups (6 moments)
    do g=15,20
        j = g-14
        y(g) = long_avr_eff10_health_edu(j,2,1)
    end do
    
    ! 21-26. Average Health Effort Profile of CL & HEALTHY by 10-year age groups (6 moments)
    do g=21,26
        j = g-20
        y(g) = long_avr_eff10_health_edu(j,1,2)
    end do
    
    ! 27-32. Average Health Effort Profile of HS & UNHEALTHY by 10-year age groups (6 moments)
    do g=27,32
        j = g-26
        y(g) = long_avr_eff10_health_edu(j,2,2)
    end do
    
    
    ! 33-38. Median wealth age group 1-6, relative to norm=winit(2)
    y(33) = wequart10y_av(1,2)
    y(34) = wequart10y_av(2,2)    
    y(35) = wequart10y_av(3,2)   
    y(36) = wequart10y_av(4,2)
    y(37) = wequart10y_av(5,2)    
    y(38) = wequart10y_av(6,2)   
        
    ! 39. Education gradients in employment       
    y(39) = avrempCL/avrempHS
    
    ! 40-42. non-adjuster shares by long age groups (with val_inaction = 0.0)  
    y(40) = long_eff_noadj(1)
    y(41) = long_eff_noadj(2)
    y(42) = long_eff_noadj(3)
    
    ! 43. VSL multiple
    y(43) = vsl2y/(avgcons)   

    ! 44. Standard Deviation of Effort
    y(44) = stdeff
    
    ! 45. Wealth Gini
    y(45) = Sginiwe
    
    ! 46. c(healthy)/c(unhealthy) for number of household = 1 or 2    
    y(46) = avr_cons_work(1)/avr_cons_work(2)       
    
    ! 47-62. vector of earnings profiles in 1,000 EUR (4x4)
    do g=47,50
        ! HS-Healthy
        j = g-46
        y(g) = long_avr_earn10_health_edu(j,1,1)*norm/1000.0d0
    end do
    
    do g=51,54
        ! HS-Unhealthy
        j = g-50
        y(g) = long_avr_earn10_health_edu(j,2,1)*norm/1000.0d0
    end do
    
    do g=55,58
        ! CL-Healthy
        j = g-54
        y(g) = long_avr_earn10_health_edu(j,1,2)*norm/1000.0d0
    end do
    
    do g=59,62
        ! CL-Unhealthy
        j = g-58
        y(g) = long_avr_earn10_health_edu(j,2,2)*norm/1000.0d0
    end do
        
    ! 63. Variance of log earnings
    y(63) = Sginiy
        
    ! 64. Pension replacement rate
    y(64) = penavg
      
    
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
    
        
        write(10,*) objval
        ! parameters    
	    do g=1,nparam
            write(10,*) x_param(g)
        end do
        
        do g=1,nmom
            write(10,*) y(g)
        end do

    
    close(10)

    
    write(*, '(x)')
    write(*,"('===  Elapsed time for one evaluation (min) = ', F6.2,'    Obj value = ', F14.5,'    Eval num = ', I5,' ===')") eptime, objval, evalnum
    
    evalnum = evalnum + 1
    
    
    if (recordsimstat==1) then
     
    ! ********************************************
    ! Replication: model outputs
    
    ! ********************************************
    ! FIGURE 5
    ! ********************************************
        
    if (exper_run == 1) then    
        
    open(15, file="Figure5model.txt")
    
    ! ***********
    ! LEFT PANEL
    ! employment share among HEALTHY
    do g=1,4
        write(15,*) long_avr_work10(g,1)
    end do
    ! employment share among UNHEALTHY
    do g=1,4
        write(15,*) long_avr_work10(g,2)
    end do
    
    ! ***********
    ! RIGHT PANEL
    ! employment share among Non-CL
    do g=1,4
        write(15,*) long_avr_work10_edu(g,1)
    end do
    ! employment share among CL
    do g=1,4
        write(15,*) long_avr_work10_edu(g,2)
    end do
    
    close(15)
    
    
    
    ! ********************************************
    ! FIGURE 6
    ! ********************************************
    open(16, file="Figure6model.txt")
    
    ! ***********
    ! LEFT PANEL
    ! labor income among Non-CL and HEALTHY
    do g=1,4
        write(16,*) long_avr_earn10_health_edu(g,1,1)*norm/1000.0d0
    end do
    ! labor income among Non-CL and UNHEALTHY
    do g=1,4
        write(16,*) long_avr_earn10_health_edu(g,2,1)*norm/1000.0d0
    end do
    
    ! ***********
    ! RIGHT PANEL
    ! labor income among CL and HEALTHY
    do g=1,4
        write(16,*) long_avr_earn10_health_edu(g,1,2)*norm/1000.0d0
    end do
    ! labor income among CL and UNHEALTHY
    do g=1,4
        write(16,*) long_avr_earn10_health_edu(g,2,2)*norm/1000.0d0
    end do
    
    close(16)
    
    
    
    ! ********************************************
    ! FIGURE 7
    ! ********************************************
    open(17, file="Figure7model.txt")
    
    ! ***********
    ! LEFT PANEL
    ! health effort among Non-CL and HEALTHY
    do g=1,6
        write(17,*) long_avr_eff10_health_edu(g,1,1)
    end do
    ! health effort among Non-CL and UNHEALTHY
    do g=1,6
        write(17,*) long_avr_eff10_health_edu(g,2,1)
    end do
    
    ! ***********
    ! RIGHT PANEL
    ! health effort among CL and HEALTHY
    do g=1,6
        write(17,*) long_avr_eff10_health_edu(g,1,2)
    end do
    ! health effort among CL and UNHEALTHY
    do g=1,6
        write(17,*) long_avr_eff10_health_edu(g,2,2)
    end do
    
    close(17)
    
    
    
    ! ********************************************
    ! FIGURE 8
    ! ********************************************
    open(18, file="Figure8model.txt")
    
    ! ***********
    ! LEFT PANEL
    ! median wealth by 10-year age group
    do g=1,6
        write(18,*) wequart10y_av(g,2)*norm/1000.0d0 * 100
    end do
    
    ! ***********
    ! RIGHT PANEL
    ! health evoluation by education among Non-CL
    do g=1,6
        write(18,*) long_avr_health10(g,1)
    end do
    ! health evoluation by education among CL
    do g=1,6
        write(18,*) long_avr_health10(g,2)
    end do
    
    close(18)
    
    
    ! ********************************************
    ! Table 3
    ! ********************************************
    open(33, file="Table3model.txt")
    
    ! ***********
    ! Shares with positive changes 10%
    do g=1,3
        write(33,*) long_pos_spike10(g)
    end do
    ! shares with positive changes 20%
    do g=1,3
        write(33,*) long_pos_spike20(g)
    end do
    ! ***********
    ! Shares with negative changes 10%
    do g=1,3
        write(33,*) long_neg_spike10(g)
    end do
    ! shares with negative changes 20%
    do g=1,3
        write(33,*) long_neg_spike20(g)
    end do
    
    close(33)
    
    
    ! ********************************************
    ! FIGURE 9: LEFT PANEL ONLY
    ! ********************************************
    open(19, file="Figure9model1.txt")
    
    ! ***********
    ! median wealth among HEALTHY
    do g=2,6
        write(19,*) wequart10y(g,1,2)*norm/1000
    end do
    ! ***********
    ! median wealth among UNHEALTHY
    do g=2,6
        write(19,*) wequart10y(g,2,2)*norm/1000
    end do
    
    close(19)
    
    ! ********************************************
    ! FIGURE 10: 1
    ! ********************************************    
    open(20, file="Figure10model1.txt")
    
    ! ***********
    ! 25th percentile wealth
    do g=2,5
        write(20,*) (wequart10y(g,1,1)-wequart10y(g,2,1))/wequart10y(g,1,1) * 100
    end do
    ! ***********
    ! 50th percentile wealth
    do g=2,5
        write(20,*) (wequart10y(g,1,2)-wequart10y(g,2,2))/wequart10y(g,1,2) * 100
    end do
    ! ***********
    ! 75th percentile wealth
    do g=2,5
        write(20,*) (wequart10y(g,1,3)-wequart10y(g,2,3))/wequart10y(g,1,3) * 100
    end do
    
    close(20)
    
        ! ********************************************
    ! Table 4
    ! ********************************************        
    ! Underlying statistics that can be used to compute reported percentage changes
    open(34, file="Table4model1.txt")
    
    ! ***********
    ! Panel A: Average Effort 
    ! 1st wealth quartile
    do g=2,5
        write(34,*) efqu(g,1)
    end do
    ! 2nd wealth quartile
    do g=2,5
        write(34,*) efqu(g,2)
    end do
    ! 3rd wealth quartile
    do g=2,5
        write(34,*) efqu(g,3)
    end do
    ! 4th wealth quartile
    do g=2,5
        write(34,*) efqu(g,4)
    end do
    
    ! ***********
    ! Panel B: Share Bad Health
    ! 1st wealth quartile
    do g=2,5
        write(34,*) hequ(g,1)
    end do
    ! 2nd wealth quartile
    do g=2,5
        write(34,*) hequ(g,2)
    end do
    ! 3rd wealth quartile
    do g=2,5
        write(34,*) hequ(g,3)
    end do
    ! 4th wealth quartile
    do g=2,5
        write(34,*) hequ(g,4)
    end do
    
    ! ***********
    ! Panel C: Wealth-Health Gaps
    ! 25th percentile
    do g=2,5
        write(34,*) (wequart10y(g,1,1)-wequart10y(g,2,1))/wequart10y(g,1,1)
    end do
    ! 50th percentile
    do g=2,5
        write(34,*) (wequart10y(g,1,2)-wequart10y(g,2,2))/wequart10y(g,1,2)
    end do
    ! 75th percentile
    do g=2,5
        write(34,*) (wequart10y(g,1,3)-wequart10y(g,2,3))/wequart10y(g,1,3)
    end do
    
    close(34)

    
    ! ********************************************
    ! Table A.11
    ! ********************************************     
    open(11, file="TableA11model.txt")
    ! lifetime earnings
    write(11,*) 100.0d0*earn_contr
    ! wealth at age 65-66
    write(11,*) 100.0d0*we_contr
    ! no. healthy years
    write(11,*) 100.0d0*hy_contr
    ! share of healthy years
    write(11,*) 100.0d0*fh_contr
    
    close(11)
    
      ! ********************************************
    ! FIGURE A.7
    ! ********************************************
    open(19, file="FigureA7model.txt")
    
    ! ***********
    ! 25th Prctl. Wealth among HEALTHY
    do g=1,6
        write(19,*) wequart10y(g,1,1)*norm/1000
    end do
    ! ***********
    ! 25th Prctl. Wealth among UNHEALTHY
    do g=1,6
        write(19,*) wequart10y(g,2,1)*norm/1000
    end do
    
    ! ***********
    ! 50th Prctl. Wealth among HEALTHY
    do g=1,6
        write(19,*) wequart10y(g,1,2)*norm/1000
    end do
    ! ***********
    ! 50th Prctl. Wealth among UNHEALTHY
    do g=1,6
        write(19,*) wequart10y(g,2,2)*norm/1000
    end do
    
    ! ***********
    ! 75th Prctl. Wealth among HEALTHY
    do g=1,6
        write(19,*) wequart10y(g,1,3)*norm/1000
    end do
    ! ***********
    ! 75th Prctl. Wealth among UNHEALTHY
    do g=1,6
        write(19,*) wequart10y(g,2,3)*norm/1000
    end do
    
    close(19)
    
    
    
    ! ********************************************
    ! FIGURE A.8
    ! ********************************************
    open(19, file="FigureA8model.txt")
    
    ! ***********
    ! Left Panel: Non-CL
    ! ***********
    ! median wealth among HEALTHY
    do g=2,6
        write(19,*) wequart10y_edu_health(g,1,1,2)*norm/1000
    end do
    ! ***********
    ! median wealth among UNHEALTHY
    do g=2,6
        write(19,*) wequart10y_edu_health(g,2,1,2)*norm/1000
    end do
        
    ! ***********
    ! Right Panel: CL
    ! ***********
    ! median wealth among HEALTHY
    do g=2,6
        write(19,*) wequart10y_edu_health(g,1,2,2)*norm/1000
    end do
    ! ***********
    ! median wealth among UNHEALTHY
    do g=2,6
        write(19,*) wequart10y_edu_health(g,2,2,2)*norm/1000
    end do
    
    
    close(19)
    
    endif       ! End of Run1 ("Baseline Run")
    ! ********************************************
    
    
    ! Run 2 - 8 (Earnings Savings Channel)
    ! ********************************************
    ! FIGURE 10, FIGURE A4, and TABLE A12: 2-8 Earnings and Savings Channels
    ! ********************************************
    if (exper_run == 2) then
        ! 1. Earnings Channel: Complete        
        open(62, file="Figure10model2.txt")
    elseif (exper_run == 3) then
        ! 2. Savings Channel: Complete
        open(62, file="Figure10model3.txt")
    elseif (exper_run == 4) then
        ! 3. Earnings and Savings Channel joint
        open(62, file="FigureA4model.txt") 
    elseif  (exper_run ==5) then
        ! 4. Earnings Channel: Wage Loss Only
        open(62, file="TableA12model1.txt")
    elseif (exper_run == 6) then
        ! 5. Earnings Channel: Disutil. Only
        open(62, file="TableA12model2.txt")
    elseif (exper_run == 7) then
        ! 6. Savings Channel: Length Only
        open(62, file="TableA12model3.txt")
    elseif (exper_run == 8) then
        ! 7. Savings Channel: Quality Only
        open(62, file="TableA12model4.txt")  
    endif
    
    ! ***********
    ! 25th percentile wealth
    do g=2,5
        write(62,*) (wequart10y(g,1,1)-wequart10y(g,2,1))/wequart10y(g,1,1) * 100
    end do
    ! ***********
    ! 50th percentile wealth
    do g=2,5
        write(62,*) (wequart10y(g,1,2)-wequart10y(g,2,2))/wequart10y(g,1,2) * 100
    end do
    ! ***********
    ! 75th percentile wealth
    do g=2,5
        write(62,*) (wequart10y(g,1,3)-wequart10y(g,2,3))/wequart10y(g,1,3) * 100
    end do
    
    close(62)
    
    ! End of Runs 2 - 8 ("Earnings and Savings Channel")
    ! ********************************************
    
    
    ! Run 9 - 12 + 14-15 ("Equal Efforts")
    ! ********************************************
    ! FIGURE 11, FIGURE A5, and TABLE A13
    ! ********************************************
    if (exper_run == 9) then 
        ! 1. Equal effort always
        open(90, file="healthCF.txt")           ! Record health distribution in constant effort case
        do t=1,nt
            do is=1,ns
                write(90, *) health(is,t)
            end do
        end do    
        close(90)
        
        open(61, file="Figure11model2.txt")     

    elseif  (exper_run == 10) then
        ! 2. Equal effort: Age 25-44
        open(61, file="FigureA5model1.txt")
    elseif (exper_run==11) then
        ! 4. Equal effort: Age 45-64
        open(61, file="FigureA5model2.txt")
    elseif (exper_run==12) then
        ! 5. Equal effort: Age 65+
        open(61, file="FigureA5model3.txt")
    elseif (exper_run==14) then
        ! 6. Equal effort but BM health distribution
        open(61, file="TableA13model2.txt")             ! 2026 Update: file names swapped to match Table A.13 (run 14 = indirect, run 15 = direct effects)
    elseif (exper_run==15) then
        ! 6. BM but CF health distribution
        open(61, file="TableA13model1.txt")             ! 2026 Update: file names swapped to match Table A.13 (run 14 = indirect, run 15 = direct effects)
    end if
    
    ! ***********
    ! 25th percentile wealth
    do g=2,5
        write(61,*) (wequart10y(g,1,1)-wequart10y(g,2,1))/wequart10y(g,1,1)
    end do
    ! ***********
    ! 50th percentile wealth
    do g=2,5
        write(61,*) (wequart10y(g,1,2)-wequart10y(g,2,2))/wequart10y(g,1,2)
    end do
    ! ***********
    ! 75th percentile wealth
    do g=2,5
        write(61,*) (wequart10y(g,1,3)-wequart10y(g,2,3))/wequart10y(g,1,3)
    end do
    
    close(61)
    
    
    ! End of Runs 9-12 + 14-15("Equal Efforts")
    ! ********************************************
    
    
    ! Run 13 (Retirement Tax Counterfactual)
    ! ********************************************
    ! Table 4
    ! ********************************************        
    ! Underlying statistics that can be used to compute reported percentage changes
    if (exper_run == 13) then
        ! 2. Retirement Tax Counterfactual
        open(34, file="Table4model2.txt")
    end if
    
    ! ***********
    ! Panel A: Average Effort 
    ! 1st wealth quartile
    do g=2,5
        write(34,*) efqu(g,1)
    end do
    ! 2nd wealth quartile
    do g=2,5
        write(34,*) efqu(g,2)
    end do
    ! 3rd wealth quartile
    do g=2,5
        write(34,*) efqu(g,3)
    end do
    ! 4th wealth quartile
    do g=2,5
        write(34,*) efqu(g,4)
    end do
    
    ! ***********
    ! Panel B: Share Bad Health
    ! 1st wealth quartile
    do g=2,5
        write(34,*) hequ(g,1)
    end do
    ! 2nd wealth quartile
    do g=2,5
        write(34,*) hequ(g,2)
    end do
    ! 3rd wealth quartile
    do g=2,5
        write(34,*) hequ(g,3)
    end do
    ! 4th wealth quartile
    do g=2,5
        write(34,*) hequ(g,4)
    end do
    
    ! ***********
    ! Panel C: Wealth-Health Gaps
    ! 25th percentile
    do g=2,5
        write(34,*) (wequart10y(g,1,1)-wequart10y(g,2,1))/wequart10y(g,1,1)
    end do
    ! 50th percentile
    do g=2,5
        write(34,*) (wequart10y(g,1,2)-wequart10y(g,2,2))/wequart10y(g,1,2)
    end do
    ! 75th percentile
    do g=2,5
        write(34,*) (wequart10y(g,1,3)-wequart10y(g,2,3))/wequart10y(g,1,3)
    end do
    
    close(34)
    
    
  
    
    
   
    
    end if
    ! End of Replication outputs
    
    

    end subroutine sub_main



