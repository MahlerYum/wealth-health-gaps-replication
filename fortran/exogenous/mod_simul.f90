    module mod_simul    
      
implicit none

    contains  
    
    
    
    ! *******************************************************************************************************
    subroutine simulation(V_nadj, ap_nadj, hwgrid, agrid, zeta, yt, zgrid, pz, penre, &
        assets, earn, consum, bequests, agedeath, work, educ, dpos, thetapos, etapos, zpos, health, survival, pensionvec, indmat, tt, mincon, avrgearn)

         use mod_parameters
         use mod_solve
         use mod_utils
         use mod_spline
    
        implicit none 
        
        real(8), intent(in) :: V_nadj(:,:,:,:,:,:,:), ap_nadj(:,:,:,:,:,:,:), hwgrid(:)
        real(8), intent(in) :: agrid(:), zeta(:,:,:), yt(:,:,:), zgrid(:)
        real(8), intent(in) :: pz(:,:), penre, tt, avrgearn, mincon
        integer, intent(in) :: indmat(:,:)
    
        real(8), intent(out) :: assets(:,:), earn(:,:), consum(:,:), bequests(:,:),  agedeath(:), work(:,:), pensionvec(:) 
        integer, intent(out) :: educ(:), dpos(:), thetapos(:), etapos(:), zpos(:,:), health(:,:), survival(:,:)

    
        ! Set variables
        real(8) V_adjtemp, V_w_adj_maxtemp, V_w_temp, V_n_temp, V_nadjtemp, V_w_nadj_temp, trash, objval            
        real(8) anow, aptemp, pinow, fnow_w, fnow_pw, Vvec_nadj(nw), fvec(nw)
        integer t, i, j, ix, jx, ia, ja, ik, jk, iz, jz, ib, jb, ih, jh, ie, ind, m, g, ind2, ind3, maxind, ihw, maxind_adj, maxind_nadj, itheta, ieta, itype
        integer is, count_active, temp, ind_wet       
        real(8) xss(nz), survivalBM(ns,nt), assetsBM(ns,nt), workBM(ns,tr-1), healthprobBM(nt), healthBM(ns,nt)        
        real(8) draw(ns), draw2(ns), draw3(ns), a, b, b1, b2, cumzss(nz), onedraw, cumpz(nz,nz), tempind(ns), cdfbeta(nb)
        real(8) states_mat(ntype,3), init_states(ntype*ne,3)       ! initial distribution over binary states + average health and effort
    
        
        type(VSL_STREAM_STATE) :: stream
        integer :: errcode
        
                
        earn = 0.0d0
        work = 0.0d0
        pensionvec = 0.0d0
        agedeath = 0.0d0    ! 2026 Update: initialise agedeath; survivors are identified by agedeath==0 below (life expectancy, Table A.11)

       ! Initial Distribution across fixed types + average health and average initial effort
        
        open(14, file="init_distr_2b2t2h.txt") 
            do i = 1,(ntype*ne)
                read(14, *) init_states(i,:)    
            end do                
        close(14)
                
        write(*,'(x)')
        write(*,'(a)') '  Initial States ' 
        write(*,'(a)') '  Prob.   Health' 
        do i = 1,(ntype*ne)
            write(*,'(2f10.4)') init_states(i,1), init_states(i,2)!, init_states(i,3)
        enddo
        
        errcode = vslnewstream(stream, brng, rnseed)
        
        ! DRAW of INTIAL SAMPLE to mirror covariance across binary education, preference (patience), health, 
        ! and productivity types and initial health effort and initial health        
        !           Type                Index       Variable Name
        ! TYPES:    Education:          ie          educ
        !           Discount factor:    ib          dpos
        !           Productivity type:  itheta      thetapos
        !           Health type:        ieta        etapos
        
        a = 0.0d0
        b = 1.0d0
        errcode = vdrnguniform(method, stream, ns, draw, a, b)
        errcode = vdrnguniform(method, stream, ns, draw2, a, b)


        do is = 1,ns
            if (draw(is) .le. init_states(ntype,1)) then      ! HS education (first ntype of init_states are HS)
                educ(is) = 1       
                
                states_mat = init_states(1:ntype,:)
                ! Get index of itype
                call getindex(states_mat(:,1),draw(is),ntype,itype)
                
                ! Get position of discount type, productivity type, and health type
                dpos(is) = indmat(1,itype)
                thetapos(is) = indmat(2,itype)
                etapos(is) = indmat(3,itype)
                
                ! Draw random realization of initial health conditional on binary states
                if (draw2(is) .le. states_mat(itype,2)) then
                    health(is,1) = 1                            ! Good Health
                else
                    health(is,1) = 2                            ! Bad Health
                endif
                
            else
                educ(is) = 2                                ! CL education
                
                states_mat = init_states(ntype+1:ntype*ne,:)
                ! Get index of itype
                call getindex(states_mat(:,1),draw(is),ntype,itype)
                
                ! Get position of discount type, productivity type, and health type 
                dpos(is) = indmat(1,itype)
                thetapos(is) = indmat(2,itype)
                etapos(is) = indmat(3,itype)
               
                ! Draw random realization of initial health conditional on binary states
                if (draw2(is) .le. states_mat(itype,2)) then
                    health(is,1) = 1                            ! Good Health
                else
                    health(is,1) = 2                            ! Bad Health
                endif
                
            endif
              
        enddo
                
        ! Fill in survival initial period
        survival(:,1) = 1
        assets(:,1) = 0.0d0          
        
        ! Print some statistics of initial distribution        
        write(*,'(a)') ' Correlations of initial types (rescaled consistently with data) '
        write(*,'(a)') '                 Educ    Pref    Health    Prod  '
        write(*,'(a,6f10.4)') ' Education   ', calc_r(dble(educ)-1,dble(educ)-1)
        write(*,'(a,6f10.4)') ' Pref Type   ', calc_r(dble(dpos)-1,dble(educ)-1), calc_r(dble(dpos)-1,dble(dpos)-1)
        write(*,'(a,6f10.4)') ' Health Type ', calc_r(dble(etapos)-1,dble(educ)-1), calc_r(dble(etapos)-1,dble(dpos)-1), calc_r(dble(etapos)-1,dble(etapos)-1)
        write(*,'(a,6f10.4)') ' Prod Type   ', calc_r(dble(thetapos)-1,dble(educ)-1), calc_r(dble(thetapos)-1,dble(dpos)-1), calc_r(dble(thetapos)-1,dble(etapos)-1), calc_r(dble(thetapos)-1,dble(thetapos)-1)
       

        ! Initial Productivity
        cumzss=0.0d0
        
        call markovss(pz,nz,tol,xss)
        do iz = 1,nz
            cumzss(iz)=sum(xss(1:iz))       ! take the stationary prod. distribution here
        end do
            
        
        errcode = vdrnguniform(method, stream, ns, draw, a, b)  

        do is=1,ns
            call getindex(cumzss,draw(is),nz,zpos(is,1))           
        end do
        

        ! Get cumulative sum of pz
        cumpz=0.0d0        
        do iz = 1,nz
            do jz=1,nz
                cumpz(iz,jz)=sum(pz(iz,1:jz))    
            end do
        end do
        
    
        ! Draw z over lifecycle
        do t = 1,tr-2
            errcode = vdrnguniform(method, stream, ns, draw, a, b) 
            do is=1,ns
                call getindex(cumpz(zpos(is,t),:),draw(is),nz,zpos(is,t+1))                
            end do
        end do


        do t = tr-1,nt
            zpos(:,t) = zpos(:,tr-1)
        end do
        
        
        ! Initialize bequests
        bequests = 0.0d0
        

        ! =================================================================================================
        ! ==========================
        ! Simulation Loop over age
        ! ==========================
                       
        
        do t = 1,nt
            if (mod(t,diagnum)==0 ) then
                write(*, "('  Simulation: age =', I4)") t
            end if
    
            if (t>=26) then
                j = 6
            elseif (t>=21) then
                j = 5
            elseif (t>=16) then
                j = 4
            elseif (t>=11) then
                j = 3
            elseif (t>=6) then
                j = 2
            else
                j = 1
            end if
                   
            do is=1,ns
        
                if (survival(is,t)==1) then
                    
                    ! individual states
        
                    iz = zpos(is,t)        
                    ie = educ(is)
                    ib = dpos(is)     
                    itheta = thetapos(is)
                    ieta = etapos(is)
                    
                    ! obtain itype
                    ! 2026 Update: type index ordering now identical to the solution (was reversed: ieta, itheta, ib)
                    itype = (ib-1)*ntheta*neta + (itheta-1)*neta + ieta
                    
                    ih = health(is,t)                    
                    
                    if (use_BM_health == 1) then        
                        ih = healthBM(is,t) 
                    end if
                    
                    anow = assets(is,t)
                        
                    if (t>=tr) then
                        aptemp = lini(agrid,ap_nadj(:,iz,ie,itype,ih,1,t),anow)
                                                                  
                    else    
                        ! ***********************
                        ! Working periods
                        ! ***********************
                        
                        Vvec_nadj = -1d+20                        
                        
                        ! obtain conditional values depending on work choice
                        do ihw = 1,nw
                            Vvec_nadj(ihw) = lini(agrid,V_nadj(:,iz,ie,itype,ih,ihw,t),anow)
                            
                        end do
                        
                        ! compare values                        
                        maxind_nadj = maxloc(Vvec_nadj,1)
                        
                        if (use_BM_work == 1) then
                            work(is,t) = workBM(is,t)
                        else
                            work(is,t) = hwgrid(maxind_nadj)
                        end if
                        
                        aptemp = lini(agrid,ap_nadj(:,iz,ie,itype,ih,maxind_nadj,t),anow)
                        
                        earn(is,t) = yt(t,educ(is),health(is,t))*zgrid(zpos(is,t))*theta_val(itheta) * work(is,t)
            
                    end if   ! (t<tr)
        
                    if (t<nt) then
                        assets(is,t+1) = aptemp                                      
                    end if
                
                    ! Record consumption (needed for SVL)
                    if (t<tr) then                        
                        if ((work(is,t)==0.0d0 ).and.(ih==2)) then  
                            consum(is,t) = max( anow*(1.0d0+r) + tt - aptemp, mincon)                                
                        else
                            consum(is,t) = max( anow*(1.0d0+r) + earn(is,t) - tax(earn(is,t),avrgearn,taul,lambda) - aptemp, mincon)
                        end if                           
                    else
                        ! Record pension vector for those that are alive in first retirement period only
                        if (t==tr) then
                            pensionvec(is) = pension(ie,iz,ih,yt(tr-1,:,:)*theta_val(itheta),zgrid,penre)      
                        end if                        
                        consum(is,t) = max( anow*(1.0d0+r) + pensionvec(is) - aptemp, mincon)
                    endif    
        
                endif ! survival == 1
        
                ! Record bequests
                if (t>1) then
                    if (survival(is,t) == 0 .and. survival(is,t-1) == 1) then
                        bequests(is,t) = assets(is,t)
                        agedeath(is) = t
                    elseif (t==nt) then
                        bequests(is,t) = aptemp        ! everyody automatically dies after nt
                    else
                        bequests(is,t) = 0.0d0
                    endif
                end if            
                
            end do ! is
            
            ! ------------------------------------------------------------------------------
        
            
            if (t<nt) then
                
                ! Draw random variables for survival to next period
                errcode = vdrnguniform(method, stream, ns, draw, a, b) 
                        
                do is=1,ns
                    if (survival(is,t) == 1) then ! if currently alive
                        ! 2026 Update: survival uses the individual's education educ(is) 
                        if (draw(is) < zeta(t,educ(is),health(is,t))) then    
                            survival(is,t+1) = 1
                        else
                            survival(is,t+1) = 0
                        end if
                    else
                        survival(is,t+1) = 0
                    end if
                    
                    if (use_BM_survival == 1) then
                        survival(is,t+1) = survivalBM(is,t+1) 
                    endif
                    
                
                end do
         
        
                ! Draw random variables for next periods health realization
        
                errcode = vdrnguniform(method, stream, ns, draw, a, b) 
            
                do is=1,ns  
            
                    if (survival(is,t+1) == 1) then ! if alive tomorrow

                        pinow = tr2yp(j,educ(is),health(is,t),etapos(is))
                        
                        if (draw(is) < pinow) then
                            health(is,t+1) = 1
                        else
                            health(is,t+1) = 2
                        end if
                        
                    else        ! not alive tomorrow
                        health(is,t+1) = 0            
                    end if
        
                end do
        
            end if ! (t<nt)
            
        end do ! t

    
        ! De-initialize
        errcode = vsldeletestream(stream)        

        ! Record death at t=76 if no death
        do is=1,ns
            if (agedeath(is) == 0) then
                agedeath(is) = nt+1
            endif
        end do
            
    end subroutine simulation
        
    
    subroutine inequ_decomp(assets, earn, consum, agedeath, educ, dpos, thetapos, etapos, health, survival, zpos, earn_contr, hy_contr, le_contr, fh_contr, we_contr)
    
        ! Calculate fraction of variance in lifetime outcomes that is due to initial conditions as in Hugget et al. (2011) and Lee and Seshadri (2018)
    
        use mod_parameters
        use mod_utils
    
        implicit none
    
        real(8), intent(in) :: assets(:,:), earn(:,:), consum(:,:), agedeath(:)
        integer, intent(in) :: educ(:), dpos(:), thetapos(:), etapos(:), health(:,:), survival(:,:), zpos(:,:)
    
        real(8), intent(out) :: earn_contr, hy_contr, le_contr, fh_contr, we_contr
    
        integer :: ib, ie, ih, is, t, count, fhabit_group(ns), iz, id, itheta, ieta
        real(8) :: pv_earn(ns), var_pv_earn, var_wealth_ret, sum, sum1, sum2, var_cond_pv_earn, sum1_hy, sum1_le, sum2_hy, sum2_le, sum1_fh, sum2_fh, sum1_we, sum2_we 
        real(8) :: ftercile(2), temp(ns), healthy_years(ns), years_alive(ns), sumhy, sumle, sumfh, sumwe, var_hy, var_le, var_cond_hy, var_cond_le, var_fh, var_cond_fh, var_cond_we
        real(8) :: frac_hy(ns)
        real(8), allocatable ::  cond_pv_earn(:,:,:,:,:,:), cond_hy(:,:,:,:,:,:), cond_le(:,:,:,:,:,:), cond_fh(:,:,:,:,:,:), cond_we(:,:,:,:,:,:)
    
        ! 1. PV of lifetime earnings
        do is = 1,ns
            pv_earn(is) = 0.0d0
            do t = 1,tr-1
                if (survival(is,t) == 1) then
                    pv_earn(is) = pv_earn(is) + earn(is,t)*(1/((1d0+r)**(dble(t)-1)))
                endif
            enddo
        enddo
    
        var_pv_earn = varsim(pv_earn)
        
        ! 1a. Wealth at age 65
        var_wealth_ret = varsim(assets(:,tr))
    
        ! 2. Healthy years and life expectancy (actual years not model periods)
        do is = 1,ns
            healthy_years(is) = 0
            years_alive(is) = 25        ! Start with 25 years of life
            do t=1,nt
                if (health(is,t) == 1) then
                    healthy_years(is) = healthy_years(is) + 2       ! Add 2 because actual years
                endif
                if (survival(is,t) == 1) then
                    years_alive(is) = years_alive(is) + 2
                endif         
            enddo
            frac_hy(is) = dble(healthy_years(is))/dble(years_alive(is))
        enddo
    
        var_hy = varsim(healthy_years)
        var_le = varsim(years_alive)
        var_fh = varsim(frac_hy)
    
        ! Calculate variance due to different initial conditions     
        
        allocate(cond_pv_earn(ne,nb,nh,nz,ntheta,neta),cond_hy(ne,nb,nh,nz,ntheta,neta),cond_le(ne,nb,nh,nz,ntheta,neta),cond_fh(ne,nb,nh,nz,ntheta,neta))
        allocate(cond_we(ne,nb,nh,nz,ntheta,neta))
        sum1 = 0.0d0
        sum2 = 0.0d0
        sum1_hy = 0.0d0
        sum1_le = 0.0d0
        sum2_hy = 0.0d0
        sum2_le = 0.0d0
        sum1_fh = 0.0d0
        sum2_fh = 0.0d0
        sum1_we = 0.0d0
        sum2_we = 0.0d0
    
        do ie = 1,ne
            do ib = 1,nb
                do ih = 1,nh
                    do iz = 1,nz
                    
                    do itheta = 1,ntheta
                    do ieta = 1,neta
                        
                    sum = 0.0d0
                    sumhy = 0.0d0
                    sumle = 0.0d0
                    sumfh = 0.0d0
                    sumwe = 0.0d0
                    count = 0
                
                    do is=1,ns
                        if ((educ(is) == ie) .and. (dpos(is) == ib) .and. (health(is,1) == ih) .and. (zpos(is,1) == iz) .and. (thetapos(is) == itheta) .and. etapos(is) == ieta) then
                            count = count + 1
                            sum = sum + pv_earn(is)
                            sumhy = sumhy + healthy_years(is)
                            sumle = sumle + years_alive(is)
                            sumfh = sumfh + dble(healthy_years(is))/dble(years_alive(is))
                            sumwe = sumwe + assets(is,tr)
                        endif
                    enddo
                
                    if (count .gt. 0) then
                        cond_pv_earn(ie,ib,ih,iz,itheta,ieta) = sum/dble(count)
                        cond_hy(ie,ib,ih,iz,itheta,ieta) = sumhy/dble(count)
                        cond_le(ie,ib,ih,iz,itheta,ieta) = sumle/dble(count)
                        cond_fh(ie,ib,ih,iz,itheta,ieta) = sumfh/dble(count)
                        cond_we(ie,ib,ih,iz,itheta,ieta) = sumwe/dble(count)
                    else
                        cond_pv_earn(ie,ib,ih,iz,itheta,ieta) = 0.0d0
                        cond_hy(ie,ib,ih,iz,itheta,ieta) = 0.0d0
                        cond_le(ie,ib,ih,iz,itheta,ieta) = 0.0d0
                        cond_fh(ie,ib,ih,iz,itheta,ieta) = 0.0d0
                        cond_we(ie,ib,ih,iz,itheta,ieta) = 0.0d0
                    endif
                
                    sum1 = sum1 + dble(count)*cond_pv_earn(ie,ib,ih,iz,itheta,ieta)/dble(ns)
                    sum2 = sum2 + dble(count)*(cond_pv_earn(ie,ib,ih,iz,itheta,ieta)**2.0d0)/dble(ns)
                    sum1_hy = sum1_hy + dble(count)*cond_hy(ie,ib,ih,iz,itheta,ieta)/dble(ns)
                    sum2_hy = sum2_hy + dble(count)*(cond_hy(ie,ib,ih,iz,itheta,ieta)**2.0d0)/dble(ns)
                    sum1_le = sum1_le + dble(count)*cond_le(ie,ib,ih,iz,itheta,ieta)/dble(ns)
                    sum2_le = sum2_le + dble(count)*(cond_le(ie,ib,ih,iz,itheta,ieta)**2.0d0)/dble(ns)
                    sum1_fh = sum1_fh + dble(count)*cond_fh(ie,ib,ih,iz,itheta,ieta)/dble(ns)
                    sum2_fh = sum2_fh + dble(count)*(cond_fh(ie,ib,ih,iz,itheta,ieta)**2.0d0)/dble(ns)
                    sum1_we = sum1_we + dble(count)*cond_we(ie,ib,ih,iz,itheta,ieta)/dble(ns)
                    sum2_we = sum2_we + dble(count)*(cond_we(ie,ib,ih,iz,itheta,ieta)**2.0d0)/dble(ns)
                    enddo
                    enddo
                    
                    enddo
                enddo
            enddo
        enddo
    
        var_cond_pv_earn = sum2 - sum1**2.0d0
        var_cond_hy = sum2_hy - sum1_hy**2.0d0
        var_cond_le = sum2_le - sum1_le**2.0d0
        var_cond_fh = sum2_fh - sum1_fh**2.0d0
        var_cond_we = sum2_we - sum1_we**2
    
        earn_contr = var_cond_pv_earn / var_pv_earn
        hy_contr = var_cond_hy / var_hy
        le_contr = var_cond_le / var_le
        fh_contr = var_cond_fh / var_fh
        we_contr = var_cond_we / var_wealth_ret
    
    end subroutine inequ_decomp
    
end module mod_simul