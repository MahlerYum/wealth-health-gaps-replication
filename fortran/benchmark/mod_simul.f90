    module mod_simul    
      
implicit none

    contains  
    
    
    
    ! *******************************************************************************************************
    subroutine simulation(V_adj, V_nadj, ap_adj, ap_nadj, f_adj, prxi, hwgrid, agrid, EV, pi0, l1, l2, fgrid, zeta, chigrid, yt, zgrid, pz, penre, &
        assets, effort, earn, consum, bequests,  effhabit, agedeath, work, educ, dpos, thetapos, etapos, zpos, health, survival, pensionvec, indmat, tt, mincon, adj, spike10, spike20, fmeanvec, avrgearn, xipos, specs)

         use mod_parameters
         use mod_solve
         use mod_utils
         use mod_spline
    
        implicit none 
        
        real(8), intent(in) :: V_adj(:,:,:,:,:,:,:,:,:), V_nadj(:,:,:,:,:,:,:,:,:), ap_adj(:,:,:,:,:,:,:,:,:), ap_nadj(:,:,:,:,:,:,:,:,:), f_adj(:,:,:,:,:,:,:,:,:), hwgrid(:)
        real(8), intent(in) :: agrid(:), EV(:,:), pi0(:,:,:,:), l1(:,:,:,:), l2(:,:,:,:), fgrid(:), zeta(:,:,:), chigrid(:,:), yt(:,:,:),zgrid(:)
        real(8), intent(in) :: pz(:,:), penre, tt, prxi(:), fmeanvec(:), avrgearn, mincon
        integer, intent(in) :: indmat(:,:), specs(:)
    
        real(8), intent(out) :: assets(:,:), effort(:,:), earn(:,:), consum(:,:), bequests(:,:),  effhabit(:,:), agedeath(:), work(:,:), pensionvec(:) 
        integer, intent(out) :: educ(:), dpos(:), thetapos(:), etapos(:), zpos(:,:), health(:,:), survival(:,:), adj(:,:), spike10(:,:), spike20(:,:), xipos(:,:)

    
        ! Set variables
        real(8) V_adjtemp, V_w_adj_maxtemp, V_w_temp, V_n_temp, V_nadjtemp, V_w_nadj_temp, chitemp, trash, objval            
        real(8) anow, aptemp, pinow, fnow, fpast, fnow_w, fnow_pw, Vvec_adj(nw), Vvec_nadj(nw), fvec(nw)
        integer t, i, j, indf, ix, jx, ia, ja, ik, jk, iz, jz, ib, jb, ih, jh, ie, ichi, bf, ind, m, g, ind2, ind3, maxind, ihw, maxind_adj, maxind_nadj, itheta, ieta, itype
        integer is, count_active, temp, ind_wet, chipos(ns,nt), ixi  
        real(8) xss(nz), survivalBM(ns,nt), assetsBM(ns,nt), workBM(ns,tr-1), effortBM(ns,nt), healthprobBM(nt), healthBM(ns,nt), cumchi(nchi), mean_anow, healthCF(ns,nt)                  
        real(8) draw(ns), draw2(ns), draw3(ns), a, b, b1, b2, cumzss(nz), onedraw, cumpz(nz,nz), tempind(ns), cdfbeta(nb), cdfxi(nxi) 
        real(8) states_mat(ntype,3), init_states(ntype*ne,3)       ! initial distribution over binary states + average health and effort
        
        integer recordsimdata,write_fmeanvec,exp_prod,exp_lab,equalsurv,exp_no_consp,exp_wealthtax_pension_ret,use_BM_health 
        integer use_BM_assets,use_BM_work,use_BM_survival,use_BM_effort,alleffort_young,alleffort_mid,alleffort_retire,use_CF_health  
    
        
        type(VSL_STREAM_STATE) :: stream
        integer :: errcode
        
    
        ! Set experiment specifications
!        recordsimdata = specs(1)
!        write_fmeanvec = specs(2)
!        exp_prod= specs(3)
!        exp_lab= specs(4)
!        equalsurv= specs(5)
!        exp_no_consp = specs(6)
!        exp_wealthtax_pension_ret = specs(7)
        use_BM_health = specs(8)
        use_BM_assets = specs(9)
        use_BM_work = specs(10)
        use_BM_survival = specs(11)
        use_BM_effort = specs(12)
        alleffort_young = specs(13)
        alleffort_mid= specs(14)
        alleffort_retire = specs(15)
        use_CF_health = specs(16)
        
    
        
        
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
        write(*,'(a)') '  Prob.   Health   Effort ' 
        do i = 1,(ntype*ne)
            write(*,'(3f10.4)') init_states(i,1), init_states(i,2), init_states(i,3)
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
                
                ! Fill in average health effort
                effhabit(is,1) = states_mat(itype,3)
                
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
                
                ! Fill in average health effort
                 effhabit(is,1) = states_mat(itype,3)
                
                ! Draw random realization of initial health conditional on binary states
                if (draw2(is) .le. states_mat(itype,2)) then
                    health(is,1) = 1                            ! Good Health
                else
                    health(is,1) = 2                            ! Bad Health
                endif
                
            endif
              
        enddo
                
        ! Fill in assets and survival initial period
        survival(:,1) = 1
        assets(:,1) = 0.0d0       
                
        ! Print some statistics of initial distribution
        write(*,"('  Initial Average Effort Habit cond. on HS =', f10.4)") sum(effhabit(:,1), mask=educ(:)==1) / sum(educ(:), mask=educ(:) == 1)
        write(*,"('  Initial Average Effort Habit cond. on CL =', f10.4)") sum(effhabit(:,1), mask=educ(:)==2) / (sum(educ(:), mask=educ(:) == 2)/2d0)
        write(*,"('  Initial Average Effort Habit cond. on Healthy =', f10.4)") sum(effhabit(:,1), mask=health(:,1)==1) / sum(health(:,1), mask=health(:,1) == 1)
        write(*,"('  Initial Average Effort Habit cond. on Unhealthy =', f10.4)") sum(effhabit(:,1), mask=health(:,1)==2) / (sum(health(:,1), mask=health(:,1) == 2)/2d0)
        write(*,"('  Initial Average Effort Habit cond. on HS AND Healthy =', f10.4)") sum(effhabit(:,1), mask=health(:,1)==1 .and. educ(:) == 1) / sum(health(:,1), mask=health(:,1) == 1 .and. educ(:) == 1)
     
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
        


        ! Get cumulative sum of chi draws
        cumchi(1) = 1.0d0/dble(nchi)
        do ichi=2,nchi
            cumchi(ichi) = cumchi(ichi-1) + 1.0d0/dble(nchi)
        end do
        
        ! nxi
        if (nxi>1) then       
            cdfxi=0.0d0
            do ixi = 1,nxi
                cdfxi(ixi) = sum(prxi(1:ixi))
            end do
        else
            cdfxi=1.0d0
        end if
        

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
        
        ! Draw adjustment cost chi
        do t = 1,nt
            errcode = vdrnguniform(method, stream, ns, draw, a, b) 
            do is=1,ns
                ! Get current realization of stochastic adjustment shock
                call getindex(cumchi,draw(is),nchi,chipos(is,t))
            end do
        end do
        
        ! Draw xi draw
        do t = 1,nt
            errcode = vdrnguniform(method, stream, ns, draw, a, b) 
            do is=1,ns
                ! Get current realization of xi utility shock
                call getindex(cdfxi,draw(is),nxi,xipos(is,t))
            end do
        end do
        
        
        ! Initialize bequests
        bequests = 0.0d0
        
        ! Load Baseline simulated data if necessary
        if (use_BM_survival == 1) then
            open(41, file="survivalBM.txt")
            do t=1,nt
                do is=1,ns
                    read(41, *) survivalBM(is,t)
                end do
            end do 
            close(41)
        endif
            
        if (use_BM_health == 1) then
            open(40, file="healthBM.txt")
            do t=1,nt
                do is=1,ns
                    read(40, *) healthBM(is,t)
                end do
            end do 
            close(40)            
        endif
        
        if (use_BM_assets == 1) then
            open(42, file="assetsBM.txt")
            do t=1,nt
                do is=1,ns
                    read(42, *) assetsBM(is,t)
                end do
            end do 
            close(42)            
        endif
        
        if (use_BM_work == 1) then
            open(43, file="workBM.txt")
            do t=1,tr-1
                do is=1,ns
                    read(43, *) workBM(is,t)
                end do
            end do 
            close(43)            
        endif
        
        if (use_BM_effort == 1) then
            open(44, file="effortBM.txt")
            do t=1,nt
                do is=1,ns
                    read(44, *) effortBM(is,t)
                end do
            end do 
            close(44)            
        endif
        
        if (exp_hesav == 1) then
            open(41, file="healthprobBM.txt")
            do t=1,nt
                read(41, *) healthprobBM(t)
            enddo
            close(41)
            write(*,"('  healthprobBM =', F10.5)") healthprobBM    
        endif

        if (use_CF_health == 1) then
            open(40, file="healthCF.txt")
            do t=1,nt
                do is=1,ns
                    read(40, *) healthCF(is,t)
                end do
            end do 
            close(40)            
        endif

        ! =================================================================================================
        ! ==========================
        ! Simulation Loop over age
        ! ==========================
        spike10 = 0                     
        spike20 = 0                     
        
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
                  
                    if (exp_nowe_eff == 1) then         ! Experiment where effort decision is based on average wealth
                            
                        ! Calculate average current assets in that period
                        mean_anow = sum(assets(:,t), mask=survival(:,t) == 1)/sum(survival(:,t))
                        !mean_anow = sum(assetsBM(:,t), mask=survival(:,t) == 1)/sum(survival(:,t))
                    endif
                    
            
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
                    
                    ! 2026 Update: set ih after the BM/CF health overrides so decision rules use the imposed health state
                    if (use_BM_health == 1) then        
                        health(is,t) = int(healthBM(is,t))
                    end if  
                    
                    if (use_CF_health == 1) then        
                        health(is,t) = int(healthCF(is,t))
                    end if  
                    
                    ih = health(is,t)        
                    
                    anow = assets(is,t)
                    fpast = effhabit(is,t)
                    chitemp = chigrid(chipos(is,t),t)
                    ixi = xipos(is,t)
                                        
                    if (t>=tr) then
                        ! Retirement periods
                        V_nadjtemp = bilini(agrid,fgrid,V_nadj(:,iz,ie,itype,ih,:,ixi,1,t),anow,fpast)
                        V_adjtemp = bilini(agrid,fgrid,V_adj(:,iz,ie,itype,ih,:,ixi,1,t),anow,fpast)                  
                        
                        if (exp_nowe_eff == 1) then
                             fnow = bilini(agrid,fgrid,f_adj(:,iz,ie,itype,ih,:,ixi,1,t),mean_anow,fpast)
                        else
                             fnow = bilini(agrid,fgrid,f_adj(:,iz,ie,itype,ih,:,ixi,1,t),anow,fpast)
                        endif 
                        
                        V_adjtemp = V_adjtemp - chitemp
                        
                        if (V_adjtemp > V_nadjtemp) then           
                            ! adjust    
                            if (fnow>fpast) then
                                adj(is,t) = 2
                                if (fnow/fpast>1.0d0+0.1d0) then
                                    spike10(is,t) = 2
                                end if
                                if (fnow/fpast>1.0d0+0.2d0) then
                                    spike20(is,t) = 2
                                end if
                                
                            else
                                adj(is,t) = 1                                                                
                                if (fnow/fpast<1.0d0-0.1d0) then
                                    spike10(is,t) = 1
                                end if                                
                                if (fnow/fpast<1.0d0-0.2d0) then
                                    spike20(is,t) = 1
                                end if             
                            end if
                            
                            if (exp_hesav == 1) then
                                aptemp = (1d0-healthprobBM(t))* bilini(agrid,fgrid,ap_adj(:,iz,ie,itype,1,:,ixi,1,t),anow,fpast) + healthprobBM(t)* bilini(agrid,fgrid,ap_adj(:,iz,ie,itype,2,:,ixi,1,t),anow,fpast)      
                            else
                                aptemp = bilini(agrid,fgrid,ap_adj(:,iz,ie,itype,ih,:,ixi,1,t),anow,fpast)
                            endif
                            
                            if (use_BM_effort == 1) then
                                effort(is,t) = effortBM(is,t)
                            else                  
                                effort(is,t) = fnow
                            endif
                            
                            
                        else
                            ! no adjust                                
                            adj(is,t) = 0
                            
                            if (exp_hesav == 1) then
                                aptemp = (1d0-healthprobBM(t))* bilini(agrid,fgrid,ap_nadj(:,iz,ie,itype,1,:,ixi,1,t),anow,fpast) + healthprobBM(t)* bilini(agrid,fgrid,ap_nadj(:,iz,ie,itype,2,:,ixi,1,t),anow,fpast)      
                            else
                                aptemp = bilini(agrid,fgrid,ap_nadj(:,iz,ie,itype,ih,:,ixi,1,t),anow,fpast)
                            endif
                            
                            if (alleffort_retire==1) then
                                effort(is,t) = fmeanvec(t)
                            else
                                if (use_BM_effort == 1) then
                                    effort(is,t) = effortBM(is,t)
                                else
                                    effort(is,t) = fpast
                                endif
                            end if
                        end if
                                                                  
                    else    
                        ! ***********************
                        ! Working periods
                        ! ***********************
                        
                        Vvec_adj = -1d+20
                        Vvec_nadj = -1d+20                        
                        
                        ! obtain conditional values depending on work choice
                        do ihw = 1,nw
                            
                            Vvec_adj(ihw) = bilini(agrid,fgrid,V_adj(:,iz,ie,itype,ih,:,ixi,ihw,t),anow,fpast)                                           
                            Vvec_nadj(ihw) = bilini(agrid,fgrid,V_nadj(:,iz,ie,itype,ih,:,ixi,ihw,t),anow,fpast)                                
                            
                            if (exp_nowe_eff == 1) then
                                fnow = bilini(agrid,fgrid,f_adj(:,iz,ie,itype,ih,:,ixi,ihw,t),mean_anow,fpast)
                            else
                                fnow = bilini(agrid,fgrid,f_adj(:,iz,ie,itype,ih,:,ixi,ihw,t),anow,fpast)
                            endif
                            
                            Vvec_adj(ihw) = Vvec_adj(ihw) - chitemp
                            
                        end do
                        
                        ! compare values
                        maxind_adj = maxloc(Vvec_adj,1)
                        maxind_nadj = maxloc(Vvec_nadj,1)
                        
                        if (Vvec_adj(maxind_adj)>Vvec_nadj(maxind_nadj)) then
                            ! adjust
                            
                            if (exp_nowe_eff == 1) then
                                fnow = bilini(agrid,fgrid,f_adj(:,iz,ie,itype,ih,:,ixi,maxind_adj,t),mean_anow,fpast)
                            else 
                                fnow = bilini(agrid,fgrid,f_adj(:,iz,ie,itype,ih,:,ixi,maxind_adj,t),anow,fpast)
                            endif
                            
                            if (fnow>fpast) then
                                adj(is,t) = 2    ! upward
                                if (fnow/fpast>1.0d0+0.1d0) then
                                    spike10(is,t) = 2
                                end if
                                if (fnow/fpast>1.0d0+0.2d0) then
                                    spike20(is,t) = 2
                                end if
                            else
                                adj(is,t) = 1    ! downward
                                if (fnow/fpast<1.0d0-0.1d0) then
                                    spike10(is,t) = 1
                                end if                           
                                if (fnow/fpast<1.0d0-0.2d0) then
                                    spike20(is,t) = 1
                                end if
                            end if
                            
                            if (use_BM_effort == 1) then
                                effort(is,t) = effortBM(is,t)
                            else
                                effort(is,t) = fnow
                            endif
                            
                            if (use_BM_work == 1) then
                                work(is,t) = workBM(is,t)
                            else
                                work(is,t) = hwgrid(maxind_adj)
                            endif
                            
                            if (exp_hesav == 1) then
                                aptemp = (1d0-healthprobBM(t))* bilini(agrid,fgrid,ap_adj(:,iz,ie,itype,1,:,ixi,maxind_adj,t),anow,fpast) + healthprobBM(t)* bilini(agrid,fgrid,ap_adj(:,iz,ie,itype,2,:,ixi,maxind_adj,t),anow,fpast)      
                            else
                                aptemp = bilini(agrid,fgrid,ap_adj(:,iz,ie,itype,ih,:,ixi,maxind_adj,t),anow,fpast)
                            endif
                            if (fnow>1.0d0) then
                                print *, fnow
                            end if
                            
                        else
                            ! not adjust
                            adj(is,t) = 0
                            
                            if (((alleffort_young == 1).and.(j==1 .or. j==2)).or.((alleffort_mid == 1).and.(j==3 .or. j==4))) then
                                effort(is,t) = fmeanvec(t)
                            else
                                if (use_BM_effort == 1) then
                                    effort(is,t) = effortBM(is,t)
                                else
                                    effort(is,t) = fpast    ! Effort stays constant while not working
                                endif
                                
                            end if
                            
                            if (use_BM_work == 1) then
                                work(is,t) = workBM(is,t)
                            else
                                work(is,t) = hwgrid(maxind_nadj)
                            endif
                            
                            if (exp_hesav == 1) then
                                aptemp = (1d0-healthprobBM(t))* bilini(agrid,fgrid,ap_nadj(:,iz,ie,itype,1,:,ixi,maxind_nadj,t),anow,fpast) + healthprobBM(t)* bilini(agrid,fgrid,ap_nadj(:,iz,ie,itype,2,:,ixi,maxind_nadj,t),anow,fpast)      
                            else
                                aptemp = bilini(agrid,fgrid,ap_nadj(:,iz,ie,itype,ih,:,ixi,maxind_nadj,t),anow,fpast)
                            endif
                        end if
                        
                        earn(is,t) = yt(t,educ(is),health(is,t))*zgrid(zpos(is,t))*theta_val(itheta) * work(is,t)
            
                    end if   ! (t<tr)
        
      
                    if (t<nt) then
                        if (use_BM_assets == 1) then
                            assets(is,t+1) = assetsBM(is,t+1)
                        else
                            assets(is,t+1) = aptemp              
                        endif
                    end if
                
                    ! Record consumption (needed for SVL)
                    if (t<tr) then
                        if ((work(is,t)==0.0d0 ).and.(ih==2)) then  
                            consum(is,t) = max( anow*(1.0d0+r) + tt - aptemp, mincon)                                
                        else
                            consum(is,t) = max( anow*(1.0d0+r) + earn(is,t) - tax(earn(is,t),avrgearn,taul,lambda) - aptemp, mincon)
                        end if
                    else
                        ! Record pension vector for those that are alive in first retirement period only -> calibrate that to average pension replacement rate of average lifetime earnings
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
            !!$omp end parallel do
     
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

                        effhabit(is,t+1) = effort(is,t)                    
                        pinow = tr2yp(j,effort(is,t),effhabit(is,t),educ(is),health(is,t),etapos(is),pi0,l1,l2)
                                            
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
        enddo
        
    end subroutine simulation
        
    
    subroutine inequ_decomp(assets, effort, earn, consum, effhabit, agedeath, educ, dpos, thetapos, etapos, health, survival, zpos, earn_contr, hy_contr, le_contr, fh_contr, we_contr)

        ! Calculate fraction of variance in lifetime outcomes that is due to initial conditions as in Hugget et al. (2011) and Lee and Seshadri (2018)
        ! 2026 Update: mean and variance of wealth at 65-66 computed over survivors only (Table A.11)
        use mod_parameters
        use mod_utils
    
        implicit none
    
        real(8), intent(in) :: assets(:,:), effort(:,:), earn(:,:), consum(:,:), effhabit(:,:), agedeath(:)
        integer, intent(in) :: educ(:), dpos(:), thetapos(:), etapos(:), health(:,:), survival(:,:), zpos(:,:)
    
        real(8), intent(out) :: earn_contr, hy_contr, le_contr, fh_contr, we_contr
    
        integer :: ib, ie, ih, is, t, count, fhabit_group(ns), iz, nd, id, itheta, ieta, count_we, ns_we
        real(8) :: pv_earn(ns), var_pv_earn, var_wealth_ret, sumea, sum1, sum2, var_cond_pv_earn, sum1_hy, sum1_le, sum2_hy, sum2_le, sum1_fh, sum2_fh, sum1_we, sum2_we, sumef, sum1_ef, sum2_ef
        real(8) :: ftercile(2), temp(ns), healthy_years(ns), years_alive(ns), sumhy, sumle, sumfh, sumwe, var_hy, var_le, var_cond_hy, var_cond_le, var_fh, var_cond_fh, var_cond_we, var_cond_ef, var_ef
        real(8) :: frac_hy(ns), avrg_effort(ns), ef_contr, mean_we 
        real(8), allocatable ::  cond_pv_earn(:,:,:,:,:,:,:), cond_hy(:,:,:,:,:,:,:), cond_le(:,:,:,:,:,:,:), cond_fh(:,:,:,:,:,:,:), cond_we(:,:,:,:,:,:,:), cond_ef(:,:,:,:,:,:,:)
            
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
        
       ! 1a. Wealth at age 65 (survivors only)          
        ns_we = 0
        mean_we = 0.0d0
        do is = 1,ns
            if (survival(is,tr) == 1) then
                ns_we = ns_we + 1
                mean_we = mean_we + assets(is,tr)
            endif
        enddo
        mean_we = mean_we/dble(ns_we)
        var_wealth_ret = 0.0d0
        do is = 1,ns
            if (survival(is,tr) == 1) then
                var_wealth_ret = var_wealth_ret + (assets(is,tr) - mean_we)**2
            endif
        enddo
        var_wealth_ret = var_wealth_ret/dble(ns_we)
    
        ! 2. Healthy years and life expectancy (actual years not model periods)
        do is = 1,ns
            healthy_years(is) = 0
            years_alive(is) = 0        ! Start with 25 years of life
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
        
        ! 2a. Average Lifestyle Behaviors over Individuals
        do is = 1,ns
            avrg_effort(is) = sum(effort(is,:),mask=survival(is,:) == 1) / (dble(years_alive(is))/2d0)
        enddo
        
    
        var_hy = varsim(healthy_years)
        var_le = varsim(years_alive)
        var_fh = varsim(frac_hy)
        var_ef = varsim(avrg_effort)

        ! Calculate variance due to different initial conditions
        
        ! Sort individuals into 3 equally sized groups by inital health effort habit
        nd = 3
        temp = effhabit(:,1)
        call quicksort(temp)
        ftercile(1) = temp(nint(dble(ns)/3.0d0))
        ftercile(2) = temp(nint(dble(ns)*2d0/3d0))
        do is = 1,ns
            if (effhabit(is,1) <= ftercile(1)) then
                fhabit_group(is) = 1
            elseif ((effhabit(is,1) > ftercile(1)) .and. (effhabit(is,1) <= ftercile(2))) then
                fhabit_group(is) = 2
            else
                fhabit_group(is) = 3
            endif
        enddo
    
        allocate(cond_pv_earn(ne,nb,nh,nz,nd,ntheta,neta),cond_hy(ne,nb,nh,nz,nd,ntheta,neta),cond_le(ne,nb,nh,nz,nd,ntheta,neta),cond_fh(ne,nb,nh,nz,nd,ntheta,neta))
        allocate(cond_we(ne,nb,nh,nz,nd,ntheta,neta),cond_ef(ne,nb,nh,nz,nd,ntheta,neta))
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
        sum1_ef = 0.0d0
        sum2_ef = 0.0d0
    
        do ie = 1,ne
            do ib = 1,nb
                do ih = 1,nh
                    do iz = 1,nz
                        do id = 1,nd         ! also condition on initial effort habit group
                            do itheta = 1,ntheta
                                do ieta = 1,neta
                        
                                    sumea = 0.0d0
                                    sumhy = 0.0d0
                                    sumle = 0.0d0
                                    sumfh = 0.0d0
                                    sumwe = 0.0d0
                                    sumef = 0.0d0
                                    count = 0
                                    count_we = 0
                
                                    do is=1,ns
                                        if ((educ(is) == ie) .and. (dpos(is) == ib) .and. (health(is,1) == ih) .and. (zpos(is,1) == iz) .and. (fhabit_group(is) ==id) .and. (thetapos(is) == itheta) .and. etapos(is) == ieta) then
                                            count = count + 1
                                            sumea = sumea + pv_earn(is)
                                            sumhy = sumhy + healthy_years(is)
                                            sumle = sumle + years_alive(is)
                                            sumfh = sumfh + dble(healthy_years(is))/dble(years_alive(is))
                                              if (survival(is,tr) == 1) then
                                                  count_we = count_we + 1
                                                  sumwe = sumwe + assets(is,tr)
                                              endif
                                            sumef = sumef + avrg_effort(is)
                                        endif
                                    enddo
                
                                    ! Calculate conditional mean for each state at age 25
                                    if (count .gt. 0) then
                                        cond_pv_earn(ie,ib,ih,iz,id,itheta,ieta) = sumea/dble(count)
                                        cond_hy(ie,ib,ih,iz,id,itheta,ieta) = sumhy/dble(count)
                                        cond_le(ie,ib,ih,iz,id,itheta,ieta) = sumle/dble(count)
                                        cond_fh(ie,ib,ih,iz,id,itheta,ieta) = sumfh/dble(count)
                                          if (count_we .gt. 0) then
                                              cond_we(ie,ib,ih,iz,id,itheta,ieta) = sumwe/dble(count_we)
                                          else
                                              cond_we(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                          end if
                                        cond_ef(ie,ib,ih,iz,id,itheta,ieta) = sumef/dble(count)
                                    else
                                        cond_pv_earn(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                        cond_hy(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                        cond_le(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                        cond_fh(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                        cond_we(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                        cond_ef(ie,ib,ih,iz,id,itheta,ieta) = 0.0d0
                                    end if
 
                                    sum1 = sum1 + dble(count)*cond_pv_earn(ie,ib,ih,iz,id,itheta,ieta)/dble(ns)
                                    sum2 = sum2 + dble(count)*(cond_pv_earn(ie,ib,ih,iz,id,itheta,ieta)**2.0d0)/dble(ns)
                                    sum1_hy = sum1_hy + dble(count)*cond_hy(ie,ib,ih,iz,id,itheta,ieta)/dble(ns)
                                    sum2_hy = sum2_hy + dble(count)*(cond_hy(ie,ib,ih,iz,id,itheta,ieta)**2.0d0)/dble(ns)
                                    sum1_le = sum1_le + dble(count)*cond_le(ie,ib,ih,iz,id,itheta,ieta)/dble(ns)
                                    sum2_le = sum2_le + dble(count)*(cond_le(ie,ib,ih,iz,id,itheta,ieta)**2.0d0)/dble(ns)
                                    sum1_fh = sum1_fh + dble(count)*cond_fh(ie,ib,ih,iz,id,itheta,ieta)/dble(ns)
                                    sum2_fh = sum2_fh + dble(count)*(cond_fh(ie,ib,ih,iz,id,itheta,ieta)**2.0d0)/dble(ns)
                                    sum1_we = sum1_we + dble(count_we)*cond_we(ie,ib,ih,iz,id,itheta,ieta)/dble(ns_we)
                                    sum2_we = sum2_we + dble(count_we)*(cond_we(ie,ib,ih,iz,id,itheta,ieta)**2.0d0)/dble(ns_we)
                                    sum1_ef = sum1_ef + dble(count)*cond_ef(ie,ib,ih,iz,id,itheta,ieta)/dble(ns)
                                    sum2_ef = sum2_ef + dble(count)*(cond_ef(ie,ib,ih,iz,id,itheta,ieta)**2.0d0)/dble(ns)
                                enddo
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
        var_cond_ef = sum2_ef - sum1_ef**2
  
        earn_contr = var_cond_pv_earn / var_pv_earn
        hy_contr = var_cond_hy / var_hy
        le_contr = var_cond_le / var_le
        fh_contr = var_cond_fh / var_fh
        we_contr = var_cond_we / var_wealth_ret
        ef_contr = var_cond_ef / var_ef
    
    end subroutine inequ_decomp
    
end module mod_simul