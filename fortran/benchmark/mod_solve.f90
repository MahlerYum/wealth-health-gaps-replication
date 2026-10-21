module mod_solve
    

    use mod_parameters
    use mod_spline
    use mod_utils


implicit none

    contains  
    
    
    
! *******************************************************************************************************


subroutine exp_specs(specs,exper_run)
    
    use mod_parameters
    
    implicit none
    
    integer, intent(out) :: specs(nspecs)
    integer, intent(in) :: exper_run
    
    ! Order of specs
    !recordsimdata 
    !write_fmeanvec 
    !exp_prod
    !exp_lab
    !equalsurv
    !exp_no_consp 
    !exp_wealthtax_pension_ret 
    !use_BM_health 
    !use_BM_assets 
    !use_BM_work 
    !use_BM_survival
    !use_BM_effort
    !alleffort_young 
    !alleffort_mid
    !alleffort_retire 
    !use_CF_health 

    specs(1:16) = 0
    if (exper_run == 1) then
        specs(1) = 1
        specs(2) = 1
    elseif (exper_run == 2) then
        specs(3) = 1
        specs(4) = 1
        specs(8) = 1
    elseif (exper_run == 3) then
        specs(5) = 1
        specs(6) = 1
        specs(8) = 1
        specs(11) = 1     ! 2026 Update: baseline survival realizations, as stated in footnote 41 (Fig. 10, Table A.12)
    elseif (exper_run == 4) then
        specs(3) = 1
        specs(4) = 1
        specs(5) = 1
        specs(6) = 1
        specs(8)  = 1     ! 2026 Update: Fig. A.4 holds health and survival realizations at baseline (as in Fig. 10)
        specs(11) = 1     
    elseif (exper_run == 5) then
        specs(3) = 1
        specs(8) = 1
    elseif (exper_run == 6) then
        specs(4) = 1
        specs(8) = 1
    elseif (exper_run == 7) then
        specs(5) = 1
        specs(8) = 1
        specs(11) = 1     ! 2026 Update: baseline survival realizations (Table A.12, length only)
    elseif (exper_run == 8) then
        specs(6) = 1
        specs(8) = 1
    elseif (exper_run == 9) then
        specs(13) = 1
        specs(14) = 1
        specs(15) = 1
    elseif (exper_run == 10) then
        specs(13) = 1
    elseif (exper_run == 11) then
        specs(14) = 1
    elseif (exper_run == 12) then
        specs(15) = 1
    elseif (exper_run == 13) then
        specs(7) = 1
        specs(9) = 1
        specs(10) = 1
        specs(11) = 1          ! 2026 Update: baseline survival realizations (Table 4) 
    elseif (exper_run == 14) then
        specs(8) = 1        
        specs(13) = 1
        specs(14) = 1
        specs(15) = 1
    elseif (exper_run == 15) then
        specs(16) = 1
    endif
    
        
end subroutine exp_specs
    
    
subroutine gss(xmin,xmax,x,critg,flag,agrid,EV,beta,mucon,bb,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)

	real(8), intent(in) :: xmin, xmax, critg, beta, mucon, bb, sigma, fnow, spnow, pinow, xinow, phinow, psi, hours,mincon,haddft
    real(8), intent(in) :: agrid(:), EV(:,:)
    integer, intent(in) :: flag,ih,t,ie
	real(8), intent(out) :: x

	real(8) rg, a, b, c, d, z, fc, fd, crit, diff
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    fc = objneg_w(c,xmax,agrid,EV,beta,mucon,bb,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)
    
	d = a + (1.0d0-rg)*(b-a)
    fd = objneg_w(d,xmax,agrid,EV,beta,mucon,bb,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)
    
    diff = 1d+4
	iter = 0

	do while (diff>critg)

 		if (fc>=fd) then

        	z = c + (1.0d0-rg)*(b-c)
        	a = c
        	c = d
        	fc = fd
        	d = z
            fd = objneg_w(d,xmax,agrid,EV,beta,mucon,bb,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)
            
	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            fc = objneg_w(c,xmax,agrid,EV,beta,mucon,bb,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft)
            
    	end if

	    diff = d-c
    	iter = iter+1

	end do

	x = (c+d)/2.0d0


end subroutine gss




function objneg_w(ap,apmax,agrid,EV,beta,mucon,bb,sigma,fnow,spnow,pinow,xinow,phinow,psi,ih,t,ie,hours,mincon,haddft) result(f)
    
    implicit none 

    real(8), intent(in) :: ap, apmax, beta, bb, mucon, sigma, fnow, spnow, pinow, xinow, phinow, psi, hours,mincon,haddft
    real(8), intent(in) :: agrid(:), EV(:,:)
    integer, intent(in) :: ih, t, ie
    real(8) cnow, eval(nh), edv, ed2v, f, disutil 
    integer ra, h
    ra = size(agrid,1)-2
   
    cnow = max(apmax - ap, mincon)
    
    do h = 1,nh
        eval(h) = lini(agrid,EV(:,h),ap)
    end do
    
    
    disutil = phinow*lsutil(hours,gamma)
        
    if (sigma == 1.0d0) then
        f = mucon*log(cnow) + bb - disutil - xinow*fcost(fnow,psi) &
            + beta*spnow*(pinow*eval(1) +(1.0d0-pinow)*eval(2))      
    else
        
        f = mucon*((cnow)**(1.0d0-sigma))/(1.0d0-sigma) + mucon*bb - disutil - xinow*fcost(fnow,psi) &
            + beta*spnow*(pinow*eval(1) +(1.0d0-pinow)*eval(2))        
    end if
    
    
    
    f = -f

end function objneg_w




function pension(ie,iz,ih,ytr,zgridtr,penre)     ! function for pension benefits
    
    implicit none 
    integer, intent(in) :: ie, iz, ih
    real(8), intent(in) :: ytr(:,:), zgridtr(:), penre
    real(8) :: pension
    
    pension = penre*ytr(ie,1) 
    
end function



function fcost(fnow,psi)     ! function for effort costs using constant frisch elasticity formula 
    
    implicit none 
    real(8), intent(in) :: fnow, psi
    real(8) :: fcost
    
    fcost = (fnow**(1.0d0+(1.0d0/psi)))/(1.0d0+(1.0d0/psi))            
    

end function


function lsutil(hours,gamma)     ! function for labor supply
    
    implicit none 
    real(8), intent(in) :: hours, gamma
    real(8) :: lsutil
    
    lsutil = (hours**(1.0d0+1.0d0/gamma)) / (1.0d0+1.0d0/gamma)
    
    
end function


    
function tr2yp(j,fnow,fpast,ie,ih,ieta,pi0,l1,l2)
    
    implicit none
    
    real(8), intent(in) :: fnow,fpast
    real(8), intent(in) :: pi0(:,:,:,:), l1(:,:,:,:), l2(:,:,:,:)
    integer, intent(in) :: j, ie, ih, ieta
    real(8) :: tr2yp, y, chealth, chealthtype, ccollege, trash
    
    if (ie == 1) then           ! ie = 1 is HS
        ccollege = 0.0d0
    else
        ccollege = college_dummy
    endif
        
    if (ih == 1) then           ! ih = 1 is currently healthy
        chealth = healthy_dummy
    else 
        chealth = 0d0
    endif
            
    if (ieta == 2) then          ! eta = 2 is healthy health type
        chealthtype = htype_dummy
    else
        chealthtype = 0d0
    endif
        
    y = const_healthtr + age_const(j) + ccollege + chealth + chealthtype + fnow*eff_param(1) + fpast*eff_param(2) + fnow**2 * eff_sq
        
    ! logit estimates
    tr2yp = exp(y) / (1.0d0 + exp(y))
    
end function

    
function tax(ynow,avrgearn,taul,lambda)
    
    implicit none
    real(8), intent(in) :: ynow,avrgearn,taul,lambda
    real(8):: tax
    
    tax = ynow - lambda*(ynow**(1.0d0-taul))*(avrgearn**taul)
    
end function


! 2026 Update: E[V_{j+1}] after adjustment evaluated at the newly chosen effort f_j 
subroutine calc_EV_adjust(EV_adj,Vmatnext,fgrid,pz,iz,ie,itype,fnow,prod_uncertainty)

    implicit none

    real(8), intent(out) :: EV_adj(:,:)
    real(8), intent(in)  :: Vmatnext(:,:,:,:,:,:)
    real(8), intent(in)  :: fgrid(:), pz(:,:)
    real(8), intent(in)  :: fnow
    integer, intent(in)  :: iz, ie, itype
    logical, intent(in)  :: prod_uncertainty

    integer :: jf, jz, jh
    real(8) :: wf

    ! Find habit-grid bracket only ONCE
    call findiw(fnow,fgrid,jf,wf)

    if (prod_uncertainty) then

        EV_adj = 0.0d0

        do jz = 1,size(pz,2)
            do jh = 1,size(EV_adj,2)

                EV_adj(:,jh) = EV_adj(:,jh) + pz(iz,jz) * &
                    ( wf*Vmatnext(:,jz,ie,itype,jh,jf) + &
                    (1.0d0-wf)*Vmatnext(:,jz,ie,itype,jh,jf+1) )

            end do
        end do

    else

        do jh = 1,size(EV_adj,2)

            EV_adj(:,jh) = &
                wf*Vmatnext(:,iz,ie,itype,jh,jf) + &
                (1.0d0-wf)*Vmatnext(:,iz,ie,itype,jh,jf+1)

        end do

    endif

end subroutine calc_EV_adjust

end module mod_solve