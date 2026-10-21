module mod_solve
    

    use mod_parameters
    use mod_spline
    use mod_utils


implicit none

    contains  
    
    
    
! *******************************************************************************************************

    

    
subroutine gss(xmin,xmax,x,critg,flag,agrid,EV,beta,mucon,bb,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)

	real(8), intent(in) :: xmin, xmax, critg, beta, mucon, bb, sigma, spnow, pinow, phinow, psi, hours,mincon
    real(8), intent(in) :: agrid(:), EV(:,:)
    integer, intent(in) :: flag,ih,t,ie
	real(8), intent(out) :: x

	real(8) rg, a, b, c, d, z, fc, fd, crit, diff
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    
        fc = objneg_w(c,xmax,agrid,EV,beta,mucon,bb,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)
    

	d = a + (1.0d0-rg)*(b-a)
    
        fd = objneg_w(d,xmax,agrid,EV,beta,mucon,bb,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)
    

    diff = 1d+4
	iter = 0

	do while (diff>critg)

 		if (fc>=fd) then

        	z = c + (1.0d0-rg)*(b-c)
        	a = c
        	c = d
        	fc = fd
        	d = z
            fd = objneg_w(d,xmax,agrid,EV,beta,mucon,bb,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)
           
	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            fc = objneg_w(c,xmax,agrid,EV,beta,mucon,bb,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon)
            

    	end if

	    diff = d-c
    	iter = iter+1

	end do

	x = (c+d)/2.0d0


end subroutine gss




function objneg_w(ap,apmax,agrid,EV,beta,mucon,bb,sigma,spnow,pinow,phinow,psi,ih,t,ie,hours,mincon) result(f)
    
    implicit none 

    real(8), intent(in) :: ap, apmax, beta, bb, mucon, sigma, spnow, pinow, phinow, psi, hours, mincon
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
        
    f = mucon*((cnow)**(1.0d0-sigma))/(1.0d0-sigma) + mucon*bb - disutil  &
        + beta*spnow*(pinow*eval(1) +(1.0d0-pinow)*eval(2))        
    
    
    f = -f

end function objneg_w




function pension(ie,iz,ih,ytr,zgridtr,penre)     ! function for pension benefits
    
    implicit none 
    integer, intent(in) :: ie, iz, ih
    real(8), intent(in) :: ytr(:,:), zgridtr(:), penre
    real(8) :: pension
    
    pension = penre*ytr(ie,1)
    
end function




function lsutil(hours,gamma)     ! function for labor supply
    
    implicit none 
    real(8), intent(in) :: hours, gamma
    real(8) :: lsutil
    
    lsutil = (hours**(1.0d0+1.0d0/gamma)) / (1.0d0+1.0d0/gamma)
    
    
end function


    
function tr2yp(j,ie,ih,ieta)
    
    implicit none
      
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
            
    if (ieta == 2) then          ! Eta = 2 is healthy health type
        chealthtype = htype_dummy
    else
        chealthtype = 0d0
    endif
        
    y = const_healthtr + age_const(j) + ccollege + chealth + chealthtype 
        
    ! logit estimates
    tr2yp = exp(y) / (1.0d0 + exp(y))
        
end function

    
function tax(ynow,avrgearn,taul,lambda)
    
    implicit none
    real(8), intent(in) :: ynow,avrgearn,taul,lambda
    real(8):: tax
    
    tax = ynow - lambda*(ynow**(1.0d0-taul))*(avrgearn**taul)
    
end function




end module mod_solve