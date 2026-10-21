Subroutine transfrm(dim,dtr,x,tx)

    implicit none
    Integer, intent(in) :: dim, dtr
    Real(8), intent(in) :: x(dim)
    Real(8), intent(out) :: tx(dim)
    Integer :: i1, i2
   
    If (dtr==1) then
        
        ! for all parameters
        tx = log(x)
        
        ! for parameters that need special treatment    
        ! parameters can take negative values                
        tx(17) = log(0.01d0+x(17))
        tx(19) = log(0.01d0+x(19))
        
    else if (dtr==2) then
        
        tx = exp(x)
       
        tx(17) = -0.01d0+exp(x(17))
        tx(19) = -0.01d0+exp(x(19))
        
        
        
    end if


end subroutine transfrm