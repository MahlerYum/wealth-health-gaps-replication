module mod_utils
    
   !use mkl_vsl     ! 07.09.2021
   
implicit none

contains  
    
    
    
subroutine rouwenhorst(rho,sigma_eps,n,zgrid,P)

    ![zgrid, P] = rouwenhorst(rho, sigma_eps, n)
    !
    ! rho is the 1st order autocorrelation
    ! sigma_eps is the standard deviation of the error term
    ! n is the number of points in the discrete approximation
    !
    ! see "Finite State Markov-chain Approximations to Highly Persistent
    ! Processes"

    implicit none

    integer, intent(in) :: n
    real(8), intent(in) :: rho,sigma_eps
    real(8), dimension(n,n), intent(out) :: P
    real(8), dimension(n), intent(out) :: zgrid

    integer :: i, status
    real(8) :: mu_eps, q, nu
    real(8), allocatable, dimension(:,:) :: P11, P12, P21, P22

    mu_eps = 0.d0

    q = (rho+1.d0)/2.d0
    nu = sqrt((dble(n)-1.d0)/(1.d0-rho**2.d0)) * sigma_eps
    P = 0.d0

    P(1,1) = q
    P(1,2) = 1-q
    P(2,1) = 1-q
    P(2,2) = q

    do i=2,n-1
        allocate(P11(i+1,i+1), P12(i+1,i+1), P21(i+1,i+1), P22(i+1,i+1), STAT=status)
        P11 = 0.d0
        P12 = 0.d0
        P21 = 0.d0
        P22 = 0.d0

        P11(1:i,1:i) = P(1:i,1:i)
        P12(1:i,2:i+1) = P(1:i,1:i)
        P21(2:i+1,1:i) = P(1:i,1:i)
        P22(2:i+1,2:i+1) = P(1:i,1:i)

        P(1:i+1,1:i+1) = q*P11 + (1-q)*P12 + (1-q)*P21 + q*P22
        P(2:i,:) = P(2:i,:)/2.d0

        deallocate(P11,P12,P21,P22, STAT = status)
    end do

    zgrid = linspace(mu_eps/(1.0d0-rho)-nu,mu_eps/(1.0d0-rho)+nu,n)

end subroutine rouwenhorst



subroutine tauchen(N,mu,rho,sigma,m,Z,Zprob)

    ! declare variables
    integer, intent(in) :: N
    real(8), intent(in) :: mu, rho, sigma, m
    real(8), intent(out), dimension(N) :: Z
    real(8), intent(out), dimension(N,N) :: Zprob

    real(8) c, w
    integer i, j, k

    c = (1.0d0-rho)*mu

    Z(N) = m*sqrt(sigma**2/(1.0d0-rho**2))
    Z(1) = -Z(N)
    w = (Z(N)-Z(1))/(N-1)

    do i = 2,(N-1)
        Z(i) = Z(1) + w*(i-1)
    end do

    do i = 1,N
        Z(i) = Z(i) + mu
    end do

    do j = 1,N
        do k = 1,N
            if (k==1) then
                Zprob(j,k) = cdf_normal((Z(1) - c - rho*Z(j) + w/2.0d0)/sigma)
            else if (k==N) then
                Zprob(j,k) = 1.0d0 - cdf_normal((Z(N) - c - rho*Z(j) - w/2.0d0)/sigma)
            else
                Zprob(j,k) = cdf_normal((Z(k) - c - rho*Z(j) + w/2.0d0)/sigma) - cdf_normal((Z(k) - c - rho*Z(j) - w/2.0d0)/sigma)
            end if
        end do
    end do

    return


end subroutine tauchen


real(8) function cdf_normal(x)

    real(8), intent(in) :: x
    cdf_normal = 0.5d0 * erfc(-x/sqrt(2.0d0))

end function cdf_normal




function linspace(a,b,m) result(x)

    integer, intent(in) :: m
    real(8), intent(in) :: a, b
    real(8), dimension(m) :: x
    integer i
    x(1) = a
    do i = 2,m
        x(i) = x(i-1) + (b-a)/(m-1)
    end do

end function linspace


function logspace(a,b,m) result(x)

    integer, intent(in) :: m
    real(8), intent(in) :: a, b
    real(8), dimension(m) :: x
    integer i

    x = linspace(log(10.0d0**a),log(10.0d0**b),m)
    x = exp(x)

end function logspace



function lini(t,ft,z) result(fz)


    real(8), intent(in) :: t(:), ft(:), z
    real(8) w, fz
    integer it, nt, ind

    nt = size(t)

    if (z<=t(1)) then
        fz = ft(1)
    elseif (z>=t(Nt)) then
        fz = ft(Nt)
    else ! t(1) < z < t(nt)
        ind = 0
        do it = 1,nt
            if (t(it)>=z) exit
            ind = ind+1
        end do
        w = (t(ind+1)-z)/(t(ind+1)-t(ind))
        fz = w*ft(ind) + (1.0d0-w)*ft(ind+1)

    end if

end function lini

! Bi-linear interpolation
function bilini(t1,t2,ft,z1,z2) result(fz)


    real(8), intent(in) :: t1(:), t2(:), ft(:,:), z1, z2
    real(8) w, w1, w2, fz, x1, x2
    integer it, nt1, nt2, ind, ind1, ind2, i

    nt1 = size(t1)
    nt2 = size(t2)
        
    if (z1<=t1(1)) then
        x1 = t1(1)
    elseif (z1>=t1(nt1)) then
        x1 = t1(nt1)
    else ! t(1) < z < t(nt)
        x1 = z1
    end if
    
    ind1 = 0    
    do it = 1,nt1
        if (t1(it)>=x1) exit
        ind1 = ind1+1
    end do
    ind1 = max(1,ind1)
    
    w1 = (t1(ind1+1)-x1)/(t1(ind1+1)-t1(ind1))
    
    
    if (z2<=t2(1)) then
        x2 = t2(1)
    elseif (z2>=t2(nt2)) then
        x2 = t2(nt2)
    else ! t(1) < z < t(nt)
        x2 = z2
    end if
    
    ind2 = 0
    do it = 1,nt2
        if (t2(it)>=x2) exit
        ind2 = ind2+1
    end do    
    ind2 = max(1,ind2)
    
    w2 = (t2(ind2+1)-x2)/(t2(ind2+1)-t2(ind2))
    
    fz = w1*w2*ft(ind1,ind2) + w1*(1.0d0-w2)*ft(ind1,ind2+1) + (1.0d0-w1)*w2*ft(ind1+1,ind2) + (1.0d0-w1)*(1.0d0-w2)*ft(ind1+1,ind2+1)


end function bilini


subroutine ginicoef(dist,grid,ngrid,ns,gini,share,negW)

    ! Calculate Gini and Lorenz curve

    ! Inputs
    ! dist: measure
    ! grid: grid
    ! ngrid: number of grid points
    ! ns: sharenum

    ! Output
    ! gini: Gini Index
    ! share: Share in bin
    ! negW: fraction with non-positive
    implicit none

    integer, intent(in) :: ngrid, ns
    real(8), intent(in) :: dist(:), grid(:)
    real(8), intent(out) :: gini, negW
    real(8), intent(out), dimension(ns) :: share

    integer :: i1, j1
    real(8) :: aaa, bbb, ccc, ddd
    real(8), dimension(ngrid) :: xaxis, a1, yaxis, distn
    real(8), dimension(ngrid-1) :: dc, db
    real(8), dimension(ns) :: sharecum
    
    ddd = ns
    
    ! normalize dist if necessary
    distn = dist / sum(dist,1)           ! should be element-wise dividsion
    
    ! Cacluate the cdf
    xaxis=0.0d0
    do i1 = 1,ngrid
        xaxis(i1)=sum(distn(1:i1))
    end do

    ! Generate the Lorenz curve (yaxis)
    a1=0.0d0
    yaxis=0.0d0
    a1=grid*distn
    aaa=sum(a1)
    do i1 = 1,ngrid
        bbb=sum(a1(1:i1))
        yaxis(i1)=bbb/aaa
    end do
    ! Calcualte the gini coefficient of wealth
    do i1 = 1,ngrid-1
        dc(i1)=xaxis(i1+1)-xaxis(i1)
        db(i1)=(yaxis(i1+1)+yaxis(i1))*0.5d0
    end do

    bbb=sum(dc*db)
    aaa=0.5d0-bbb

    gini=2.0d0*aaa

    sharecum=0.0d0
    do j1 = 1,ns
        ccc = j1/ddd
        sharecum(j1) = lini(xaxis,yaxis,ccc)
    end do

    share = 0.0d0
    share(1) = sharecum(1)
    do j1 = 1,ns-2
        share(j1+1) = sharecum(j1+1) - sharecum(j1)
    end do
    share(ns) = 1 - sharecum(ns-1)

    ! share of nonpositive holders
    negW = 0.0d0
    do i1 = 1,ngrid
        if (grid(i1) <= 0) then
            negW = negW + distn(i1)
        end if
        if (grid(i1) > 0) then
            exit
        end if        
    end do
    
    

end subroutine ginicoef

subroutine compressdist(grid,ngrid,nn,distcomp,gridcomp)
    
        ! get an equally weighted original grid of size ngrid and compress into one
        ! of size nn and recalculate the new distribution
        implicit none

        integer, intent(in) :: ngrid, nn        ! original and target size of grid
        real(8), intent(in) :: grid(:)
        real(8), intent(out), dimension(nn) :: distcomp(nn), gridcomp(nn)
    
        integer :: pos, i1, i2
        real(8)  maxv, minv
        
        maxv = maxval(grid) + 0.001d0   ! to make sure loop always exits
        minv = minval(grid)
        
        gridcomp = linspace(minv,maxv,nn)
        distcomp = 0.0d0
        
        do i1=1,ngrid
            do i2=1,nn
                if(grid(i1) <= gridcomp(i2)) exit
            enddo
            distcomp(i2) = distcomp(i2) + 1
        enddo
        
        distcomp = distcomp / sum(distcomp)
        
    end subroutine compressdist

subroutine justgini(dist,grid,ngrid,gini)

    ! Calculate Gini and Lorenz curve
    ! just the gini and sorts and normalizes distribution and grid if necessary
    
    ! Inputs
    ! dist: measure
    ! grid: grid
    ! ngrid: number of grid points

    ! Output
    ! gini: Gini Index

    implicit none

    integer, intent(in) :: ngrid
    real(8), intent(in) :: dist(:), grid(:)
    real(8), intent(out) :: gini

    integer :: i1, j1
    real(8) :: aaa, bbb, ccc
    real(8), dimension(ngrid) :: xaxis, a1, yaxis, distn
    real(8), dimension(ngrid-1) :: dc, db
    
    
    ! normalize dist if necessary
    distn = dist / sum(dist,1)           ! should be element-wise dividsion
    
    ! Cacluate the cdf
    xaxis=0.0d0
    do i1 = 1,ngrid
        xaxis(i1)=sum(distn(1:i1))
    end do

    ! Generate the Lorenz curve (yaxis)
    a1=0.0d0
    yaxis=0.0d0
    a1=grid*distn
    aaa=sum(a1)
    do i1 = 1,ngrid
        bbb=sum(a1(1:i1))
        yaxis(i1)=bbb/aaa
    end do
    ! Calcualte the gini coefficient
    do i1 = 1,ngrid-1
        dc(i1)=xaxis(i1+1)-xaxis(i1)
        db(i1)=(yaxis(i1+1)+yaxis(i1))*0.5d0
    end do

    bbb=sum(dc*db)
    aaa=0.5d0-bbb 
    
    gini=2.0d0*aaa

end subroutine justgini


subroutine markovss(pij,sizex,tol,xss)

    integer, intent(in) :: sizex
    real(8), intent(in), dimension(sizex,sizex) :: pij
    real(8), intent(in) :: tol
    real(8), intent(out), dimension(sizex) :: xss

    integer :: ix, jx
    real(8) :: dist
    real(8), dimension(sizex) :: xold, xnew

    xold = (1.0d0/dble(sizex))
    xnew = xold
    dist = 1000.0d0

    do while (dist > tol*2)
        do jx = 1,sizex
            xnew(jx) = 0.0d0
            do ix = 1,sizex
                xnew(jx) = xnew(jx) + pij(ix,jx)*xold(ix)
            end do
        end do
        dist = maxval(abs(log(xnew)-log(xold)))
        xold = xnew
    end do
    xss = xnew

end subroutine markovss



    
Subroutine gridlookup(r,grid,valplace,iloc)
! Looking for index of grid about a value
! Code from Aubhik's matlab code
! R(input) :  number of inner grids 
! grid(input) : grid we have
! valplace(input) : value we have
! iloc(output) : index number
    Implicit none
    integer, intent(in) :: r
    Real(8), intent(in) :: valplace, grid(r+2)
    integer, intent(out) :: iloc
    integer :: ilow, ihigh, distance, inow
    Real(8) :: valnow

    ilow = 1
    ihigh = r+2
    distance = 2

Do while (distance > 1)

    inow = floor((ilow + ihigh)/2.);
    valnow = grid(inow);


    if (valnow > valplace) then
        ihigh = inow;
    else
        ilow = inow;
    end if

    distance = ihigh - ilow;

end do

iloc = ilow;

end subroutine gridlookup 

subroutine bsort(N,X,Y,Ind)
implicit none
integer, intent(in) :: N
real(8), intent(in) :: X(N)
real(8), intent(out) :: Y(N)
integer, intent(out) :: Ind(N)

real(8) ::temp,XX(N)
integer :: I, J, K,tempind(N),tempi

Logical :: sorted

sorted = .false.
tempind=0

do J=1,N
    tempind(J)=J
end do

K=0
XX= X
do while(.not.Sorted)
    Sorted = .true.
    K = K + 1
    do J = 1, N-K
        if (XX(J)>XX(J+1)) then
            temp = XX(J)
            XX(J)=XX(J+1)
            XX(J+1)=temp
            tempi=tempind(J)
            tempind(J)=tempind(J+1)
            tempind(J+1) = tempi
            sorted = .false.
        end if
    end do
end do
Y= XX
Ind=tempind

end subroutine bsort


subroutine sort(x, size)
    ! This subroutine receives an array x() and sorts it into ascending
    implicit none
    integer, intent(in) :: size
    real(8), intent(inout) :: x(size)
    integer :: i, location
      do i=1, size-1			! except for the last
         location = findminimum(x, i, size)	! find min from this to last
         call swap(x(i), x(location))	! swap this and the minimum
      end do
end subroutine  sort

! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a)
  implicit none
  real(8) :: a(:)
  real(8) x, t
  integer :: first = 1, last
  integer i, j

  last = size(a, 1)
  x = a( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) call quicksort(a(first : i - 1))
  if (j + 1 < last)  call quicksort(a(j + 1 : last))
end subroutine quicksort

subroutine  swap(a, b)
   !This subroutine swaps the values of its two formal arguments.
    implicit none
    real(8), intent(inout) :: a, b
    integer :: temp

      temp = a
      a    = b
      b    = temp
end subroutine swap

subroutine quart(dist,grid,ngrid,qq)
    ! calculates the 25, 50, and 75th prctl of dist
    
    ! Inputs
    ! dist: weights
    ! grid: values (sorted)?
    ! ngrid: number of grid points

    ! Output
    ! qq: vector containing percentiles
    
    implicit none
    integer, intent(in) :: ngrid
    real(8), intent(in) :: dist(:), grid(:)
    real(8), intent(out) :: qq(3)

    integer :: i1, i2
    real(8) :: total, q1, q2, q3
    real(8) :: cumdist(ngrid)

    total = sum(dist)
    q1 = total*0.25d0
    q2 = total*0.50d0
    q3 = total*0.75d0
    
    cumdist=0.0d0
    do i1 = 1,ngrid
        cumdist(i1)=sum(dist(1:i1))
    end do
    
    do i2=1,ngrid
        if (cumdist(i2) > q1) exit
    end do
    qq(1) = grid(i2)
    
    do i2=1,ngrid
        if (cumdist(i2) > q2) exit
    end do
    qq(2) = grid(i2)
    
    do i2=1,ngrid
        if (cumdist(i2) > q3) exit
    end do
    qq(3) = grid(i2)    

    
end subroutine quart



function  findminimum(x, s, e)
    ! This function returns the location of the minimum in the section
    ! between s and e.
    implicit none
    real(8), intent(in) :: x(:)
    integer, intent(in) :: s, e
    integer :: findminimum
    integer :: minimum, location, i


      minimum  = x(s)		! assume the first is the min
      location = s			! record its position
      do i = s+1, e		! start with next elements
         if (x(i) < minimum) then	!   if x(i) less than the min
            minimum  = x(i)		!      Yes, a new minimum found
            location = i                !      record its position
         end if 
      end do
      findminimum = location        	! return the position
end function  findminimum

    function loggrid(xlow,xhigh,xnum,deggrid) result(xgrid)
        ! custormized non-evenly spaced grud

        integer, intent(in) :: xnum
        real(8), intent(in) :: xlow, xhigh, deggrid
        real(8), dimension(xnum) :: xgrid
    
        integer i
    
        if (deggrid == 0.0d0) then
            xgrid = logspace(log(xlow - 1.0d0*xlow + 1.0d0)/log(10.0d0), log(xhigh - 1.0d0*xlow + 1.0d0)/log(10.0d0), xnum)
            xgrid = xgrid + xlow - 1.0d0
        else
            xgrid = logspace(0.0d0,deggrid,xnum) - 10.0d0**(0.0d0)
            xgrid = xgrid/(10.0d0**deggrid - 10.0d0**0.0d0)
            xgrid = xgrid*(xhigh-xlow) + xlow 
        end if
        
    end function loggrid
    
    
    function invmat(A) result(Ainv)
    
        implicit none 
        
        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(size(A,1),size(A,2)) :: Ainv
    
        real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
    
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
    
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
    
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
    
        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if
    
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
    
    end function invmat
    
    
    function meansim(simdata) result(mean)
    
        implicit none
        real(8) :: mean, num
        real(8), intent(in), dimension(:) :: simdata
    
        num = dble(size(simdata))
        mean = sum(simdata)/num    
    
        end function meansim
    
    
    function varsim(simdata) result(var)
    
        implicit none
        real(8) :: mean, var, num
        real(8), intent(in), dimension(:) :: simdata
    
        num = dble(size(simdata))
        mean = sum(simdata)/num
        var = sum((simdata-mean)**2.0d0)/num
    
    end function varsim
    
    
    subroutine findiw(x,xx,ix,wx)

        real(8), intent(in) :: x, xx(:)
        integer, intent(out) :: ix
        real(8), intent(out) :: wx
        integer nx


        nx = size(xx,1)

        if (x<=xx(1)) then
            ix = 1
        elseif (x>=xx(nx)) then
            ix = nx-1
        else
            do ix=1,nx-1

                if (x<xx(ix+1)) exit

            end do
        
        end if

        wx = (xx(ix+1)-x)/(xx(ix+1)-xx(ix))

    end subroutine findiw
    
        
    function dimreduce3_new(n1,n2,n3) result(indmat)
           
        integer, intent(in) :: n1, n2, n3        
        real(8), allocatable :: indmat(:,:)

        integer :: i1, i2, i3, size, index        
        
        size = n1*n2*n3
        allocate(indmat(3,size))

        do i1 = 1,n1
            do i2 = 1,n2
                do i3 = 1,n3
                    index = (i1-1)*n3*n2 + (i2-1)*n3 + i3
                    indmat(1,index) = i1
                    indmat(2,index) = i2                
                    indmat(3,index) = i3                
                end do                
            end do
        end do
    

    end function dimreduce3_new
    
    
    subroutine getindex(cumdist,rdraw,ngrid,index)
    
    ! subroutine finds the position where the random draw is located on grid with cumulative distribution cumdist
    integer, intent(in) :: ngrid
    real(8), intent(in) :: rdraw
    real(8), intent(in) :: cumdist(:)
    integer, intent(out) :: index
    
    integer :: i1, pos
                        
    index = 1
    do  i1=1,ngrid
        if (rdraw < cumdist(i1)) exit
        index = index + 1
    end do
    
    if (index >ngrid) then
        write(*, '(a)') ' Error: Index is larger than grid size! '
    end if
    
    
    end subroutine getindex
    
    
    
function calc_r(x, y) result(r_out)
    implicit none
    real(8), intent(in) :: x(:), y(:)
    real(8)    :: r_out
    integer :: size_x, size_y, i
    real(8)    :: mean_x, mean_y, cov, var_x, var_y

    size_x = size(x)
    size_y = size(y)
    if (size_x == 0 .or. size_y == 0) then
      print *, "[ERROR] array size == 0"
      stop
    end if
    if (size_x /= size_y) then
      print *, "[ERROR] size(X) != size(Y)"
      stop
    end if

    mean_x = sum(x) / size_x
    mean_y = sum(y) / size_y
    cov = sum((x(1:size_x) - mean_x) * (y(1:size_x) - mean_y))
    var_x = sum((x(1:size_x) - mean_x) * (x(1:size_x) - mean_x))
    var_y = sum((y(1:size_x) - mean_y) * (y(1:size_x) - mean_y))
    r_out = (cov / sqrt(var_x)) / sqrt(var_y)
  end function calc_r  

        
    
end module mod_utils
    