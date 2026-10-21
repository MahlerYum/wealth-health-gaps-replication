module sort_all

contains
SUBROUTINE bsort(N,X,Y,Ind)
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



SUBROUTINE bsort2(ini,lst,X,Y,Ind)
implicit none
integer, intent(in) :: ini,lst
real(8), intent(in) :: X(ini:lst)
real(8), intent(out) :: Y(1:lst-ini+1)
integer, intent(out) :: Ind(1:lst-ini+1)
real(8) ::temp,XX(1:lst-ini+1)
integer :: I, J, K,tempind(1:lst-ini+1),tempi
Logical :: sorted

sorted = .false.

tempind=0
do J=1,lst-ini+1
    tempind(J)=J
end do
K=0
XX(1:lst-ini+1)= X(ini:lst)
do while(.not.Sorted)
    Sorted = .true.
    K = K + 1
    do J = 1, lst-ini+1-K
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
Ind=tempind+ini-1
end subroutine bsort2




end module sort_all

                