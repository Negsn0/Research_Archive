! �萔���g�����߂�module �{���x�Œ�`�ρB
! PI,

module constants
    implicit none
    
    ! (�{���x)
    integer,parameter :: dp = selected_real_kind(15, 307)
    
    !�~����pi
    real(kind=dp),parameter :: PI = 4.0_dp * atan(1.0_dp)
end module constants

!���w��(����)�֐����`
!�K��A��`���A���􉽊֐��A���Ώ̃|�e���V�����̊��
module math_unit
    implicit none
    contains
    
    !�K��
    function factorial(N) result(retval)
        integer, intent(in) :: N
        integer :: i
        real:: retval
        
        retval = 1.0
        if ( N /= 0 ) then
            do i = 1, N
                retval = retval * real(i)
            end do
        end if
    end function factorial

    ! ��`��
    function daikei(f,a,b,n) result(result)
        interface
            function f(x)
                real, intent(in) :: x
                real :: f
            end function f
        end interface

        real, intent(in) :: a, b
        real :: result, x, h
        integer :: i, n
        h=(b - a) / n
        result = 0.0
        x = a
        do i = 0, n
            result = result + f(x)
            x = x + h
        end do
        result = result * h
    end function daikei

    !���􉽊֐�
    function hyper(n,gamma,x) result(Mn)
        
        integer, intent(in) :: n
        real, intent(in) :: gamma,x
        real :: Mn, Mn2 = 1.0 , Mn1
        integer :: i
        
        if ( n==0 ) then
            Mn = 1.0
        else if ( n==1 ) then
            Mn = 1.0 - x / gamma
        else
            Mn = 0.0
            Mn1 = 1.0 - x / gamma
            do i = 2, n
                Mn = ((2.0 * real(i) + gamma - x - 2.0) * Mn1 - (real(i) - 1.0) * Mn2 ) / (real(i) + gamma - 1.0)
                Mn2 = Mn1
                Mn1 = Mn
            end do
        end if
    end function hyper

    !���
    function radial(n0,l0,r,a) result(retval)
        integer, intent(in) :: n0, l0
        real, intent(in) :: r, a
        real :: retval, nn, q, n, l
        n = real(n0)
        l = real(l0)
        nn = sqrt((2 * gamma(l + 1.5 + n))/(factorial(n0) * gamma(l+1.5) * gamma(l+1.5)))
        q = a * r
        retval = sqrt(a) * N * q ** (l+1.0) * exp(-0.5*q*q) * hyper(n0,l+1.5,q * q)
    end function radial

end module math_unit

! �����w�I�Ȋ֐��⎮���`�B
module physics_unit
    implicit none
    contains
    ! single - particle kinetic energy 

    
    !Wood - Saxon
    function WS(r,R0,a) result(retval)
        real, intent(in) :: r,R0,a
        real :: retval
        retval = 1.0 / (1.0 + exp((r - R0) / a))
    end function WS

    !WS ����
    function dWS(r,R0,a) result(retval)
        real, intent(in) :: r,R0,a
        real :: retval,X
        X = exp((r - R0) / a)
        retval = - 1.0 * (1.0 / (1.0 + X) * (1.0 + X)) * X / a
    end function dWS


end module physics_unit