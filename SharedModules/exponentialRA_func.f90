module exponentialRA_func
  !! This module contains a function for efficiently computing an exponential
  !! for use in MoC implementations. RA stands for rational approximation.
  !! This is based on the implementation given in Minray:
  !! github.com/jtramm/minray/blob/master/cpu_srce/flux_attenuation_kernel.c
  !! I believe this originates from the M&C 2019 publication:
  !! "Adding a third level of parallelism to OpenMOC"
  !!

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private
  
  public :: exponential, expTau, expG

  ! Numerator coefficients in rational approximation for 1 - exp(-tau)
  real(defFlt), parameter :: c1n = -1.0000013559236386308, c2n = 0.23151368626911062025,&
          c3n = -0.061481916409314966140, c4n = 0.0098619906458127653020, c5n = -0.0012629460503540849940, &
          c6n = 0.00010360973791574984608, c7n = -0.000013276571933735820960

  ! Denominator coefficients in rational approximation for 1 - exp(-tau)
  real(defFlt), parameter :: c0d = ONE, c1d = -0.73151337729389001396, c2d = 0.26058381273536471371, &
          c3d = -0.059892419041316836940, c4d = 0.0099070188241094279067, c5d = -0.0012623388962473160860, &
          c6d = 0.00010361277635498731388, c7d = -0.000013276569500666698498


  ! Numerator coefficients in rational approximation for 1/x-(1-exp(-x))/x**2
        real(defFlt), parameter :: d0n = 0.5, d1n = 0.176558112351595,&
        d2n = 0.04041584305811143, d3n = 0.006178333902037397, d4n = 0.0006429894635552992 , d5n = 0.00006064409107557148

  ! Denominator coefficients in rational approximation for 1/x-(1-exp(-x))/x**2
        real(defFlt), parameter :: d0d = 1.0, d1d = 0.6864462055546078,&
        d2d = 0.02263358514260129, d3d = 0.04721469893686252, d4d = 0.006883236664917246, &
        d5d = 0.0007036272419147752 , d6d = 0.00006064409107557148 
      
  ! Coefficients for numerator in rational approximation
        real(defFlt), parameter :: h0n = 0.5, h1n = 0.05599412483229184, &
        h2n = 0.01294939509305754, h3n = 0.002341166644220405, &
        h4n = 0.00003686858969421769, h5n = 0.00004220477028150503

        ! Coefficients for denominator in rational approximation
        real(defFlt), parameter :: h0d = 1.0, h1d = 0.7787274561075199, &
        h2d = 0.2945145030273455, h3d = 0.07440380752801196, &
        h4d = 0.01220791761275212, h5d = 0.002354181374425252, &
        h6d = 0.00003679462493221416, h7d = 0.00004220477028150503

contains

  !!
  !! Computes x = 1 - exp(-tau) for use in MoC calcs
  !! Tau is the optical distance
  !!
  elemental function exponential(tau) result(x)
    real(defFlt), intent(in)    :: tau
    real(defFlt)                :: x
    real(defFlt)                :: den, num

    x = -tau
    den = c7d
    den = den * x + c6d
    den = den * x + c5d
    den = den * x + c4d
    den = den * x + c3d
    den = den * x + c2d
    den = den * x + c1d
    den = den * x + c0d

    num = c7n
    num = num * x + c6n
    num = num * x + c5n
    num = num * x + c4n
    num = num * x + c3n
    num = num * x + c2n
    num = num * x + c1n
    num = num * x

    x = num / den

  end function exponential

!!
!! Computes x = (1 - exp(-tau))/tau for use in MoC calcs
!! Tau is the optical distance
!!

elemental function expTau(tau) result(x)
real(defFlt), intent(in)    :: tau
real(defFlt)                :: x
real(defFlt)                :: den, num

  x = -tau
  den = c7d
  den = den * x + c6d
  den = den * x + c5d
  den = den * x + c4d
  den = den * x + c3d
  den = den * x + c2d
  den = den * x + c1d
  den = den * x + c0d

  num = c7n
  num = num * x + c6n
  num = num * x + c5n
  num = num * x + c4n
  num = num * x + c3n
  num = num * x + c2n
  num = num * x + c1n
  num = -num 

  x = num / den

end function expTau

!This method computes 1/x-(1-exp(-x))/x**2 using a 5/6th order rational approximation.
!FROM: OpenMoC https://github.com/mit-crpg/OpenMOC/blob/7c8c9460c1c95f68dae102a402a39afa233a0b8c/src/exponentials.h#L9

elemental function expG(tau) result(x)
  real(defFlt), intent(in)    :: tau
  real(defFlt)                :: x
  real(defFlt)                :: den, num

  x = tau

  den = d6d * x + d5d
  den = den * x + d4d
  den = den * x + d3d
  den = den * x + d2d
  den = den * x + d1d
  den = den * x + d0d
  !den = 1/den

  num = d5n * x + d4n
  num = num * x + d3n
  num = num * x + d2n
  num = num * x + d1n
  num = num * x + d0n

  x = num / den

end function expG

! Computes H exponential term using a rational approximation (1-exp(-x)*(1+x))/x**2 
! using a 5/7th order rational approximation. 
! FROM: OpenMoC https://github.com/mit-crpg/OpenMOC/blob/7c8c9460c1c95f68dae102a402a39afa233a0b8c/src/exponentials.h#L9

elemental subroutine expH_fractional(x)
  real(defFlt), intent(in)    :: tau
  real(defFlt)                :: x
  real(defFlt)                :: den, num

  x = tau

  num = h5n*x + h4n
  num = num*x + h3n
  num = num*x + h2n
  num = num*x + h1n
  num = num*x + h0n

  den = h7d*x + h6d
  den = den*x + h5d
  den = den*x + h4d
  den = den*x + h3d
  den = den*x + h2d
  den = den*x + h1d
  den = den*x + h0d

  x = num / den
end subroutine expH_fractional



    
end module exponentialRA_func
