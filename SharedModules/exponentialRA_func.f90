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
  
  public :: exponential, expTau, expG, expH, expG2, expF2

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
        d2d = 0.2263358514260129, d3d = 0.04721469893686252, d4d = 0.006883236664917246, &
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




  ! Coefficients for numerator in rational approximation
        real(defFlt), parameter :: g1n = -0.08335775885589858, g2n = -0.003603942303847604, &
        g3n = 0.0037673183263550827, g4n = 0.00001124183494990467, g5n = 0.00016837426505799449

        ! Coefficients for denominator in rational approximation
        real(defFlt), parameter :: g1d = 0.7454048371823628, g2d = 0.23794300531408347, &
        g3d = 0.05367250964303789, g4d = 0.006125197988351906, g5d = 0.0010102514456857377

  ! Coefficients for numerator
        real(defFlt), parameter :: f1n5 = 0.000136757570702, &
        f1n4 = 0.000640849175509, f1n3 = 0.007675127136944, &
        f1n2 = 0.035904163235632, f1n1 = 0.166666147003676

  ! Coefficients for denominator
          real(defFlt), parameter :: f1d0 = 1.000000000000000, &
          f1d1 = 0.715333312893290, f1d2 = 0.254155566312370, &
          f1d3 = 0.056133925714270, f1d4 = 0.009476002327853, &
          f1d5 = 0.000914563747782, f1d6 = 0.000136757570702

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

! Computes y = 1/x-(1-exp(-x))/x**2 using a 5/6th order rational approximation.
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

! Computes H : y = (1-exp(-x)*(1+x))/x**2 
! using a 5/7th order rational approximation. 
! FROM: OpenMoC https://github.com/mit-crpg/OpenMOC/blob/7c8c9460c1c95f68dae102a402a39afa233a0b8c/src/exponentials.h#L9

elemental function expH(tau) result(x)
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
end function expH

! Computes G2 : y = 2/3 - (1 + 2/x) * (1/x + 0.5 - (1 + 1/x) * (1-exp(-x)) / x) 
! using a 5/5th order rational approximation,
! FROM: OpenMoC https://github.com/mit-crpg/OpenMOC/blob/7c8c9460c1c95f68dae102a402a39afa233a0b8c/src/exponentials.h#L9

elemental function expG2(tau) result(x)
  real(defFlt), intent(in)    :: tau
  real(defFlt)                :: x
  real(defFlt)                :: den, num

  x = tau 
  ! Calculate numerator
  num = g5n*x + g4n
  num = num*x + g3n
  num = num*x + g2n
  num = num*x + g1n
  num = num*x

  ! Calculate denominator
  den = g5d*x + g4d
  den = den*x + g3d
  den = den*x + g2d
  den = den*x + g1d
  den = den*x + 1.0_defFlt

  ! Compute final value
  x = num / den

end function expG2

! Computes F2 : y = (x-2+exp(-x)*(2+x))/x**2
! using a 5/6th order rational approximation,
! FROM: OpenMoC https://github.com/mit-crpg/OpenMOC/blob/7c8c9460c1c95f68dae102a402a39afa233a0b8c/src/exponentials.h#L9

!!!!not current correct form we use.

elemental function expF2(tau) result(x)
  real(defFlt), intent(in)    :: tau
  real(defFlt)                :: x
  real(defFlt)                :: den, num

  x = tau

  ! Compute numerator
  num = f1n5*x + f1n4
  num = num*x + f1n3
  num = num*x + f1n2
  num = num*x + f1n1
  num = num*x

  ! Compute denominator
  den = f1d6*x + f1d5
  den = den*x + f1d4
  den = den*x + f1d3
  den = den*x + f1d2
  den = den*x + f1d1
  den = den*x + f1d0

  x = num / den
end function expF2



    
end module exponentialRA_func
