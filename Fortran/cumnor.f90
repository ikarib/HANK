subroutine cumnor ( arg, cum, ccum )

!*****************************************************************************80
!
!! CUMNOR computes the cumulative normal distribution.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine dependent
!    constants.
!
!  Author: 
!
!    William Cody
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    Argonne, IL 60439
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation, 
!    1969, pages 631-637.
!
!    William Cody, 
!    Algorithm 715: 
!    SPECFUN - A Portable FORTRAN Package of Special Function Routines 
!    and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, 1993, pages 22-32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the upper limit of integration.
!
!    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
!    complementary CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) EPS, the argument below which anorm(x) 
!    may be represented by 0.5 and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0D+00 + X = 1.0D+00   to machine precision.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, &
    1.6102823106855587881D+02, &
    1.0676894854603709582D+03, &
    1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, &
    9.7609855173777669322D+02, &
    1.0260932208618978205D+04, &
    4.5507789335026729956D+04 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, &
    8.8831497943883759412D+00, &
    9.3506656132177855979D+01, &
    5.9727027639480026226D+02, &
    2.4945375852903726711D+03, &
    6.8481904505362823326D+03, &
    1.1602651437647350124D+04, &
    9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, &
    2.3538790178262499861D+02, &
    1.5193775994075548050D+03, &
    6.4855582982667607550D+03, &
    1.8615571640885098091D+04, &
    3.4900952721145977266D+04, &
    3.8912003286093271411D+04, &
    1.9685429676859990727D+04 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  integer i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, &
    1.274011611602473639D-01, &
    2.2235277870649807D-02, &
    1.421619193227893466D-03, &
    2.9112874951168792D-05, &
    2.307344176494017303D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, &
    4.68238212480865118D-01, &
    6.59881378689285515D-02, &
    3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) xsq
!
!  Machine dependent constants
!
  eps = epsilon ( 1.0D+00 ) * 0.5D+00

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    if ( eps < y ) then
      xsq = x * x
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    cum = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = cum
    cum = 0.5D+00 + temp
    ccum = 0.5D+00 - temp
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
  else if ( y <= root32 ) then

    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    cum = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call r8_swap ( cum, ccum )
    end if
!
!  Evaluate ANORM for sqrt(32) < |X|.
!
  else

    cum = 0.0D+00
    xsq = 1.0D+00 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do

    cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    cum = ( sqrpi - cum ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) &
      * exp ( - del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call r8_swap ( cum, ccum )
    end if

  end if

  if ( cum < tiny ( cum ) ) then
    cum = 0.0D+00
  end if

  if ( ccum < tiny ( ccum ) ) then
    ccum = 0.0D+00
  end if

  return
end

subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8 values.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
