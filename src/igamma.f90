!=========================================================
! subroutine: digamiv()
!   vector version of digami()
!   subroutine of digami copied from r. j. moore algorithm
!   as 187 applied statistics (1982) vol 31, no 3
!   d(1): 1st derivative of i(x,p) w.r.t. x
!   d(2): 2nd derivative of i(x,p) w.r.t. x**2
!   d(4): 2nd derivative of i(x,p) w.r.t. p**2
!   d(3): 1st derivative of i(x,p) w.r.t. p
!   d(5): 2nd derivative of i(x,p) w.r.t. x and p
!   d(6): value of i(x,p)
!=========================================================
subroutine digamiv(num,d,x,p,plimit,ifault)
!dec$ attributes dllexport::digamiv
!dec$ attributes c, reference, alias:'digamiv_'::digamiv
  implicit real(kind=8)(a-h,o-z)
  integer i,j,index0,num
  dimension d(6*num),x(num),p(num),dd(6)
  index0=0
  do i=1,num
    xx=x(i)
	pp=p(i)
    call digami(dd,xx,pp,plimit,ifault)
    do j=1,6
      d(index0+j)=dd(j)
    end do
    index0=index0+6
  end do
  return
end subroutine

!=========================================================
! subroutine: digami()
!   subroutine of digami copied from r. j. moore algorithm
!   as 187 applied statistics (1982) vol 31, no 3
!   d(1): 1st derivative of i(x,p) w.r.t. x
!   d(2): 2nd derivative of i(x,p) w.r.t. x**2
!   d(3): 1st derivative of i(x,p) w.r.t. p
!   d(4): 2nd derivative of i(x,p) w.r.t. p**2
!   d(5): 2nd derivative of i(x,p) w.r.t. x and p
!   d(6): value of i(x,p)
!=========================================================
subroutine digami(d,x,p,plimit,ifault)
!!!dec$ attributes dllexport::digami
!!!dec$ attributes c, reference, alias:'digami_'::digami
  implicit real(kind=8)(a-h,o-z)
  dimension pn(6),d(6),dp(6),dpp(6)
  data e,oflo,tmax,zero/1.0d-12,1.0d30,100.0d0,1.0d-30/

  d(1:6)=0.0d0
  ! Check that we have valid value for x and p
  if (p <= 0.0d0 .or. x < 0.0d0) then
    ifault=1
    return
  end if
  if (x == 0.0d0) then
    return
  end if

  ! Use a normal approximation if p > plimit
  if(p > plimit) then
    return
  end if

  ! gplog: the value of log(gamma(p))
  ! gp1log: the log(gamma(p+1)) = log(p) + log(gamma(p))
  ! psip: digamma(p) ----- 1st derivative of log(gamma(p))
  ! psip1: digamma(p+1) = 1/p + digamma(p)
  ! psidp: trigamma(p) --- derivative of digamma(p)
  ! psidp1: trigamma(p+1) = trigamma(p) - 1/(p*p)

  ifault=0
  ifail=1
  gplog=alngam(p)
  gp1log=dlog(p)+gplog
  psip=digama(p,ifault)
  psip1=(1.0d0/p)+psip
  psidp=trigam(p,ifault)
  psidp1=psidp-1.0d0/(p**2.0d0)

  ! derivative with respect to x
  pm1=p-1.0d0
  xlog=dlog(x) ! must restrict to x
  d(1)=dexp(-gplog+pm1*xlog-x)
  d(2)=d(1)*(pm1/x-1.0d0)
  d(5)=d(1)*(xlog-psip)

  ! derivatives with respect to p
  if (x > 1.0d0 .and. x >= p) goto 30

  ! series expansion
  f=dexp(p*xlog-gp1log-x)
  dfp=f*(xlog-psip1)
  dfpp=dfp*dfp/f-f*psidp1

  tmaxp=tmax+p
  c=1.0d0
  s=1.0d0
  cp=0.0d0
  cpp=0.0d0
  dsp=0.0d0
  dspp=0.0d0
  a=p
1 a=a+1.0d0
  cpc=cp/c
  cp=cpc-1.0d0/a
  cpp=cpp/c-cpc*cpc+1.0d0/(a*a)
  c=c*x/a
  cp=cp*c
  cpp=cpp*c+cp*cp/c
  s=s+c
  dsp=dsp+cp
  dspp=dspp+cpp
  if (a > tmaxp) goto 1001
  if (c > e*s) goto 1
  d(6)=s*f
  d(3)=s*dfp+f*dsp
  d(4)=s*dfpp+2.0d0*dfp*dsp+dspp*f
  return

  ! continued fraction expansion
  30 f=dexp(p*xlog-gplog-x)
  dfp=f*(xlog-psip)
  dfpp=dfp*dfp/f-f*psidp

  a=pm1
  b=x+1.0d0-a
  term=0.0d0
  pn(1)=1.0d0
  pn(2)=x
  pn(3)=x+1.0d0
  pn(4)=x*b
  s0=pn(3)/pn(4)
  do 31 i = 1, 4
    dp(i) = 0.0d0
    dpp(i) = 0.0d0
  31 continue
  dp(4) = -x

32 a=a-1.0d0
  b=b+2.0d0
  term=term+1.0d0
  an=a*term
  pn(5)=b*pn(3)+an*pn(1)
  pn(6)=b*pn(4)+an*pn(2)
  dp(5)=b*dp(3)-pn(3)+an*dp(1)+pn(1)*term
  dp(6)=b*dp(4)-pn(4)+an*dp(2)+pn(2)*term
  dpp(5)=b*dpp(3)+an*dpp(1)+2.0d0*(term*dp(1)-dp(3))
  dpp(6)=b*dpp(4)+an*dpp(2)+2.0d0*(term*dp(2)-dp(4))

  if (dabs(pn(6)) < zero) goto 35
  s=pn(5)/pn(6)
  c=dabs(s-s0)
  if (c*p > e) goto 34
  if (c <= e*s) goto 42

  34 s0 = s
  35 do 36 i=1,4
    i2=i+2
    dp(i)=dp(i2)
    dpp(i)=dpp(i2)
    pn(i)=pn(i2)
  36 continue

  if (term > tmax) goto 1001
  if (dabs(pn(5)) < oflo) goto 32
  do 41 i=1,4
    dp(i)=dp(i)/oflo
    dpp(i)=dpp(i)/oflo
    pn(i)=pn(i)/oflo
  41 continue
  goto 32

  42 d(6)=1.0d0-f*s
  dsp=(dp(5)-s*dp(6))/pn(6)
  dspp=(dpp(5)-s*dpp(6)-2.0d0*dsp*dp(6))/pn(6)
  d(3)=-f*dsp-s*dfp
  d(4)=-f*dspp-2.0d0*dsp*dfp-s*dfpp
  return

  ! set fault indicator
  1001 ifault=1
  return
end subroutine

!===================
! function: digama()
!===================
function digama(ww,ifault)
  implicit real(kind=8)(a-h,o-z)
  integer::ifault
  data ss,c,s3,s4,s5,d1/1.0d-5,8.50d0,8.3333333d-2,8.333333d-3,3.968253968d-3,-0.5772156649d0/

  ! Check argument is positive
  digama=0.0d0
  yy=ww
  ifault=1
  if (yy <= 0.0d0) return
  ifault = 0

  ! Use approximation if argument <= ss
  if (yy > ss) goto 600
  digama=d1-1.0d0/yy
  return

  ! Reduce to digamma(ww+n), (ww+n) >= c
600 continue
  if (yy >= c) goto 601
  digama=digama-1.0d0/yy
  yy=yy+1.0d0
  goto 600

  ! Use stirling if argument >= c
601 continue
  r=1.0d0/yy
  digama=digama+dlog(yy)-0.5d0 * r
  r=r*r
  digama=digama-r*(s3-r*(s4-r*s5))
  return
end function

!===================
! function: trigam()
!===================
function trigam(ww,ifault)
  implicit real(kind=8)(a-h,o-z)
  integer::ifault
  data a,b,one,half/1.0d-4,5.0d0,1.0d0,0.50d0/
  ! b2,b4,b6 and b6 are Bernoulli numbers
  data b2,b4,b6,b8/0.16666667d0,-0.03333333d0,0.02380952381d0, -0.03333333d0/

  ! Check for positive value of ww
  trigam=0.0d0
  ifault=1
  if (ww <= 0.0d0) return
  ifault=0
  z=ww

  ! Use small value approximation if x <= a
  if (z > a) goto 700
  trigam=1.0d0/(z*z)
  return

  ! Increase argument to (x + 1) >= b
700 continue
  if (z >= b) goto 701
  trigam=trigam+1.0d0/(z*z)
  z=z+1.0d0
  goto 700

  ! Apply asymptotic formula if argument >= b
701 continue
  yy=1.0d0/(z*z)
  trigam=trigam+half*yy+(one+yy*(b2+yy*(b4+yy*(b6+yy*b8))))/z
  return
end function

!===================
! function: alngam()
!===================
function alngam(xvalue) result(fn_val)

!     ALGORITHM AS245  APPL. STATIST. (1989) VOL. 38, NO. 2

!     Calculation of the logarithm of the gamma function

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

implicit none
real(kind=8), intent(in) :: xvalue
real(kind=8)             :: fn_val

! Local variables
real(kind=8) :: x, x1, x2, y

!     Coefficients of rational functions

real(kind=8), parameter :: r1(9) = (/ -2.66685511495d0, -2.44387534237d1,  &
                                      -2.19698958928d1,  1.11667541262d1,  &
                                       3.13060547623d0,  6.07771387771d-1, &
                                       1.19400905721d1,  3.14690115749d1,  &
                                       1.52346874070d1 /)
real(kind=8), parameter :: r2(9) = (/ -7.83359299449d1, -1.42046296688d2,  &
                                       1.37519416416d2,  7.86994924154d1,  &
                                       4.16438922228d0,  4.70668766060d1,  &
                                       3.13399215894d2,  2.63505074721d2,  &
                                       4.33400022514d1 /)
real(kind=8), parameter :: r3(9) = (/ -2.12159572323d5,  2.30661510616d5,  &
                                       2.74647644705d4, -4.02621119975d4,  &
                                      -2.29660729780d3, -1.16328495004d5,  &
                                      -1.46025937511d5, -2.42357409629d4,  &
                                      -5.70691009324d2 /)
real(kind=8), parameter :: r4(5) = (/ 2.79195317918525d1, 4.917317610505968d-1, &
                                      6.92910599291889d-2, 3.350343815022304d0, &
                                      6.012459259764103d0 /)

!     Fixed constants

real(kind=8), parameter :: alr2pi = 9.18938533204673d-1, four = 4.0d0,  &
                           half = 5.0d-1, one = 1.0d0, onep5 = 1.5d0,    &
                           twelve = 1.2d1, zero = 0.0d0

!     Machine-dependant constants.
!     A table of values is given at the top of page 399 of the paper.
!     These values are for the IEEE double-precision format for which
!     B = 2, t = 53 and U = 1023 in the notation of the paper.

real(kind=8), parameter :: xlge = 5.10d6, xlgst = huge(1.0d0)

x = xvalue
fn_val = zero

!     Test for valid function argument

! AS 245: Argument x too large
if (x >= xlgst) return

! AS 245: Argument x <= 0
if (x <= zero) return

!     Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined

if (x < onep5) then
  if (x < half) then
    fn_val = -dlog(x)
    y = x + one

!     Test whether X < machine epsilon

    if (y == one) return
  else
    fn_val = zero
    y = x
    x = (x - half) - half
  end if
  fn_val = fn_val + x * ((((r1(5)*y + r1(4))*y + r1(3))*y + r1(2))*y + r1(1)) / &
                    ((((y + r1(9))*y + r1(8))*y+ r1(7))*y + r1(6))
  return
end if

!     Calculation for 1.5 <= X < 4.0

if (x < four) then
  y = (x - one) - one
  fn_val = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x + r2(1)) /  &
               ((((x + r2(9))*x + r2(8))*x + r2(7))*x+ r2(6))
  return
end if

!     Calculation for 4.0 <= X < 12.0

if (x < twelve) then
  fn_val = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) /  &
           ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
  return
end if

!     Calculation for X >= 12.0

y = dlog(x)
fn_val = x * (y - one) - half * y + alr2pi
if (x > xlge) return
x1 = one / x
x2 = x1 * x1
fn_val = fn_val + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) /  &
         ((x2 + r4(5))*x2 + r4(4))
return
end function
