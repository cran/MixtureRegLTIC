module turnbull_est
  implicit none
  ! Global variables
  integer,save::nobs,njump,iterative,istruncation
  real(kind=8),save::tolerance=1.0d-06,llf,inf=999999.0d0
  real(kind=8),allocatable,save::yc(:,:),yt(:,:),weight(:)
  real(kind=8),allocatable,save::q_jump(:),p_jump(:),s_jump(:),alpha(:,:),beta(:,:)
contains
!==============================================================================================
! subroutine: turnbull_jump()
! Find the jump intervals of Turnbull (1976) method corrected by Frydram(1994) and Alioum(1996)
!==============================================================================================
subroutine turnbull_jump()
  ! Global variables: istruncation,nobs,njump,yc(:,:),yt(:,:),q_jump(:),p_jump
  implicit none
  integer::i,dimLR,numjump,leftcensordata=0,rightcensordata=0
  real(kind=8),allocatable::LR(:,:),q(:),p(:)

  ! Calculate the needed dimensions for LR
  dimLR=0
  do i=1,nobs
    if (yc(i,3)==2.0d0) then
      dimLR=dimLR+1
    else if (yc(i,3)==4.0d0) then
      dimLR=dimLR+1
    else if (yc(i,3)==3.0d0) then
      dimLR=dimLR+2
    else if (yc(i,3)==1.0d0) then
      dimLR=dimLR+2
    end if
    if (istruncation == 1) then
      if (yt(i,3)==2.0d0) then
        dimLR=dimLR+1
      else if (yt(i,3)==4.0d0) then
        dimLR=dimLR+1
      else if (yt(i,3)==3.0d0) then
        dimLR=dimLR+2
      end if
    end if
  end do
  do i=1,nobs
    if (yc(i,3)==2.0d0) then
      rightcensordata=1
    end if
    if (yc(i,3)==4.0d0) then
      leftcensordata=1
    end if
  end do
  if(leftcensordata==1) then
    dimLR=dimLR+1  ! For -inf end point of left censored data
  end if
  if(rightcensordata==1) then
    dimLR=dimLR+1  ! For inf end point of right censored data
  end if

  ! Create matrix LR with the time and corresponding index
  allocate(LR(dimLR,2))
  dimLR=0
  do i=1,nobs
    if (yc(i,3)==2.0d0) then
      dimLR=dimLR+1
      LR(dimLR,1)=yc(i,1)
      LR(dimLR,2)=3.0d0
    else if (yc(i,3)==4.0d0) then
      dimLR=dimLR+1
      LR(dimLR,1)=yc(i,2)
      LR(dimLR,2)=2.0d0
    else if (yc(i,3)==3.0d0) then
      dimLR=dimLR+1
      LR(dimLR,1)=yc(i,1)
      LR(dimLR,2)=3.0d0
      dimLR=dimLR+1
      LR(dimLR,1)=yc(i,2)
      LR(dimLR,2)=2.0d0
    else if (yc(i,3)==1.0d0) then
      dimLR=dimLR+1
      LR(dimLR,1)=yc(i,1)
      LR(dimLR,2)=1.0d0
      dimLR=dimLR+1
      LR(dimLR,1)=yc(i,2)
      LR(dimLR,2)=2.0d0
    end if
    if (istruncation == 1) then
      if (yt(i,3)==2.0d0) then
        dimLR=dimLR+1
        LR(dimLR,1)=yt(i,1)
        LR(dimLR,2)=2.0d0
      else if (yt(i,3)==4.0d0) then
        dimLR=dimLR+1
        LR(dimLR,1)=yt(i,2)
        LR(dimLR,2)=3.0d0
      else if (yt(i,3)==3.0d0) then
        dimLR=dimLR+1
        LR(dimLR,1)=yt(i,1)
        LR(dimLR,2)=2.0d0
        dimLR=dimLR+1
        LR(dimLR,1)=yt(i,2)
        LR(dimLR,2)=3.0d0
      end if
    end if
  end do
  if(leftcensordata==1) then
    dimLR=dimLR+1
    LR(dimLR,1)=-inf
    LR(dimLR,2)=3.0d0
  end if
  if(rightcensordata==1) then
    dimLR=dimLR+1
    LR(dimLR,1)=inf
    LR(dimLR,2)=2.0d0
  end if

  ! Sort the above matrix by 2(Increasing) and 1(Increasing) columns
  call matrix_sort(LR,dimLR,2,2,1)
  call matrix_sort(LR,dimLR,2,1,1)

  ! Adjustment for exact data
  do i=1,dimLR
    if(LR(i,2) == 1.0d0) then
      LR(i,2)=3.0d0
    end if
  end do

  ! Calculate jump number and create jump intervals
  allocate(q(dimLR))
  allocate(p(dimLR))
  numjump=0
  do i=1,(dimLR-1)
    if (LR(i,2) > LR(i+1,2)) then
      numjump=numjump+1
      q(numjump)=LR(i,1)
      p(numjump)=LR(i+1,1)
    end if
  end do
  deallocate(LR)

  ! Assign numjump, q_jump(:) and p_jump(:)
  njump=numjump
  allocate(q_jump(njump),p_jump(njump))
  q_jump(:)=q(1:njump)
  p_jump(:)=p(1:njump)

  deallocate(q)
  deallocate(p)

  return
end subroutine
!==============================
! subroutine: matrix_sort()
! To sort matrix by some column
!==============================
subroutine matrix_sort(matrix,n,m,sortcol,isincrease)
  implicit none
  integer::i,j,n,m,sortcol,isincrease
  real(kind=8)::matrix(n,m),tmp(m)

  if (isincrease==1) then ! increase
    do i=n,2,-1
      do j=1,i-1,1
        if (matrix(j,sortcol) > matrix(j+1,sortcol)) then
          tmp=matrix(j,:)
          matrix(j,:)=matrix(j+1,:)
          matrix(j+1,:)=tmp
        end if
      end do
    end do
  else ! decrease
    do i=n,2,-1
      do j=1,i-1,1
        if (matrix(j,sortcol) < matrix(j+1,sortcol)) then
          tmp=matrix(j,:)
          matrix(j,:)=matrix(j+1,:)
          matrix(j+1,:)=tmp
        end if
      end do
    end do
  end if

  return
end subroutine

!==================================================================================
! subroutine: turnbull_self()
! Find the mass of jump intervals with self-consistent algorithm of Turnbull (1976)
!==================================================================================
subroutine turnbull_self()
  ! Global Variables: tolerance,nobs,njump,yc(:,:),yt(:,:),q_jump(:),p_jump(:),s_jump(:)
  implicit none
  integer::i,j,k,m,isconvergence
  real(kind=8)::s_front(njump),s_rear(njump),ups(njump),downs(njump),llf_front,llf_rear
  real(kind=8)::buffer_alpha,buffer_beta,diffnorm,frontnorm,diffllf,difference
  character(len=10)::typeconvergence

  ! Calculate initial value of s_jump(:) => s_front(:)
  allocate(s_jump(njump))
  s_front=1.0d0/njump

  ! Calculate alpha(:,:) and beta(:,:)
  allocate(alpha(nobs,njump))
  alpha(:,:)=0.0d0
  call alpha_beta(yc,alpha)
  allocate(beta(nobs,njump))
  if (istruncation == 1) then
    beta(:,:)=0.0d0
    call alpha_beta(yt,beta)
  else
    beta(:,:)=1.0d0
  end if

  ! Caculate rear s(:) from former s(:)
  m=njump
  isconvergence=0
  iterative=0
  llf_front=0.0d0
  typeconvergence='norm' ! 'llf' or 'norm'
  do 99999 while(isconvergence /= 1)
    iterative=iterative+1
    ups(:)=0.0d0 ! used to estimate s_i
    downs(:)=0.0d0 ! used to estimate s_i

    ! Calculate rear s(:) from former s(:)
    do i=1,nobs
      ! Calculate sum(alpha(i,k)*s_front(k))
      buffer_alpha=0.0d0
      do k=1,m
        buffer_alpha=buffer_alpha+alpha(i,k)*s_front(k)
      end do
      ! Calculate sum(beta(i,k)*s_front(k))
      buffer_beta=0.0d0
      do k=1,m
        buffer_beta=buffer_beta+beta(i,k)*s_front(k)
      end do
      do k=1,m
        ups(k)=ups(k)+weight(i)*alpha(i,k)*s_front(k)/buffer_alpha &
                     +(1.0d0-beta(i,k))*weight(i)*s_front(k)/buffer_beta
        downs(k)=downs(k)+weight(i)*1.0d0/buffer_beta
      end do
    end do
    do k=1,m
      s_rear(k)=ups(k)/downs(k)
    end do

    ! Calculate log-likelihood function
    llf_rear=0.0d0
    do i=1,nobs
      ! Calculate sum(alpha(i,k)*s_rear(k))
      buffer_alpha=0.0d0
      do k=1,m
        buffer_alpha=buffer_alpha+alpha(i,k)*s_rear(k)
      end do
      ! Calculate sum(beta(i,k)*s_rear(k))
      buffer_beta=0.0d0
      do k=1,m
        buffer_beta=buffer_beta+beta(i,k)*s_rear(k)
      end do

      llf_rear=llf_rear+log(buffer_alpha)-log(buffer_beta)
    end do

    ! To see if convergence
    if(typeconvergence == 'llf') then
      diffllf=abs(llf_rear-llf_front)
      difference=diffllf
    else
      diffnorm=0.0d0
      frontnorm=0.0d0
      do k=1,m
        diffnorm=diffnorm+(s_rear(k)-s_front(k))**2.0d0
        frontnorm=frontnorm+s_front(k)**2.0d0
      end do
      diffnorm=dsqrt(diffnorm)
      frontnorm=dsqrt(frontnorm)
      difference=diffnorm/frontnorm
    end if
    if (difference < tolerance) then
      isconvergence=1
    else
      s_front=s_rear
      llf_front=llf_rear
    end if

  99999 continue
  s_jump=s_rear
  llf=llf_rear

  deallocate(alpha,beta,yc,yt,weight)
  return
end subroutine
!======================================
! subroutine: alpha_beta()
! To calculate alpha(:,:) and beta(:,:)
! To be called by initial()
!======================================
subroutine alpha_beta(y,index_matrix)
  ! Global Variables: nobs,njump,q_jump(:),p_jump(:)
  implicit none
  integer::i,j,front,rear
  real(kind=8)::y(nobs,3),sq(njump),sp(njump),index_matrix(nobs,njump)

  ! Assign Variables sq(:) and sp(:) from q_jump(:) and p_jump(:)
  sq(:)=q_jump(:)
  sp(:)=p_jump(:)

  ! Calculate index_matrix(:,:)
  do i=1,nobs
    front=0
    rear=0
    j=1
    do while(front==0 .and. j<=njump)
      if (y(i,1) <= sq(j) .and. sp(j) <= y(i,2)) then
        if ((sq(j) == sp(j)) .and. (sq(j)==y(i,1)) .and. (sq(j) /= y(i,2))) then
          j=j+1
        else
          front=j
        end if
      else
        j=j+1
      end if
    end do
    j=njump
    do while(rear == 0 .and. j>=1)
      if (y(i,1) <= sq(j) .and. sp(j) <= y(i,2)) then
        rear=j
      else
        j=j-1
      end if
    end do
    if (front/=0 .and. rear/=0) then
      index_matrix(i,front:rear)=1
    end if
  end do
  return
end subroutine

!======================
! subroutine: meanvar()
!======================
subroutine meanvar(ey,vary,skewy,p0,p2)
  ! Global Variables: iterative,tolerance,njump,q_jump(:),p_jump(:),s_jump(:)
  implicit none
  real(kind=8)::sur(njump),smin,smax
  integer::i,k,njumpd,first,last
  real(kind=8),allocatable::surd(:)
  real(kind=8)::ey,ey2,ey3,vary,skewy,p0,p2

  ! Calculate mixture survival function and cured susceptibility probability
  p0=s_jump(1)
  p2=s_jump(njump)
  sur(1)=1-s_jump(1)
  do i=2,njump
    sur(i)=sur(i-1)-s_jump(i)
  end do
  if (p_jump(1) /= 0.0d0 .and. p_jump(1) /= -inf) then
    smin=1.0d0
  else
    smin=sur(1)
  end if
  if (p_jump(njump) /= inf) then
    smax=0.0d0
  else
    if (njump /= 1) then
      smax=sur(njump-1)
    else
      smax=0.0d0
    end if
  end if
  first=1
  njumpd=njump
  if (smin /= 1.0d0) then
    njumpd=njumpd-1
    first=2
  end if
  last=njump
  if (smax /= 0.0d0) then
    njumpd=njumpd-1
    last=njump-1
  end if

  ! Calculate susceptibility survival function
  allocate(surd(njumpd))
  surd(1:njumpd)=(sur(first:last)-smax)/(smin-smax)

  ! Calculate ey, ey2, ey3
  ey=0.0d0
  ey2=0.0d0
  ey3=0.0d0
  if (njumpd == 1) then
    if (first /= last) then
      k=first
      if ((p_jump(k) <= 0.0d0)) then
        ey=ey-1.0d0*(0.0d0-p_jump(k)**1.0d0)
        ey2=ey2-1.d0*(0.0d0-p_jump(k)**2.0d0)
        ey3=ey3-1.0d0*(0.0d0-p_jump(k)**3.0d0)
      else if  ((p_jump(k) >= 0.0d0)) then
        ey=ey+1.0d0*(p_jump(k)**1.0d0)
        ey2=ey2+1.0d0*(p_jump(k)**2.0d0)
        ey3=ey3+1.0d0*(p_jump(k)**3.0d0)
      end if
    end if
  else
    do i=1,(njumpd-1)
      k=i+first-1
      if ((p_jump(k) <= 0.0d0) .and. (p_jump(k+1) <= 0.0d0)) then
        ey=ey-(1.0d0-surd(i))*(p_jump(k+1)**1.0d0-p_jump(k)**1.0d0)
        ey2=ey2-(1.0d0-surd(i))*(p_jump(k+1)**2.0d0-p_jump(k)**2.0d0)
        ey3=ey3-(1.0d0-surd(i))*(p_jump(k+1)**3.0d0-p_jump(k)**3.0d0)
      else if ((p_jump(i) >= 0.0d0) .and. (p_jump(i+1) >= 0.0d0)) then
        if (i == 1) then
          ey=ey+1.0d0*(p_jump(k)**1.0d0)
          ey2=ey2+1.0d0*(p_jump(k)**2.0d0)
          ey3=ey3+1.0d0*(p_jump(k)**3.0d0)
          ey=ey+surd(i)*(p_jump(k+1)**1.0d0-p_jump(k)**1.0d0)
          ey2=ey2+surd(i)*(p_jump(k+1)**2.0d0-p_jump(k)**2.0d0)
          ey3=ey3+surd(i)*(p_jump(k+1)**3.0d0-p_jump(k)**3.0d0)
        else
          ey=ey+surd(i)*(p_jump(k+1)**1.0d0-p_jump(k)**1.0d0)
          ey2=ey2+surd(i)*(p_jump(k+1)**2.0d0-p_jump(k)**2.0d0)
          ey3=ey3+surd(i)*(p_jump(k+1)**3.0d0-p_jump(k)**3.0d0)
        end if
      else
        ey=ey-(1.0d0-surd(k))*(0.0d0-p_jump(k)**1.0d0)
        ey2=ey2-(1.0d0-surd(k))*(0.0d0-p_jump(k)**2.0d0)
        ey3=ey3-(1.0d0-surd(k))*(0.0d0-p_jump(k)**3.0d0)
        ey=ey+surd(k)*(p_jump(k+1)**1.0d0-0.0d0)
        ey2=ey2+surd(k)*(p_jump(k+1)**2.0d0-0.0d0)
        ey3=ey3+surd(k)*(p_jump(k+1)**3.0d0-0.0d0)
      end if
    end do
  end if

  ! Calculate vary And skewy
  vary=ey2-(ey**2.0d0)
  skewy=ey3-3.0d0*ey*vary-(ey**3.0d0)

  return
end subroutine

!=====================
! subroutine: output()
!=====================
subroutine output(dist,ndist,mdist)
  ! Global Variables: iterative,tolerance,njump,q_jump(:),p_jump(:),s_jump(:)
  implicit none
  character(len=20)::out
  integer::i,ierror,ndist,mdist
  real(kind=8)::sur(njump),dist(ndist,mdist)

  ! Calculate survival function and write to file
  sur(1)=1.0d0-s_jump(1)
  do i=2,njump
    sur(i)=sur(i-1)-s_jump(i)
  end do

  ! Output with right point of intervals
!  ierror=-1
!  out="$selfname$.out"
!  open(unit=12,file=out,status='unknown',iostat=ierror)
!  do i=1,njump
!    write(12,3) q_jump(i),p_jump(i),p_jump(i),sur(i)
!  end do
!  3 format(1x,f15.4,1x,f15.4,1x,f15.4,1x,f15.4)
!  close(12)

  do i=1,njump
    dist(i,1)=q_jump(i)
    dist(i,2)=p_jump(i)
    dist(i,3)=p_jump(i)
    dist(i,4)=sur(i)
  end do

  deallocate(q_jump,p_jump,s_jump)

!  ierror=-1
!  open(unit=2,file="turnbullsA.txt",status='unknown',iostat=ierror)
!  write(2,*) njump,ndist
!  close(2)

  return
end subroutine

!=======================
! subroutine: turnbull()
!=======================
subroutine turnbull(ey,vary,skewy,p0,p2,survtime,nsurvtime,msurvtime,istruncation0,weight0,dist,ndist,mdist)
  ! Global variables: istruncation,nobs,njump,yc(:,:),yt(:,:),weight(:)
  ! Functions: turnbull_jump(), turnbull_self(), meanvar(), output()
  implicit none
  integer::nsurvtime,msurvtime,istruncation0,ndist,mdist
  real(kind=8)::survtime(nsurvtime,msurvtime),weight0(nsurvtime),dist(ndist,mdist)
  real(kind=8)::ey,vary,skewy,p0,p2
  integer::i,nvar,negativey=0
  real(kind=8),allocatable::origdata(:,:)

  allocate(origdata(nsurvtime,msurvtime))
  origdata=survtime

  istruncation=istruncation0
  nobs=nsurvtime
  nvar=msurvtime

  ! Assign observation responses to yc(:,:), yt(:,:) (assign 'yc' and 'yt')
  allocate(yc(nobs,3),yt(nobs,3))
  yc(:,1:3)=origdata(:,1:3)
  if (istruncation==1) then
    yt(:,1:3)=origdata(:,4:6)
  end if
  do i=1,nobs
    if (yc(i,3)==2.0d0) then
      yc(i,2)=inf
    else if (yc(i,3)==4.0d0) then
      yc(i,1)=-inf
    end if
    if (istruncation==1) then
      if (yt(i,3)==2.0d0) then
        yt(i,2)=inf
      else if (yt(i,3)==4.0d0) then
        yt(i,1)=-inf
      end if
    end if
  end do
  deallocate(origdata)

  allocate(weight(nobs))
  weight=weight0

! For nagative yc(:),yt(:)
  do i=1,nobs
    if ((yc(i,3)==1 .or. yc(i,3)==2 .or. yc(i,3)==3) .and. yc(i,1)<0) then
      negativey=1
      goto 3
    end if
  end do
3 do i=1,nobs
    if ((yc(i,3)==1 .or. yc(i,3)==3 .or. yc(i,3)==4) .and. yc(i,2)<0) then
      negativey=1
      goto 6
    end if
  end do
6 call turnbull_jump()
  call turnbull_self()

! For nagative yc(:),yt(:)
  if (negativey==0) then
    if (q_jump(1) == -inf) then
      q_jump(1)=0.0d0
    end if
  end if

  call meanvar(ey,vary,skewy,p0,p2)

  call output(dist,ndist,mdist) ! For turnbull estimator only

  return
end subroutine
end module

!========================
! subroutine: turnbulls()
!========================
! For DLL
subroutine turnbulls(ey,vary,skewy,p0,p2,survtime,nsurvtime,msurvtime,istruncation0,weight0,dist,ndist,mdist,tolerance0)
!DEC$ ATTRIBUTES DLLEXPORT:: turnbulls
!DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'turnbulls_' :: turnbulls
  use turnbull_est
  implicit none
  integer::nsurvtime,msurvtime,istruncation0,ndist,mdist
  real(kind=8)::survtime(nsurvtime,msurvtime),weight0(nsurvtime),dist(ndist,mdist)
  real(kind=8)::ey,vary,skewy,p0,p2,tolerance0

  tolerance=tolerance0

  call turnbull(ey,vary,skewy,p0,p2,survtime,nsurvtime,msurvtime,istruncation0,weight0,dist,ndist,mdist)

  return
end subroutine
