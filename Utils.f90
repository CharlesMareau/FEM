! **********************************************************************
! *** module math                                                    ***
! *** contains the subroutines for basic operations                  ***
! **********************************************************************
module utils
 implicit none
 real(8), parameter, public :: pi=3.141592653589793d0
contains
! **********************************************************************
! *** function radian                                                ***
! *** converts an angle from degree to radian                        ***
! ***                                                                ***
! *** input :  degree, angle to be converted into radian             ***
! **********************************************************************
function radian(degree)
 implicit none
 real(8) :: degree,radian
 radian=degree*pi/1.8d2
end function
! **********************************************************************
! *** function degree                                                ***
! *** converts an angle from radian to degree                        ***
! ***                                                                ***
! *** input :  radian, angle to be converted into degree             ***
! **********************************************************************
function degree(radian)
 implicit none
 real(8) :: degree,radian
 degree=radian*1.8d2/pi
end function
! **********************************************************************
! *** subroutine ludcmp                                              ***
! *** performs the LU decomposition of a matrix                      ***
! ***                                                                ***
! *** input :  np, physical size of the matrix                       ***
! ***          n, actual size                                        ***
! ***                                                                ***
! *** in/output : a, matrix before/after LU decomposition            ***
! ***                                                                ***
! *** output : indx, permutation vector                              ***
! ***          d, ±1 depending on whether the number of row          ***
! ***             interchanges was even or odd                       ***
! **********************************************************************
subroutine ludcmp(a,n,np,indx,d)
 implicit none
 integer, intent(in) :: n,np
 integer, intent(out) :: indx(n)
 integer, parameter :: nmax=500
 real(8), intent(inout) :: a(np,np)
 real(8), intent(out) :: d
 real(8), parameter :: tiny=1.d-20
 integer :: i,imax,j,k
 real(8) :: aamax,dum,sum,vv(nmax)
 d=1.d+0
 do i=1,n
  do j=1,n
   if (abs(a(i,j)).gt.aamax) then
    aamax=abs(a(i,j))
   endif
  enddo
  if (aamax.eq.0.d+0) then
   write(*,*) 'singular matrix in ludcmp'
   stop
  endif
  vv(i)=1.d+0/aamax
 enddo
 do j=1,n
  do i=1,j-1
   sum=a(i,j)
   do k=1,i-1
    sum=sum-a(i,k)*a(k,j)
   enddo
   a(i,j)=sum
  enddo
  aamax=0.d+0
  do i=j,n
   sum=a(i,j)
   do k=1,j-1
    sum=sum-a(i,k)*a(k,j)
   enddo
   a(i,j)=sum
   dum=vv(i)*abs(sum)
   if (dum.ge.aamax) then
    imax=i
    aamax=dum
   endif
  enddo
  if (j.ne.imax)then
   do k=1,n
    dum=a(imax,k)
    a(imax,k)=a(j,k)
    a(j,k)=dum
   enddo
   d=-d
   vv(imax)=vv(j)
  endif
 indx(j)=imax
 if(a(j,j).eq.0.d+0) then
  a(j,j)=tiny
 endif
 if(j.ne.n)then
  dum=1.d+0/a(j,j)
  do i=j+1,n
   a(i,j)=a(i,j)*dum
  enddo
 endif
 enddo
 return
end subroutine ludcmp
! **********************************************************************
! *** subroutine lubksb                                              ***
! *** gives the solution of a linear system of equations after a LU  ***
! *** decomposition                                                  ***
! ***                                                                ***
! *** input :  a, square matrix (after LU decomposition)             ***
! ***          np, physical size of the matrix                       ***
! ***          n, actual size                                        ***
! ***          indx, permutation vector returned by ludcmp           ***
! ***                                                                ***
! *** in/output : b, problem/solution vector                         ***
! **********************************************************************
subroutine lubksb(a,n,np,indx,b)
 implicit none
 integer, intent(in) :: n,np,indx(n)
 real(8), intent(in) :: a(np,np)
 real(8), intent(inout) :: b(n)
 integer :: i,ii,j,ll
 real(8) :: sum
 ii=0 
 do i=1,n
  ll=indx(i)
  sum=b(ll)
  b(ll)=b(i)
  if (ii.ne.0)then
   do j=ii,i-1
    sum=sum-a(i,j)*b(j)
   enddo
  else if (sum.ne.0.d+0) then
   ii=i
  endif 
  b(i)=sum
 enddo
 do i=n,1,-1
  sum=b(i)
  do j=i+1,n
   sum=sum-a(i,j)*b(j)
  enddo
  b(i)=sum/a(i,i)
 enddo
 return 
end subroutine lubksb
! **********************************************************************
! *** function invmat                                                ***
! *** inverts a square matrix by using a LU decomposition            ***
! ***                                                                ***
! *** input :  a, matrix to be inverted                              ***
! ***          np, physical size of the matrix                       ***
! ***          n, actual size of the matrix                          ***
! **********************************************************************
function invmat(a,n,np)
 implicit none
 integer, intent(in) :: n,np
 real(8), intent(in) :: a(np,np)
 real(8) :: invmat(np,np)
 integer :: indx(np),i,j
 real(8) :: y(np,np),d,amax,dum

 amax=0.d0
 do i=1,n
  do j=1,n
   dum=abs(a(i,j))
   if(dum.gt.amax) then
    amax=dum
   endif
   y(i,j)=0.d+0
  enddo
  y(i,i)=1.d+0
 enddo
 invmat=a/amax
 call ludcmp(invmat,n,np,indx,d)
 do j=1,n
  call lubksb(invmat,n,np,indx,y(1,j))
 enddo
 invmat=y/amax
 
 return
end function invmat
! **********************************************************************
! *** function pseudoinverse                                         ***
! *** inverts a rectangular matrix by using a SVD decomposition      ***
! ***                                                                ***
! *** input :  sym, tag to indicate whether the matrix is symmetric  ***
! ***               or not                                           ***
! ***          a, matrix to be inverted                              ***
! ***          m, number of lines                                    ***
! ***          n, number of columns                                  ***
! ***          rcond, threshold for singular value detection         ***
! **********************************************************************
function pseudoinverse(a,m,n,rcond)
 implicit none
 integer, intent(in) :: m,n
 integer :: k,lwork
 real(8), intent(in)  :: a(m,n)
 real(8), intent(in)  :: rcond
 real(8) :: pseudoinverse(n,m)
 integer :: rank, info
 integer :: i, j, l
 real(8) :: work_query(1), tol
 real(8), allocatable :: ac(:,:),u(:,:),vt(:,:),s(:),work(:)

 k=max(m,n)
 lwork=max(1,3*k+max(m,n),5*k)
 allocate(ac(m,n),u(m,k),vt(k,n),s(k),work(lwork))

 if (rcond.lt.0.0d0) then
  write(*,'("ERROR IN PSEUDOINVERSE. PROGRAM WILL STOP.")')
  return
 end if
 ac=a
 call dgesvd('S', 'S', m, n, ac, m, s, u, m, vt, k, work, lwork, info)
 if (info.ne.0) then
  rank = 0
  return
 end if
 tol = rcond * S(1)
 pseudoinverse = 0.0d0
 rank = 0
 do l = 1, k
  if (s(l) > tol) then
   rank = rank + 1
   do j = 1, m
    do i = 1, n
     pseudoinverse(i,j) = pseudoinverse(i,j)+vt(l,i)*u(j,l)/s(l)
    end do
   end do
  end if
 end do
 deallocate(ac,u,vt,s,work)
end function
! **********************************************************************
! *** function inverse                                               ***
! *** inverts a square matrix by using a LU decomposition            ***
! ***                                                                ***
! *** input :  sym, tag to indicate whether the matrix is symmetric  ***
! ***               or not                                           ***
! ***          a, matrix to be inverted                              ***
! ***          n, actual size of the matrix                          ***
! **********************************************************************
function inverse(uplo,a,n)
 implicit none
 integer, intent(in) :: n
 real(8), intent(in) :: a(n,n)
 real(8) :: inverse(n,n),work(n)
 integer :: ipiv(n),info,i,j
 character :: uplo
 inverse=a
 if (uplo.eq.'U'.or.uplo.eq.'L') then
  call dsytrf(uplo,n,inverse,n,ipiv,work,n,info)
  if (info.ne.0) then
   write(*,'("ERROR IN MATRIX FACTORIZATION. PROGRAM WILL STOP.")')
   stop
  endif
  call dsytri(uplo,n,inverse,n,ipiv,work,info)
  if (info.ne.0) then
   write(*,'("ERROR IN MATRIX INVERSION. PROGRAM WILL STOP.")')
   stop
  endif
  if (uplo.eq.'U') then
   do i=1,n
    do j=(i+1),n
     inverse(j,i)=inverse(i,j)
    enddo
   enddo
  else
   do i=1,n
    do j=(i+1),n
     inverse(i,j)=inverse(j,i)
    enddo
   enddo
  endif
 else
  call dgetrf(n,n,inverse,n,ipiv,info)
  if (info.ne.0) then
   write(*,'("ERROR IN MATRIX FACTORIZATION. PROGRAM WILL STOP.")')
   stop
  endif
  call dgetri(n,inverse,n,ipiv,work,n,info)
  if (info.ne.0) then
   write(*,'("ERROR IN MATRIX INVERSION. PROGRAM WILL STOP.")')
   stop
  endif
 endif
 return
end function inverse
! **********************************************************************
! *** function sqrmat                                                ***
! *** returns the square root of a square matrix                     ***
! ***                                                                ***
! *** input :  a, matrix                                             ***
! ***          np, physical size of the matrix                       ***
! ***          n, actual size of the matrix                          ***
! **********************************************************************
function sqrmat(a,n,np)
 implicit none
 integer n,np
 real(8) :: a(np,np)
 real(8) :: sqrmat(np,np)
 real(8) :: y(np,np),z(np,np),b(np,np),q(np,np)
 real(8), parameter :: toler=1.e-8
 integer :: iter
 integer, parameter :: niter=100
 real(8) :: error,nrm1,nrm2
 integer :: i,j,k

 y=a
 do i=1,n
  do j=1,n
   y(i,j)=a(i,j)
   z(i,j)=0.d0
  enddo
  z(i,i)=1.d0
 enddo
 error=toler+1.d0
 iter=0
 do while (error.gt.toler)
  do i=1,n
   do j=1,n
    sqrmat(i,j)=y(i,j)
   enddo
  enddo
  iter=iter+1
  a=invmat(y,n,np)
  b=invmat(z,n,np)
  do i=1,n
   do j=1,n
    y(i,j)=0.5d0*(y(i,j)+b(i,j))
    z(i,j)=0.5d0*(z(i,j)+a(i,j))
   enddo
  enddo
  nrm1=0.d0
  nrm2=0.d0
  do i=1,6
   do j=1,6
    nrm1=nrm1+(y(i,j)-sqrmat(i,j))**2.d0
    nrm2=nrm2+y(i,j)**2.d0
   enddo
  enddo
  nrm1=sqrt(nrm1)
  nrm2=sqrt(nrm2)
  error=nrm1/nrm2
  if (iter.gt.niter) then
   write(*,'("ERROR IN sqrmat. PROGRAM WILL STOP")')
   stop
  endif
 enddo

end function
! **********************************************************************
! *** function matxmat                                               ***
! *** performs the multiplication between two matrices               ***
! ***                                                                ***
! *** input :  a, right matrix                                       ***
! ***          b, left matrix                                        ***
! ***          n, left dimension of a                                ***
! ***          m, right dimension of a / left dimension of b         ***
! ***          p, right dimension of b                               ***
! **********************************************************************
function matxmat(a,b,n,m,p)
 implicit none
 integer, intent(in) :: n,m,p
 real(8), intent(in) :: a(n,m),b(m,p)
 real(8) :: matxmat(n,p)
 integer :: i,j,k

 
 do j=1,p
  do i=1,n
   matxmat(i,j)=0.d0
   do k=1,m
    matxmat(i,j)=matxmat(i,j)+a(i,k)*b(k,j)
   enddo
  enddo
 enddo
 
 return
end function matxmat
! **********************************************************************
! *** function matxvec                                               ***
! *** performs the multiplication of one matrix with a vector        ***
! ***                                                                ***
! *** input :  a dimension:(n,m) type:real(8)                        ***
! ***          b dimension:(m) type:real(8)                          ***
! ***          n  dimension:(1) type:integer                         ***
! ***          m  dimension:(1) type:integer                         ***
! ***                                                                ***
! *** output : matxvec dimension:(n) type:real(8)                    ***
! **********************************************************************
function matxvec(a,b,n,m)
 implicit none
 integer, intent(in) :: n,m
 real(8), intent(in) :: a(n,m),b(m)
 real(8) :: matxvec(n)
 integer :: i,j

 do i=1,n
  matxvec(i)=0.d0
  do j=1,m
   matxvec(i)=matxvec(i)+a(i,j)*b(j)
  enddo
 enddo
 
 return
end function matxvec
! **********************************************************************
! *** function euler                                                 ***
! *** calculates the rotation matrix from euler angles (given in     ***
! *** degrees)                                                       ***
! ***                                                                ***
! *** input :  ph, euler angle #1                                    ***
! ***          th, euler angle #2                                    ***
! ***          tm, euler angle #3                                    ***
! **********************************************************************
function euler(ph,th,tm)
 implicit none
 real(8), intent(in) :: ph,th,tm
 real(8) :: euler(3,3)
 real(8) :: cph,cth,ctm,sph,sth,stm

 sph=sin(ph*pi/180.d+0)
 cph=cos(ph*pi/180.d+0)
 sth=sin(th*pi/180.d+0)
 cth=cos(th*pi/180.d+0)
 stm=sin(tm*pi/180.d+0)
 ctm=cos(tm*pi/180.d+0)
 
 euler(1,1)=ctm*cph-sph*stm*cth
 euler(2,1)=-stm*cph-sph*ctm*cth
 euler(3,1)=sph*sth
 euler(1,2)=ctm*sph+cph*stm*cth
 euler(2,2)=-sph*stm+cph*ctm*cth
 euler(3,2)=-sth*cph
 euler(1,3)=sth*stm
 euler(2,3)=ctm*sth
 euler(3,3)=cth
 
 return
end function euler
! **********************************************************************
! *** subroutine angles                                              ***
! *** calculates the from euler angles (in degrees) from a rotation  ***
! *** matrix                                                         ***
! ***                                                                ***
! *** input :  rot, rotation matrix                                  ***
! ***                                                                ***
! *** output :  ph, euler angle #1                                   ***
! ***           th, euler angle #2                                   ***
! ***           tm, euler angle #3                                   ***
! **********************************************************************
subroutine angles(rot,ph,th,tm)
 implicit none
 real(8), intent(out) :: ph,th,tm
 real(8), intent(in) :: rot(3,3)
 real(8), parameter :: toler=1.d-7

 th=acos(rot(3,3))
 if (th.gt.toler.and.th.lt.(pi-toler)) then ! general case
  ph=atan2(rot(3,1),-rot(3,2))
  tm=atan2(rot(1,3),rot(2,3))
 endif 
 if (th.le.toler) then ! specific case theta=0
  tm=0.d0
  ph=atan2(-rot(2,1),rot(1,1))
 endif 
 if (th.ge.(pi-toler)) then ! specific case theta=pi
  tm=0.d0
  ph=-atan2(-rot(2,1),rot(1,1))
 endif
  
 th=th*180.d0/pi
 ph=ph*180.d0/pi
 tm=tm*180.d0/pi
 
 return
end subroutine angles
! **********************************************************************
! *** function  pythag                                               ***
! *** calculates (a^2+b^2)^{1/2} with limited risk of overflow       ***
! ***                                                                ***
! *** input : a, first real number                                   ***
! ***         b, second real number                                  ***
! **********************************************************************
function pythag(a,b)
 implicit none
 real(8), intent(in) :: a,b
 real(8) :: pythag
 real(8) :: absa,absb
 
 absa=abs(a)
 absb=abs(b)
 if(absa.gt.absb)then
  pythag=absa*sqrt(1.d+0+(absb/absa)**2.d+0)
 else
  if(absb.eq.0.d+0)then
   pythag=0.d+0
  else
   pythag=absb*sqrt(1.d+0+(absa/absb)**2.d+0)
  endif
 endif
 return
end function pythag
! **********************************************************************
! *** function  bump                                                 ***
! *** returns the value of the bump function                         ***
! ***                                                                ***
! *** input : x, variable                                            ***
! ***         w, parameter                                           ***
! **********************************************************************
function bump(x,w)
 implicit none
 real(8), intent(in) :: x,w
 real(8) :: bump
 
 if (abs(x).lt.w) then
  bump=exp(w**2.d0/(x**2.d0-w**2.d0))
 else
  bump=0.d0
 endif
 return
end function bump
! **********************************************************************
! *** subroutine svdcmp                                              ***
! *** performs a SVD decomposition of a matrix a=u.w.v^t             ***
! ***                                                                ***
! *** input :  mp, physical size                                     ***
! ***          np, physical size                                     ***
! ***          m, actual size                                        ***
! ***          n, actual size                                        ***
! ***                                                                ***
! *** in/output : a, matrix before SVD decomposition (becomes u)     ***
! ***                                                                ***
! *** output : w, SVD middle vector                                  ***
! ***          v, SVD right matrix                                   ***
! **********************************************************************
subroutine svdcmp(a,m,n,mp,np,w,v)
 implicit none
 integer, intent(in) :: m,mp,n,np
 real(8), intent(inout) :: a(mp,np)
 real(8), intent(out) :: v(np,np),w(np)
 integer, parameter :: nmax=500
 integer :: i,its,j,jj,k,l,nm
 real(8) :: anorm,c,f,g,h,s,scale,x,y,z,rv1(nmax)
 
 g=0.d+0
 scale=0.d+0
 anorm=0.d+0
 do 25 i=1,n
  l=i+1
  rv1(i)=scale*g
  g=0.d+0
  s=0.d+0
  scale=0.d+0
  if(i.le.m)then
   do 11 k=i,m
    scale=scale+abs(a(k,i))
   11 continue
   if(scale.ne.0.d+0)then
    do 12 k=i,m
     a(k,i)=a(k,i)/scale
     s=s+a(k,i)*a(k,i)
    12 continue
    f=a(i,i)
    g=-sign(sqrt(s),f)
    h=f*g-s
    a(i,i)=f-g
    do 15 j=l,n
     s=0.d+0
     do 13 k=i,m
      s=s+a(k,i)*a(k,j)
     13 continue
     f=s/h
     do 14 k=i,m
      a(k,j)=a(k,j)+f*a(k,i)
     14 continue
    15 continue
    do 16 k=i,m
     a(k,i)=scale*a(k,i)
    16 continue
   endif
  endif
  w(i)=scale *g
  g=0.d+0
  s=0.d+0
  scale=0.+0
  if((i.le.m).and.(i.ne.n))then
   do 17 k=l,n
    scale=scale+abs(a(i,k))
   17 continue
   if(scale.ne.0.d+0)then
    do 18 k=l,n
     a(i,k)=a(i,k)/scale
     s=s+a(i,k)*a(i,k)
    18 continue
    f=a(i,l)
    g=-sign(sqrt(s),f)
    h=f*g-s
    a(i,l)=f-g
    do 19 k=l,n
     rv1(k)=a(i,k)/h
    19 continue
    do 23 j=l,m
     s=0.d+0
     do 21 k=l,n
      s=s+a(j,k)*a(i,k)
     21 continue
     do 22 k=l,n
      a(j,k)=a(j,k)+s*rv1(k)
     22 continue
    23 continue
    do 24 k=l,n
     a(i,k)=scale*a(i,k)
    24 continue
   endif
  endif
  anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
 25 continue
 do 32 i=n,1,-1
  if(i.lt.n)then
   if(g.ne.0.d+0)then
    do 26 j=l,n
     v(j,i)=(a(i,j)/a(i,l))/g
    26 continue
    do 29 j=l,n
     s=0.d+0
     do 27 k=l,n
      s=s+a(i,k)*v(k,j)
     27 continue
     do 28 k=l,n
      v(k,j)=v(k,j)+s*v(k,i)
     28 continue
    29 continue
   endif
   do 31 j=l,n
    v(i,j)=0.d+0
    v(j,i)=0.d+0
   31 continue
  endif
  v(i,i)=1.d+0
  g=rv1(i)
  l=i
 32 continue
 do 39 i=min(m,n),1,-1
  l=i+1
  g=w(i)
  do 33 j=l,n
   a(i,j)=0.d+0
  33 continue
  if(g.ne.0.d+0)then
   g=1.d+0/g
   do 36 j=l,n
    s=0.d+0
    do 34 k=l,m
     s=s+a(k,i)*a(k,j)
    34 continue
    f=(s/a(i,i))*g
    do 35 k=i,m
     a(k,j)=a(k,j)+f*a(k,i)
    35 continue
   36 continue
   do 37 j=i,m
    a(j,i)=a(j,i)*g
   37 continue
  else
   do 38 j= i,m
    a(j,i)=0.d+0
   38 continue
  endif
  a(i,i)=a(i,i)+1.d+0
 39 continue
 do 49 k=n,1,-1
  do 48 its=1,30
   do 41 l=k,1,-1
    nm=l-1
    if((abs(rv1(l))+anorm).eq.anorm)  goto 2
    if((abs(w(nm))+anorm).eq.anorm)  goto 1
   41 continue
   1 c=0.d+0
   s=1.d+0
   do 43 i=l,k
    f=s*rv1(i)
    rv1(i)=c*rv1(i)
    if((abs(f)+anorm).eq.anorm) goto 2
     g=w(i)
     h=pythag(f,g)
     w(i)=h
     h=1.d+0/h
     c=(g*h)
     s=-(f*h)
     do 42 j=1,m
      y=a(j,nm)
      z=a(j,i)
      a(j,nm)=(y*c)+(z*s)
      a(j,i)=-(y*s)+(z*c)
     42 continue
    43 continue
    2 z=w(k)
    if(l.eq.k)then
     if(z.lt.0.d+0)then
      w(k)=-z
      do 44 j=1,n
       v(j,k)=-v(j,k)
      44 continue
     endif
     goto 3
    endif
    if(its.eq.30) write(*,*) 'no convergence in svdcmp'
    x=w(l)
    nm=k-1
    y=w(nm)
    g=rv1(nm)
    h=rv1(k)
    f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d+0*h*y)
    g=pythag(f,1.d+0)
    f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
    c=1.d+0
    s=1.d+0
    do 47 j=l,nm
     i=j+1
     g=rv1(i)
     y=w(i)
     h=s*g
     g=c*g
     z=pythag(f,h)
     rv1(j)=z
     c=f/z
     s=h/z
     f=(x*c)+(g*s)
     g=-(x*s)+(g*c)
     h=y*s
     y=y*c
     do 45 jj=1,n
      x=v(jj,j)
      z=v(jj,i)
      v(jj,j)=(x*c)+(z*s)
      v(jj,i)=-(x*s)+(z*c)
     45 continue
     z=pythag(f,h)
     w(j)=z
     if(z.ne.0.d+0)then
      z=1.d+0/z
      c=f*z
      s=h*z
     endif
     f=(c*g)+(s*y)
     x=-(s*g)+(c*y)
     do 46 jj=1,m
      y=a(jj,j)
      z=a(jj,i)
      a(jj,j)=(y*c)+(z*s)
      a(jj,i)=-(y*s)+(z*c)
     46 continue
    47 continue
    rv1(l)=0.d+0
    rv1(k)=f
    w(k)=x
   48 continue
  3 continue
 49 continue
 return
end subroutine
! **********************************************************************
! *** subroutine svbksb                                              ***
! *** performs a singular value backsubstitution to solve a linear   ***
! *** system of equations                                            ***
! ***                                                                ***
! *** input :  mp, physical size                                     ***
! ***          np, physical size                                     ***
! ***          m, actual size                                        ***
! ***          n, actual size                                        ***
! ***          u, SVD left matrix                                    ***
! ***          w, SVD middle vector                                  ***
! ***          v, SVD right matrix                                   ***
! ***          b, right hand side vector                             ***
! ***                                                                ***
! *** output : x, solution vector                                    ***
! **********************************************************************
subroutine svbksb(u,w,v,m,n,mp,np,b,x)
 integer, intent(in) :: m,mp,n,np
 real(8), intent(in) :: b(mp),u(mp,np),v(np,np),w(np)
 real(8), intent(out) :: x(np)
 integer, parameter :: nmax=500
 integer :: i,j,jj
 real(8) :: s,tmp(nmax)

 do 12 j=1,n
  s=0.d+0
  if(w(j).ne.0.d+0)then
   do 11 i=1,m
    s=s+u(i,j)*b(i)
   11 continue
   s=s/w(j)
  endif
  tmp(j)=s
 12 continue
 do 14 j=1,n
  s=0.d+0
  do 13 jj=1,n
   s=s+v(j,jj)*tmp(jj)
  13 continue
  x(j)=s
 14 continue
 return
end subroutine svbksb
! **********************************************************************
! *** subroutine jacobi                                              ***
! *** computes all eigenvalues and eigenvectors of a square matrix   ***
! ***                                                                ***
! *** input :  np, physical size                                     ***
! ***          n, actual size                                        ***
! ***          a, symmetric square matrix                            ***
! ***                                                                ***
! *** output : d, eigenvalues                                        ***
! ***          v, eigenvectors                                       ***
! **********************************************************************
subroutine jacobi(a,n,np,d,v)
 implicit none
 integer, intent(in) :: n,np
 real(8), intent(in) :: a(np,np)
 real(8), intent(out) :: d(np),v(np,np)
 integer, parameter :: nmax=500
 integer :: i,ip,iq,j,nrot
 real(8) c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),aa(np,np)
     
 do 12 ip=1,n
  do 11 iq=1,n
   v(ip,iq)=0.d+0
   aa(ip,iq)=a(ip,iq)
  11 continue
  v(ip,ip)=1.d+0
 12 continue
 do 13 ip=1,n
  b(ip)=aa(ip,ip)
  d(ip)=b(ip)
  z(ip)=0.d+0
 13 continue
 nrot=0
 do 24 i=1,50
  sm=0.d+0
  do 15 ip=1,n-1
   do 14 iq=ip+1,n
    sm=sm+abs(aa(ip,iq))
   14 continue
  15 continue
  if(sm.eq.0.d+0)return
  if(i.lt.4)then
   tresh=0.2d+0*sm/n**2.d+0
  else
   tresh=0.d+0
  endif
  do 22 ip=1,n-1
   do 21 iq=ip+1,n
    g=100.d+0*abs(aa(ip,iq))
    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
     aa(ip,iq)=0.d+0
    else if(abs(aa(ip,iq)).gt.tresh)then
     h=d(iq)-d(ip)
     if(abs(h)+g.eq.abs(h))then
      t=aa(ip,iq)/h
     else
      theta=0.5d+0*h/aa(ip,iq)
      t=1.d+0/(abs(theta)+sqrt(1.d+0+theta**2.d+0))
      if(theta.lt.0.d+0)t=-t
     endif
     c=1.d+0/sqrt(1.d+0+t**2.d+0)
     s=t*c
     tau=s/(1.d+0+c)
     h=t*aa(ip,iq)
     z(ip)=z(ip)-h
     z(iq)=z(iq)+h
     d(ip)=d(ip)-h
     d(iq)=d(iq)+h
     aa(ip,iq)=0.d+0
     do 16 j=1,ip-1
      g=aa(j,ip)
      h=aa(j,iq)
      aa(j,ip)=g-s*(h+g*tau)
      aa(j,iq)=h+s*(g-h*tau)
     16 continue
     do 17 j=ip+1,iq-1
      g=aa(ip,j)
      h=aa(j,iq)
      aa(ip,j)=g-s*(h+g*tau)
      aa(j,iq)=h+s*(g-h*tau)
     17 continue
     do 18 j=iq+1,n
      g=aa(ip,j)
      h=aa(iq,j)
      aa(ip,j)=g-s*(h+g*tau)
      aa(iq,j)=h+s*(g-h*tau)
     18 continue
     do 19 j=1,n
      g=v(j,ip)
      h=v(j,iq)
      v(j,ip)=g-s*(h+g*tau)
      v(j,iq)=h+s*(g-h*tau)
     19 continue
     nrot=nrot+1
    endif
   21 continue
  22 continue
  do 23 ip=1,n
   b(ip)=b(ip)+z(ip)
   d(ip)=b(ip)
   z(ip)=0.d+0
  23 continue
 24 continue
 write(*,'("TOO MANY ITERATIONS IN JACOBI. PROGRAM WILL STOP")') 
 stop
 return
end subroutine
! **********************************************************************
! *** function psdinvmat                                             ***
! *** inverts a matrix by using a SVD decomposition                  ***
! ***                                                                ***
! *** input :  a, matrix to be inverted                              ***
! ***          mp, physical leftt size of the matrix                 ***
! ***          np, physical right size of the matrix                 ***
! ***          m, actual left size of the matrix                     ***
! ***          n, actual right size of the matrix                    ***
! **********************************************************************
function psdinvmat(a,m,n,mp,np)
 implicit none
 integer, intent(in) :: m,n,mp,np
 real(8), intent(in) :: a(mp,np)
 real(8) :: w(np),v(np,np),u(mp,np)
 integer :: i,j,k
 real(8) :: psdinvmat(np,mp),wmin,wmax
 real(8), parameter :: small=1.d-8,little=1.d-10
 
 do i=1,m
  do j=1,n
   u(i,j)=a(i,j)
  enddo
 enddo
 
 call svdcmp(u,m,n,mp,np,w,v)
 
 wmax=-1.d32
 do j=1,n
  if (w(j).gt.wmax) then
   wmax=w(j)
  endif
 enddo
 wmin=small*wmax
 do j=1,n
  if (w(j).le.wmin.or.w(j).le.little) then
   w(j)=0.d0
  endif
 enddo
 
 do i=1,n
  do j=1,m
   psdinvmat(i,j)=0.d0
   do k=1,n
    if (w(k).gt.wmin) then
     psdinvmat(i,j)=psdinvmat(i,j)+v(i,k)*(1.d0/w(k))*u(j,k)
    endif
   enddo
  enddo
 enddo
 
end function psdinvmat
! **********************************************************************
! *** function iseven                                                ***
! *** returns true is the input integer number is even               ***
! ***                                                                ***
! *** input :  n, integer number                                     ***
! **********************************************************************
function iseven(n)
 implicit none
 integer, intent(in) :: n
 logical :: iseven

 if (mod(n,2).eq.0) then
  iseven=.true.
 else
  iseven=.false.
 endif
end function iseven
! **********************************************************************
! *** function avrg                                                  ***
! *** calculates the average value associated with a weighted list   ***
! *** of scalar values                                               ***
! ***                                                                ***
! *** input : n, number of elements                                  ***
! ***         w, list of weights                                     ***
! ***         v, list of scalar values                               ***
! **********************************************************************
function avrg(n,w,v)
 integer, intent(in) :: n
 real(8), intent(in) :: w(n),v(n)
 real(8) :: avrg
 integer :: i

 avrg=0.d0
 do i=1,n
  avrg=avrg+w(i)*v(i)
 enddo

end function avrg
! **********************************************************************
! *** function weight                                                ***
! *** calculates the weight when computing the non-local average of  ***
! *** some variableues                                               ***
! ***                                                                ***
! *** input : d, euclidean distance                                  ***
! ***         r, truncation radius                                   ***
! **********************************************************************
function weight(d,r)
 real(8) :: weight,d,r
 
 if (d.lt.r) then
  weight=15.d0/16.d0*(1-(d/r)**2.d0)**2.d0
 else
  weight=0.d0
 endif
! weight=exp(-0.5d0*(d/r)**2.d0)/r/sqrt(2.d0*pi)
end function
! **********************************************************************
! *** function fact                                                  ***
! *** return the factorial of an integer number                      ***
! ***                                                                ***
! *** input : n, integer number                                      ***
! **********************************************************************
function fact(n)
 integer :: n,fact,i
 
 if (n.eq.0) then
  fact=1
 else
  fact=1
  do i=2,n
   fact=fact*i
  enddo
 endif
 
end function
! **********************************************************************
! *** subroutine rank                                                ***
! *** returns the ranking in ascending order of a an array of real   ***
! *** numbers                                                        ***
! **********************************************************************
subroutine rank(n,arr,indx)
 integer :: n,indx(n),m,nstack
 real(8) :: arr(n)
 parameter (m=7,nstack=50)
 integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
 real(8) :: a
 
 do 11 j=1,n
  indx(j)=j
11 continue
 jstack=0
 l=1
 ir=n
1 if (ir-l.lt.m) then
  do 13 j=l+1,ir
   indxt=indx(j)
   a=arr(indxt)
   do 12 i=j-1,1,-1
    if(arr(indx(i)).le.a)goto 2
    indx(i+1)=indx(i)
12 continue
   i=0
2  indx(i+1)=indxt
13 continue
  if(jstack.eq.0)return
  ir=istack(jstack)
  l=istack(jstack-1)
  jstack=jstack-2
 else
  k=(l+ir)/2
  itemp=indx(k)
  indx(k)=indx(l+1)
  indx(l+1)=itemp
  if(arr(indx(l+1)).gt.arr(indx(ir)))then
   itemp=indx(l+1)
   indx(l+1)=indx(ir)
   indx(ir)=itemp
  endif
  if(arr(indx(l)).gt.arr(indx(ir)))then
   itemp=indx(l)
   indx(l)=indx(ir)
   indx(ir)=itemp
  endif
  if(arr(indx(l+1)).gt.arr(indx(l)))then
   itemp=indx(l+1)
   indx(l+1)=indx(l)
   indx(l)=itemp
  endif
  i=l+1
  j=ir
  indxt=indx(l)
  a=arr(indxt)
3 continue
  i=i+1
  if(arr(indx(i)).lt.a)goto 3
4 continue
  j=j-1
  if(arr(indx(j)).gt.a)goto 4
  if(j.lt.i)goto 5
  itemp=indx(i)
  indx(i)=indx(j)
  indx(j)=itemp
  goto 3
5 indx(l)=indx(j)
  indx(j)=indxt
  jstack=jstack+2
  if(jstack.gt.nstack)pause 'nstack too small in indexx'
  if(ir-i+1.ge.j-l)then
   istack(jstack)=ir
   istack(jstack-1)=i
   ir=j-1
  else
   istack(jstack)=j-1
   istack(jstack-1)=l
   l=i
  endif
 endif
 goto 1
end subroutine
! **********************************************************************
! *** subroutine polyfit                                             ***
! *** returns the coefficients of a polynomial function that fits    ***
! *** a dataset                                                      ***
! **********************************************************************
subroutine polyfit(npt,npa,x,y,w,p)
 implicit none
 integer :: npt,npa
 real(8), intent(in) :: x(npt),y(npt),w(npt)
 real(8), intent(out) :: p(npa)
 real(8) :: a(npa,npa),z(npa,npa),v(npa)
 integer :: ipa,jpa,ipt

! *** construct the matrix Ajk = Df(xi)Dpj wii Df(xi)Dpk
 do ipa=1,npa
  do jpa=1,npa
   a(ipa,jpa)=0.d0
   do ipt=1,npt
    a(ipa,jpa)=a(ipa,jpa)+(x(ipt)**real(ipa-1,8))*w(ipt)*(x(ipt)**real(jpa-1,8))
   enddo
  enddo
 enddo
! *** construct the vector Vk = Df(xi)Dpk wii yi
 do ipa=1,npa
  v(ipa)=0.d0
  do ipt=1,npt
   v(ipa)=v(ipa)+(x(ipt)**real(ipa-1,8))*w(ipt)*y(ipt)
  enddo
 enddo
! *** invert M matrix
 z=psdinvmat(a,npa,npa,npa,npa)
! *** determine optimal parameter
 do ipa=1,npa
  p(ipa)=0.d0
  do jpa=1,npa
   p(ipa)=p(ipa)+z(ipa,jpa)*v(jpa)
  enddo
 enddo
end subroutine
! **********************************************************************
! *** subroutine weighted_linear_regression                          ***
! *** returns the coefficients of an affine function that fits       ***
! *** a dataset and the corresponding uncertainties                  ***
! **********************************************************************
subroutine weighted_linear_regression(n,x,y,w,slope,intercept,slope_err,intercept_err,status)
 implicit none
 integer, intent(in) :: n
 real(8), dimension(n), intent(in) :: x,y
 real(8), dimension(n), optional, intent(in) :: w
 real(8), intent(out) :: slope
 real(8), intent(out) :: intercept
 real(8), intent(out), optional :: slope_err
 real(8), intent(out), optional :: intercept_err
 integer, intent(out), optional :: status
 real(8) :: sum_w,sum_wx,sum_wy,sum_wxx,sum_wxy
 real(8) :: delta
 integer :: i

 if (present(status)) then
  status = 0
 endif
 if (n < 2) then
  if (present(status)) then
   status = -1
  endif
  return
 endif
 sum_w=0.0d0
 sum_wx=0.0d0
 sum_wy=0.0d0
 sum_wxx=0.0d0
 sum_wxy=0.0d0
 if (present(w)) then
  do i = 1, n
   sum_w=sum_w+w(i)
   sum_wx=sum_wx+w(i)*x(i)
   sum_wy=sum_wy+w(i)*y(i)
   sum_wxx=sum_wxx+w(i)*x(i)**2
   sum_wxy=sum_wxy+w(i)*x(i)*y(i)
  enddo
 else
  do i = 1, n
   sum_w=sum_w+1.d0
   sum_wx=sum_wx+x(i)
   sum_wy=sum_wy+y(i)
   sum_wxx=sum_wxx+x(i)**2
   sum_wxy=sum_wxy+x(i)*y(i)
  enddo
 endif
 delta=sum_w*sum_wxx-sum_wx**2
 if (abs(delta) < 1.0d-12) then
  if (present(status)) then
   status = -1
  endif
  slope=0.0d0
  intercept=0.0d0
  return
 end if
 slope=(sum_w*sum_wxy-sum_wx*sum_wy)/delta
 intercept=(sum_wxx*sum_wy-sum_wx*sum_wxy)/delta
 if (present(slope_err)) then
  slope_err=sqrt(sum_w/delta)
 end if
 if (present(intercept_err)) then
  intercept_err=sqrt(sum_wxx/delta)
 end if
end subroutine weighted_linear_regression
end module
