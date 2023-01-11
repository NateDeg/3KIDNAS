cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains modified versions of the amoeba and
c       amotry routines from Press et al.  They implement the
c       Nelder and Mead method of minimizing a multi-dimensional function
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DownhillSimplexMod



      implicit none

      contains

ccccccc
      subroutine amoeba(p,y,mp,np,ndim,ftol,funk,iter,PID
     &          ,NoConvergenceFlag)
      implicit none
      integer iter,mp,ndim,np
      integer NMAX,ITMAX
      integer PID

      real ftol,p(mp,np),y(mp),TINY
      parameter(NMAX=20,ITMAX=5000,TINY=1.e-10)
      external funk

      integer i,ihi,ilo,inhi,j,m,n
      real rtol,sum,swap,ysave,ytry
      real,ALLOCATABLE :: psum(:)
      logical,intent(INOUT):: NoConvergenceFlag
c
      ALLOCATE(psum(ndim))

      NoConvergenceFlag=.False.
c      print*, "Starting downhill simplex",PID,y
      iter=0
 1    do n=1,ndim
        sum=0.
        do m=1,ndim+1
            sum=sum+p(m,n)
        enddo
        psum(n)=sum
      enddo

2     ilo=1
      if(y(1) .gt. y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do i=1,ndim+1
        if(y(i) .le. y(ilo)) ilo=i
        if(y(i) .gt. y(ihi)) then
            inhi=ihi
            ihi=i
        elseif(y(i) .gt. y(inhi)) then
            if(i .ne. ihi) inhi=i
        endif
      enddo

      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
      print*, "Current tolerance", iter,rtol,y(ihi),y(ilo)
      if(rtol .lt. ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do n=1,ndim
            swap=p(1,n)
            p(1,n)=p(ilo,n)
            p(ilo,n)=swap
        enddo
        return
      endif
c
      if(iter .ge. ITMAX) then
        NoConvergenceFlag=.True.
        return
      endif

      iter=iter+2
      call amotry(ytry, p,y,psum,mp,np,ndim,funk,ihi,-1.0)
      if(ytry .le. y(ilo)) then
        call amotry(ytry, p,y,psum,mp,np,ndim,funk,ihi,2.0)
      elseif(ytry .ge. y(inhi)) then
        ysave=y(ihi)
        call amotry(ytry, p,y,psum,mp,np,ndim,funk,ihi,0.5)
        if(ytry .ge. ysave) then
            do i=1,ndim+1
                if(i .ne. ilo) then
                    do j=1,ndim
                        psum(j)=0.5*(p(i,j)+p(ilo,j))
                        p(i,j)=psum(j)
                    enddo
                    call funk(psum,ytry)
                endif
            enddo
            iter=iter+ndim
            goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2


c
      DEALLOCATE(psum)
      return
      end subroutine
ccccccccccc


      subroutine amotry(ynew,p,y,psum,mp,np,ndim,funk,ihi,fac)
      implicit none
      integer ihi,mp,ndim,np
      real ynew,fac,p(mp,np),psum(np),y(mp)
      external funk

      integer j
      real fac1,fac2,ytry
      real,ALLOCATABLE :: ptry(:)
c
      ALLOCATE(ptry(1:ndim))
c      print*, "in amotry"
c
      fac1=(1.-fac)/real(ndim)
      fac2=fac1-fac
      do j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
      enddo
c      print*, p(ihi,:)
c      print*, ptry
c      print*, shape(ptry)

      call funk(ptry,ytry)
      if(ytry .lt. y(ihi)) then
        y(ihi)=ytry
        do j=1,ndim
            psum(j)=psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j)=ptry(j)
        enddo
      endif
      ynew=ytry

      DEALLOCATE(ptry)
      return
      end subroutine

      end module
