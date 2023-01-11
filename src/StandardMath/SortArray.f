      module sort
      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     index
c
c     This sorts an array or set of arrays by size
c
c     IT's from numerical recipies
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine indexx(n, arr, indx)
c
      implicit none
c
      integer n, indx(n), M, NSTACK
c
      real arr(n)
      parameter(M=7, NSTACK=100)
c
      integer i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK)
      real a
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      print*, 'indexx',size(arr), n
      do 11 j=1, n
         indx(j)=j
 11   enddo

      jstack=0
      l=1
      ir=n
 1    if(ir-l .lt. M) then
         do 13  j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12  i=j-1,1,-1
               if(arr(indx(i)) .le. a) goto 2
               indx(i+1)=indx(i)
 12         enddo 
            i=l-1
 2          indx(i+1)=indxt
 13      enddo 
         if(jstack .eq. 0) return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l)) .gt. arr(indx(ir))) then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)) .gt. arr(indx(ir))) then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)) .gt. arr(indx(l+1))) then
            itemp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l+1)
         a=arr(indxt)
 3       continue
         i=i+1
         if(arr(indx(i)) .lt. a) goto 3
 4       continue
         j=j-1
         if(arr(indx(j)) .gt. a) goto 4
         if (j .lt. i) goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
 5       indx(l+1)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack .gt. NSTACK) then
            print*, 'NSTACK too small in indexx'
            stop
         endif
         if(ir-i+1 .ge. j-1) then
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     index
c
c     This sorts an array or set of arrays by size
c
c     IT's from numerical recipies
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine indexxD(n, arr, indx)
c
      implicit none
c
      integer n, indx(n), M, NSTACK
c
      double precision arr(n)
      parameter(M=7, NSTACK=100)
c
      integer i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK)
      double precision a
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      print*, 'indexx', arr, n
      do 11 j=1, n
         indx(j)=j
 11   enddo

      jstack=0
      l=1
      ir=n
 1    if(ir-l .lt. M) then
         do 13  j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12  i=j-1,1,-1
               if(arr(indx(i)) .le. a) goto 2
               indx(i+1)=indx(i)
 12         enddo 
            i=l-1
 2          indx(i+1)=indxt
 13      enddo 
         if(jstack .eq. 0) return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l)) .gt. arr(indx(ir))) then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)) .gt. arr(indx(ir))) then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)) .gt. arr(indx(l+1))) then
            itemp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l+1)
         a=arr(indxt)
 3       continue
         i=i+1
         if(arr(indx(i)) .lt. a) goto 3
 4       continue
         j=j-1
         if(arr(indx(j)) .gt. a) goto 4
         if (j .lt. i) goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
 5       indx(l+1)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack .gt. NSTACK) then
            print*, 'NSTACK too small in indexx'
            stop
         endif
         if(ir-i+1 .ge. j-1) then
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end module
