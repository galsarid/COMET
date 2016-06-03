!*******************************************************
      FUNCTION RTBIS(x1,u,xi,usol,alfa,beta,eta,H,Tm,cw)
!*******************************************************

!       Function description:
!       ---------------------
!       Perform root finding by the bisection method, to the function defined by UFUNC.

!       Function mapping:
!       -----------------
!       Called in: DIVIS, SMERGE, GRIDEF.

      INTEGER JMAX,J
      REAL*8 rtbis,x1,u,ufunc,xi,usol,alfa,beta,eta,H,Tm,cw
      REAL*8 x2,dx,f,fmid,xmid,xacc,delx
      PARAMETER (JMAX=40,delx=1.D1,xacc=1.D-6)

      f=UFUNC(x1,u,xi,usol,alfa,beta,eta,H,Tm,cw)
      do j=1,JMAX
       if(f.lt.0.D0)then
        x2=x1+delx
        fmid=UFUNC(x2,u,xi,usol,alfa,beta,eta,H,Tm,cw)
         if(fmid.ge.0.D0)then
          go to 1
         else
          x1=x2
          f=fmid
         end if
       else
        x2=max(0.D0,x1-delx)
        fmid=UFUNC(x2,u,xi,usol,alfa,beta,eta,H,Tm,cw)
        if(fmid.lt.0.D0)then
         go to 1
        else
         x1=x2
         f=fmid
        end if
       end if
      end do

    1 continue
      if(f.lt.0.D0)then
       rtbis=x1
       dx=x2-x1
      else
       rtbis=x2
       dx=x1-x2
      end if
      do j=1,JMAX
       dx=dx*5.D-1
       xmid=rtbis+dx
       fmid=UFUNC(xmid,u,xi,usol,alfa,beta,eta,H,Tm,cw)
       if(fmid.le.0.D0)rtbis=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.D0) return
      end do

      pause 'too many bisections in rtbis'

      END

