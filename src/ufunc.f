!******************************************************
      FUNCTION UFUNC(T,u,xi,usol,alfa,beta,eta,H,Tm,cw)
!******************************************************

!       Function description:
!       ---------------------
!       Heat function for root finding (RTBIS).

!       Function mapping:
!       -----------------
!       Called in: RTBIS.
 
      REAL*8 ufunc,T,u,xi,usol,alfa,beta,eta,H,Tm,cw,xm

      xm=1.D0/(1.D0+exp(eta*(1.D0-T/Tm)))
      ufunc=usol*T+xi*(alfa*T*T+beta*T+
     1      xm*(H+(cw-beta)*T-alfa*T*T))-u
      END

