!******************************************
      SUBROUTINE TABLE(KDI,IL,IP,IJ,PRIMAT)
!******************************************

!       Subroutine description:
!       -----------------------
!

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: PRINT.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      DIMENSION PRIMAT(IMX,11)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
  
      DO I=IL,IP,IJ
       ITAB=I
       WRITE(2,10) I,(PRIMAT(I,K),K=1,KDI)
       DO KK=1,KDI
        PRIMAT(I,KK)=ZERO
       END DO
      END DO
      IF(ITAB.LT.IP)THEN 
       WRITE(2,10) IP,(PRIMAT(IP,K),K=1,KDI)
      END IF
  10  FORMAT(1X,I5,1P11E11.3)
      
      RETURN
      END

