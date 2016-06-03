!***********************
      SUBROUTINE PROFILE
!***********************

!       Subroutine description:
!       -----------------------
!       Prepare and print data to COMPROF - profile output file.

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      DIMENSION PRIMAT(25)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      SAVE
      
      RIP=R(IMAXP1)
      YEAR=3.15576D7

      DO I=2,IMAXP1
       PRIMAT(1)=LOG10(TIME/YEAR)
       PRIMAT(2)=LOG10(RIP-R(I-1))-TWO
       PRIMAT(3)=T(I)
       PRIMAT(4)=RHO(I)
       PRIMAT(5)=PO(I)
       PRIMAT(6)=ROA(I)/RHO(I)
       PRIMAT(7)=ROC(I)/RHO(I)
       PRIMAT(8)=ROB(I)/RHO(I)
       PRIMAT(9)=P(I)
       PRIMAT(10)=XM0(I)
       PRIMAT(11)=COND(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         PRIMAT(11+M)=ROSG(I,M)/RHO(I)
        END DO
       END IF

       WRITE(4,'(5X,1P12E12.4)')(PRIMAT(K),K=1,11+MMAX)
      END DO

      RETURN
      END

