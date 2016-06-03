!**********************************************
      SUBROUTINE TSTRES(IP1,PI43,TS,RS,TSR,TST)
!**********************************************

!       Subroutine description:
!       -----------------------
!

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: PRINT.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION T(IMX),R(IMX),V(IMX),DV(IMX),SUMI(IMX),TGAL(IMX)
      DIMENSION TS(IMX),RS(IMX),TSR(IMX),TST(IMX)
      SAVE
      DATA  SIG/0.33D0/,E/9.D10/,ALPHA/3.8D-5/

      EA3=E*ALPHA/THRE
      IP=IP1
      V(1)=ZERO
      DO 100 I=2,IP
      T(I)=TS(I)
      R(I)=RS(I)
      TSR(I)=ZERO
      TST(I)=ZERO
      V(I)=PI43*R(I)**3
 100  DV(I)=V(I)-V(I-1)
      T(IP+1)=T(IP)
      DV(IP+1)=ZERO
      SUMX=ZERO
      DO 200 I=2,IP
      TGAL(I)=(DV(I+1)*T(I)+DV(I)*T(I+1))/(DV(I)+DV(I+1))
      SUMX=SUMX+T(I)*DV(I)
 200  SUMI(I)=SUMX
      TR=T(IP)
      SOV=SUMX/V(IP)
      C1=ONE/(ONE-TWO*SIG)
      C2=ONE/(THRE*(ONE-SIG))
      C3=SIG/(ONE-SIG)
      DO 300 I=2,IP
      TSR(I)=EA3*(C1*(TGAL(I)-TR)-TWO*C2*(SUMI(I)/V(I)-SOV))
 300  TST(I)=EA3*(C1*(C3*TGAL(I)-TR)+C2*(SUMI(I)/V(I)+TWO*SOV))

      RETURN
      END

