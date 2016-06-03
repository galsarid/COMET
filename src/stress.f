!****************************************************
      SUBROUTINE STRESS(IP1,PI43,PS,RS,SRR,STT,NTIME)
!****************************************************

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
      DIMENSION  P(IMX),R(IMX),V(IMX),DV(IMX),SUMI(IMX),PGAL(IMX)
      DIMENSION  SRR(IMX),STT(IMX),PS(IMX),RS(IMX)
      SAVE
      DATA  SIG/0.33D0/,PMIN/1.D3/

      IP=IP1
      DO I=2,IP
       P(I)=PS(I)
       R(I)=RS(I)
      END DO
      DO I=2,IP
       SRR(I)=ZERO
       STT(I)=ZERO
       J=I
       IF(P(J).GT.PMIN)  GO  TO  1
      END DO
      WRITE(2,'(5X,17HNO STRESSES FOUND)')
      WRITE(17,*)NTIME
      WRITE(17,'(5X,17HNO STRESSES FOUND)')
      RETURN
   1  JM=J-1
      P(JM)=ZERO
      P(IP+1)=ZERO
      V(JM)=PI43*R(JM)**3
      DV(JM)=ZERO
      SUMI(JM)=ZERO
      DO I=J,IP
       V(I)=PI43*R(I)**3
       DV(I)=V(I)-V(I-1)
       SUMI(I)=ZERO
      END DO
      DV(IP+1)=ZERO
      SUMX=ZERO
      DO I=JM,IP
       PGAL(I)=(DV(I+1)*P(I)+DV(I)*P(I+1))/(DV(I)+DV(I+1))
       TERM=P(I)*DV(I)
       SUMX=SUMX+TERM
       SUMI(I)=SUMX
      END DO
      C1=SIG/(ONE-SIG)
      C2=TWO*(ONE-TWO*SIG)/(THRE*(ONE-SIG))
      IF(J.EQ.2)THEN
       WRITE(2,12)
       WRITE(17,*)NTIME
       WRITE(17,12)
  12   FORMAT(5X,'NON-VANISHING PRESSURES DOWN TO THE CENTER')
      ELSE
       DO I=JM,IP
        SRR(I)=PGAL(I)-C2/V(I)*(SUMI(I)-(V(I)-V(JM))
     1         /(V(IP)-V(JM))*SUM)
        STT(I)=C1*PGAL(I)+C2/V(I)*(SUMI(I)+(TWO*V(I)+V(JM))
     1         /(V(IP)-V(JM))*SUM)
       END DO
      END IF

      RETURN
      END

