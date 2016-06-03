!*************************************************
      SUBROUTINE SUBLI(AP,BP,TI,T0I,RO,QMAX,
     1                 STOV,AIB,PI,PORI,QVI,DQV1I,
     2                 DQV2I,DQVDTI)
!*************************************************

!       Subroutine description:
!       -----------------------
!

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: SOURCE.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)

      T=TI
      QV=STOV*SQRT(AIB/(TWO*PI*T))
      PVI=AP*EXP(BP/T)
      P=RO*T/AIB/PORI
      QVI=QV*(PVI-P)
         DQV1I=QV*PVI
         DQV2I=QV*T/AIB/PORI
         IF(QVI.GT.QMAX)THEN
         QVI=QMAX
         DQV1I=QMAX
         DQV2I=ZERO
         QV=ZERO
         ENDIF
      DQVDTI=-QV/T*((0.5D0+BP/T)*PVI+0.5D0*P)
      RETURN
      END

