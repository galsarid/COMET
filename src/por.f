!*********************************************************************
      FUNCTION POR(ROA,ROB,ROC,RODUST,RHOICE,RHODUS,MMAX,ROSGI,RHOSOL)
!*********************************************************************

!       Function description:
!       ---------------------
!       Calculate value of porosity.

!       Function mapping:
!       -----------------
!       Called in: DIVIS, SMERGE, GRIDEF, COMPOS, INPUT, ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION  RHOSOL(*),ROSGI(*)
      
      SUMX=ZERO
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        SUMX=SUMX+ROSGI(M)/RHOSOL(M)
       END DO
      END IF
      POR=ONE-((ROA+ROB+ROC)/RHOICE+RODUST/RHODUS+SUMX)
      IF(POR.GT.0.99D0)THEN
       POR=0.99D0
!       WRITE(2,*) 'Warning: porosity out of range', POR
      END IF
      IF(POR.LT.1.D-2)THEN
       POR=1.D-2
!       WRITE(2,*) 'Warning: porosity out of range', POR
      END IF

      RETURN
      END

