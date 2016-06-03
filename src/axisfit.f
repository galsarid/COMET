!****************************
      FUNCTION AXISFIT(TIME0)
!****************************

!       Function description:
!       ---------------------
!       Power-law fit, for A and B, within an interval [TINIT,TEND].
!       This should follow a fit to the orbital data - semi-major axis

!       Function mapping:
!       -----------------
!       Called in: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
! Fit to Tempel 1 orbital data:
      DATA TINIT/1.0470D10/, TEND/1.2141D10/
      DATA A/12.221D0/, B/-0.0584D0/

      TIME=TIME0+TINIT
      IF(TIME.GE.TEND)THEN
        WRITE(2,*)TIME,TEND
        STOP
      ENDIF
      AXISFIT=A*TIME**B
      END

