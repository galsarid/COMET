!*****************************
      FUNCTION ECCENFIT(TIME0)
!*****************************

!       Function description:
!       ---------------------
!       Power-law fit, for A and B, within an interval [TINIT,TEND].
!       This should follow fit to the orbital data - eccentricity.

!       Function mapping:
!       -----------------
!       Called in: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
! Fit to Tempel 1 orbital data:
      DATA TINIT/1.0470D10/, TEND/1.2141D10/
      DATA A/0.011D0/, B/0.1653D0/

      TIME=TIME0+TINIT
      IF(TIME.GE.TEND)THEN
        WRITE(2,*)TIME,TEND
        STOP
      ENDIF
      ECCENFIT=A*TIME**B
      END

