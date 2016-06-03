!*********************
      SUBROUTINE MODEL
!*********************

!       Subroutine description:
!       -----------------------
!       Write data into COMODO - model output file.

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: Main, ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      SAVE

      WRITE(8,*)NTIME
      WRITE(8,*)IMAXP1,JMAX,MMAX,LMAX,ICLAT
      WRITE(8,*)TIME
      DFIMOD=TWO*PI/DFI
      WRITE(8,*)FI,DFIMOD
      WRITE(8,*)AXIS,ECCEN
      DO I=2,IMAXP1
       XA(I)=ROA0(I)/RHO0(I)
       XB(I)=ROB0(I)/RHO0(I)
       XC(I)=ROC0(I)/RHO0(I)
       XD(I)=ROD0(I)/RHO0(I)
       XV(I)=ROV0(I)/RHO0(I)
       DO M=1,LG
        XSG(I,M)=ROSG0(I,M)/RHO0(I)
        XG(I,M)=ROG0(I,M)/RHO0(I)
       END DO
       WRITE(8,*)I,DM0(I),RHO0(I),T0(I),TET(I),XA(I),XB(I),XC(I),
     1          (XSG(I,M),M=1,LG),XD(I),(XG(I,M),M=1,LG),XV(I),
     2          (ZGI(I,M),M=1,LG),PCR(I),PSOL(I),ALFP(I)
       IF(LMAX.GT.0)THEN
        WRITE(8,*)I,(XRAD(I,L),L=1,LMAX)
       END IF
      END DO

      WRITE(2,890)NTIME,TIME,FI
      WRITE(17,890)NTIME,TIME,FI
 890  FORMAT(1X,'MODEL NUMBER',I10,' SAVED, AT TIME=',1PE15.5,
     1       6H FI  =,1PE15.5)

      RETURN
      END

