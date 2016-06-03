!***********************
      SUBROUTINE  IMPACT
!***********************

!       Subroutine description:
!       -----------------------
!       Calculate the effect of an impact onto the nucleus.

!       At the current scheme, the impact is simulated by removing a 
!       layer of given thickness THICK (cm) and adding a surface energy 
!       source in the form of a Gaussian, starting at a given distance 
!       DDI (AU) on the inbound leg of a given orbit NDI, for a given 
!       timescale DELTDI(sec). The onset of the impact is marked by IIMPAC>0.
!       The total energy (kinetic) of the impact is given by ENAREA (erg). 
!       The diameter of the resulting crater is CRDIAM (cm). The energy 
!       partition factor (between mechanical and heating) is ENPART.


!       Subroutine mapping:
!       -------------------
!       Calls: DIVIS.
!       Called by: ENDSTP, FLUX.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      DIMENSION  DGTOT(LG)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      SAVE
      DATA KIN/0/, KPR/0/


      IF(KPR.EQ.1)THEN
       SAREA=FACR(IMAXP1)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         DVTOT=-(VAFL(IMAXP1)-VAPFL-VAPFLA)/AMICE/SAREA
         DDOT=-XMFL(IMAXP1)/SAREA
         DGTOT(M)=-(GAFL(IMAXP1,M)-VGFL(M)-VGFLA(M))/AMG(M)/SAREA
        END DO
       END IF
       WRITE(17,'(1X,I5,1PE15.7,1PE13.5,1P8E10.2,1PE10.3,1PE12.5)')
     1             NORB,TIME-TIMPACT,RHO(IMAXP1),T(IMAXP1),DVTOT,DDOT,
     2             DGTOT(1),DGTOT(2),DGTOT(3),DGTOT(4),
     3             DGTOT(5),FL(IMAXP1)/SAREA,R(IMAXP1)
      END IF
      RETURN

      IF(KIN.EQ.0)THEN
       OPEN(17,file='deep_impact')
       KIN=1
       DDI=DDI*ONEAU
       FIDI=ACOS((ONE-DDI/AXIS)/ECCEN)
       TIMPACT=-PSICON*(FIDI-ECCEN*SIN(FIDI)-PI)+PSICON*(NDI-1)*TWO*PI
!       DELDDI=ABS(ECCEN*SQRT(GM*AXIS)*SIN(FIDI)/DDI*6.D0*DELTDI)
       DELDDI=ABS(ECCEN*SQRT(GM*AXIS)*SIN(FIDI)/DDI*24.D0*DELTDI)
       ENTOT=KENRAT*ENTOTDI
       ENAREA=ENTOT/(PI*CRDIAM**2)
       FSUR0=ENAREA/DELTDI
       IDI=0
       DO I=1,IMAXP1-2
        IF(THICK.GT.(R(IMAXP1)-R(IMAXP1-I))) IDI=I
       END DO
       WRITE(2,*)'IMPACT PARA',FSUR0,FIDI,DDI,DELDDI,TIMPACT,DELTDI,IDI
       IF(TIMPACT.LT.ZERO) STOP
      END IF

!  Check "impact mode'

      IF(NORB.EQ.NDI.AND.KPR.GE.0)THEN
       IF(DIST.LE.(DDI+DELDDI))THEN
        KPR=1
        IF(KIN.EQ.1)THEN
         KIN=2
         ITMODE0=ITMODE
         ITMODE=1
         EPSX20=EPSX2
         EPSX2=EPSX2*1.D2
         IMAXP1=IMAXP1-IDI
         MODIV=1
         IP=IMAXP1-2
         CALL DIVIS(MODIV,IP)
         IP=IMAXP1-1
         CALL DIVIS(MODIV,IP)
         DO NDIVIS=1,5
          IP=IMAXP1
          CALL DIVIS(MODIV,IP)
         END DO
        END IF
        FIMPACT=ENPART*FSUR0*EXP(-PI*((TIME-TIMPACT)/DELTDI)**2)
       END IF
       IF(DIST.LE.(DDI-DELDDI))THEN
        KPR=-1
        ITMODE=ITMODE0
        EPSX2=EPSX20
       END IF
      END IF

      RETURN
      END
