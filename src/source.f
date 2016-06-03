!**********************
      SUBROUTINE SOURCE
!**********************

!       Subroutine description:
!       -----------------------
!       Calculate vapor and gas gains or losses, QV(I) and QG(I),
!       due to sublimation and crystallization, and energy gain, S(I), 
!       due to crystallization and radioactive decay.

!       Subroutine mapping:
!       -------------------
!       Calls: SUBLI. 
!       Called by: Main.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0,FOUR=4.D0)
      DIMENSION ROCINT(IMX)
      SAVE


      TAC=TZERO
      SVCORR=ONE

!  Calculate correction for surface-to-volume ratio, assuming a 
!  power law distribution of pore sizes, with a given exponent and 
!  pore size range (ALFPS & XPS, see COMIN file).

      IF(ICORR.GT.0)THEN
       DO I=1,IMAXP1
        IF(NMODEL.LT.0) ALFP(I)=ALFPS
       END DO
       SVCORR=(ALFPS-THRE)**2/((ALFPS-TWO)*(ALFPS-FOUR))
     1         *(XPS**(ALFPS-TWO)-ONE)*(XPS**(ALFPS-FOUR)-ONE)
     2         /(XPS**(ALFPS-THRE)-ONE)**2
       COR1=(ALFPS-FOUR)/(ALFPS-THRE)*(XPS**(ALFPS-THRE)-ONE)
     1       /(XPS**(ALFPS-FOUR)-ONE)
       PSIZE=TET(IMAXP1-1)
       IF(MOD(NTIME,NPRINT).EQ.0)THEN
        WRITE(2,*) 'SVR correction: ',
     1             'ALFPS,XPS,PSIZE*COR,PSIZE*COR/XPS,SVCORR = ',
     2              ALFPS,XPS,PSIZE*COR1,PSIZE*COR1/XPS,SVCORR
        WRITE(17,*) 'SVR correction: ',
     1              'ALFPS,XPS,PSIZE*COR,PSIZE*COR/XPS,SVCORR = ',
     2               ALFPS,XPS,PSIZE*COR1,PSIZE*COR1/XPS,SVCORR
       END IF
      END IF

      STOT=ZERO

      DO I=2,IMAXP1

       TI=T(I)
       T0I=T0(I)
       RHO0I=RHO0(I)
       S(I)=ZERO
       QV(I)=ZERO
       DQV1(I)=ZERO
       DQV2(I)=ZERO
       DQVDT(I)=ZERO
       ZGAS=ZERO
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         ZGAS=ZGAS+ZGI(I,M)
         QVG(I,M)=ZERO
         DQVG1(I,M)=ZERO
         DQVG2(I,M)=ZERO
         DQVGDT(I,M)=ZERO
        END DO
       END IF

!  Phase transition: Amorphous to Crystalline ice

       IF(T0I.GT.TAC.AND.ROA0(I).GT.ZERO.AND.TI.GT.TMIN)THEN
        RATAC=ATAU*EXP(-EXPTAU/TI)
        DTTAU=DTIME*RATAC
        ROAINT(I)=ROA0(I)/(ONE+DTTAU)
        IF(ROAINT(I).LT.1.D-15) ROAINT(I)=ZERO
        ROTAU=ROAINT(I)*RATAC
       ELSE
        ROAINT(I)=ROA0(I)
        ROTAU=ZERO
        DTTAU=ZERO
       END IF
       ROCINT(I)=(ONE+(ROC0(I)+ROA0(I))-ROAINT(I))-ONE
       S(I)=ROTAU*(ONE-ZGAS)*HAMC
       DSDT(I)=S(I)*EXPTAU/TI**2/(ONE+DTTAU)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         QG(I,M)=ROTAU*ZGI(I,M)
         DQGDT(I,M)=QG(I,M)*EXPTAU/TI**2/(ONE+DTTAU)
        END DO
       END IF

!  Radioactivity

       IF(LMAX.GT.0)THEN
        DO L=1,LMAX
         XMRADIL=XMRAD0(I,L)/(ONE+DTIME/TAUR(L))
         IF(XMRADIL.LT.1.D-15) XMRADIL=ZERO
         S(I)=S(I)+XMRAD0(I,L)*HRAD(L)/(TAUR(L)+DTIME)/DV(I)
         XMRAD(I,L)=XMRADIL
        END DO
       END IF

       STOT=STOT+S(I)*DV(I)

!  Sublimation of volatiles

       SURVOL=ZERO
       IF(ICORR.GT.0)THEN
        ALFPS=ALFP(I)
        SVCORR=(ALFPS-THRE)**2/((ALFPS-TWO)*(ALFPS-FOUR))
     1         *(XPS**(ALFPS-TWO)-ONE)*(XPS**(ALFPS-FOUR)-ONE)
     2         /(XPS**(ALFPS-THRE)-ONE)**2
       END IF
       TETI=TET(I)
       IF(TETI.GT.ZERO)THEN
        PORI=PO(I)
        IF(PORI.LT.0.6D0)THEN
         SURVOL=TWO*PORI/TETI*SVCORR
        ELSE
         SURVOL=THRE*(ONE-PORI)/TETI*SVCORR
        END IF
       END IF
       IF(I.EQ.IMAXP1)THEN
        SURVB=ZERO
        IF(TETAB.GT.ZERO)THEN
         IF(TETI.GT.ZERO)THEN
          SURVB=SURVOL
         ELSE
          PORI=PO(I)
          IF(PORI.LT.0.6D0)THEN
           SURVB=TWO*PORI/TETAB*SVCORR
          ELSE
           SURVB=THRE*(ONE-PORI)/TETAB*SVCORR
          END IF
         END IF
        END IF
        PQV=SURVB*ROC0(I)/RHO0I
       END IF
       IF(TI.GT.TSUBV)THEN
        QMAXV=(ROCINT(I)+ROB0(I))/DTIME
        CALL SUBLI(PZERO,EXPT,TI,T0I,ROV(I),QMAXV,
     1             SURVOL,AIB,PI,PO(I),QV(I),DQV1(I),
     2             DQV2(I),DQVDT(I))
       END IF
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         IF(TI.GT.TSUB(M))THEN
          QMAXG=ROSG0(I,M)/DTIME
          CALL SUBLI(PZEROG(M),EXPTG(M),TI,T0I,ROG(I,M),QMAXG,
     1               SURVOL,AGB(M),PI,PO(I),QVG(I,M),DQVG1(I,M),
     2               DQVG2(I,M),DQVGDT(I,M))
         END IF
        END DO
       END IF

      END DO

      RETURN
      END
