!**********************
      SUBROUTINE COMPOS
!**********************

!       Subroutine description:
!       -----------------------
!       Calculate composition mass fractions for given T and RHO

!       Subroutine mapping:
!       -------------------
!       Calls: SGTSL.
!       Called by: Main.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION A(IMX),B(IMX),C(IMX),D(IMX), SUBB(IMX)
      DIMENSION SUMG(LG),SUMSG(LG),SUMPG(LG),SUMPSG(LG)
      DIMENSION FACNG(LG),FACNSG(LG)
      SAVE

      DELMT=ZERO
      SUMA=ZERO
      SUMB=ZERO
      SUMC=ZERO
      SUMV=ZERO
      SUMD=ZERO
      SUMPA=ZERO
      SUMPB=ZERO
      SUMPC=ZERO
      SUMPV=ZERO
      SUMPD=ZERO
      FACNA=ONE
      FACNB=ONE
      FACNC=ONE
      FACNV=ONE
      FACND=ONE
      DO I=1,IMAXP1
       SUBB(I)=ZERO
      END DO
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        SUMG(M)=ZERO
        SUMSG(M)=ZERO
        SUMPG(M)=ZERO
        SUMPSG(M)=ZERO
        FACNG(M)=ONE
        FACNSG(M)=ONE
       END DO
      END IF

! Calculate water vapor density by linear matrix inversion

      DO I=2,IMAXP1
       IM=I-1
       DVI=DV(I)
       DTOV=DTIME/DVI
       A(IM)=-DTOV*VFR(I)
       C(IM)=-DTOV*VFL(IM)
       B(IM)=ONE+DTOV*(VFL(I)+VFR(IM))+DQV2(I)*DTIME
       D(IM)=ROV0(I)+DQV1(I)*DTIME
      END DO
      IF(IB.GT.0)THEN
       I=IMAXP1
       B(I-1)=ONE+DTOV*VFR(I-1)/SQV
       D(I-1)=D(I-1)+DTOV*VAFL(I)
      ENDIF
      CALL SGTSL(IMAXP1-1,C,B,A,D,INFO)
      IF(INFO.NE.0)THEN
       WRITE(2,*)'ERROR: Matrix inversion for H2O vapor (COMPOS). ',
     1            'INFO element of the diagonal=0 (INFO=',INFO,'). ',
     2            'Exiting.'

       WRITE(17,*)NTIME
       WRITE(17,*)'ERROR: Matrix inversion for H2O vapor (COMPOS). ',
     1            'INFO element of the diagonal=0. (INFO=',INFO,') ',
     2            'Exiting.'
       STOP
      ENDIF

! Calculate water ice densities
      DO I=2,IMAXP1

       ROV(I)=D(I-1)
       ZGAS=ZERO
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         ZGAS=ZGAS+ZGI(I,M)
        END DO
       END IF

       IF((ROC0(I)+ROB0(I)).GT.ZERO)THEN
        CORRAC=ROC0(I)/(ROC0(I)+ROB0(I))
       ELSE
        CORRAC=ONE
       END IF

       ROVI=ROV(I)
       QVI=DQV1(I)-DQV2(I)*ROVI
       IF(QVI.LT.ZERO)THEN
! Condensation to form crystalline ice
        ROCI=ROC0(I)+(ONE-ZGAS)*(ROA0(I)-ROAINT(I))-QVI*DTIME
        ROBI=ROB0(I)
       ELSE
! Sublimation of crystalline and remnant amorphous ice (or clathrate)
        ROCI=ROC0(I)+(ONE-ZGAS)*(ROA0(I)-ROAINT(I))-QVI*DTIME*CORRAC
        ROBI=ROB0(I)-QVI*DTIME*(ONE-CORRAC)
        SUBB(I)=QVI*(ONE-CORRAC)
       END IF

       IF(I.EQ.IMAXP1)THEN
        ROCI=ROCI-VAPFL*DTIME/DV(I)
        ROBI=ROBI-VAPFLA*DTIME/DV(I)/(ONE-ZGAS)
       END IF

       SUMC=SUMC+ROCI*DV(I)
       SUMB=SUMB+ROBI*DV(I)
       SUMA=SUMA+ROAINT(I)*DV(I)
       SUMV=SUMV+ROVI*DV(I)
       ROC(I)=MAX(ZERO,ROCI)
       ROB(I)=MAX(ZERO,ROBI)
       ROA(I)=MAX(ZERO,ROAINT(I))
       ROV(I)=MAX(ZERO,ROVI)
       SUMPC=SUMPC+ROC(I)*DV(I)
       SUMPB=SUMPB+ROB(I)*DV(I)
       SUMPA=SUMPA+ROA(I)*DV(I)
       SUMPV=SUMPV+ROV(I)*DV(I)

      END DO

      IF(MMAX.GT.0)THEN
       DO M=1,MMAX

! Calculate volatile densities by linear matrices inversion
        DO I=2,IMAXP1
         IM=I-1
         DTOV=DTIME/DV(I)
         A(IM)=-DTOV*GFR(I,M)
         C(IM)=-DTOV*GFL(IM,M)
         B(IM)=ONE+DTOV*(GFL(I,M)+GFR(IM,M))+DQVG2(I,M)*DTIME
         D(IM)=ROG0(I,M)+(QG(I,M)+DQVG1(I,M)+SUBB(I)*ZGI(I,M))*DTIME
        END DO
        CALL SGTSL(IMAXP1-1,C,B,A,D,INFO)
        IF(INFO.NE.0)THEN
         WRITE(2,*)'ERROR: Matrix inversion for volatiles (COMPOS). ',
     1             'INFO element of the diagonal=0. (INFO=',INFO,'). ', 
     2             'Exiting.'
         WRITE(17,*)NTIME
         WRITE(17,*)'ERROR: Matrix inversion for volatiles (COMPOS). ',
     1              'INFO element of the diagonal=0. (INFO=',INFO,'). ',
     2              'Exiting.'
         STOP
        ENDIF

! Calculate volatile ice densities
        DO I=2,IMAXP1
         ROG(I,M)=D(I-1)
         ROGI=ROG(I,M)
         ROSOGI1=ROSG0(I,M)+(DQVG2(I,M)*ROGI-DQVG1(I,M))*DTIME
         ROSOGI2=ROSG0(I,M)-QVG(I,M)*DTIME
         IF(I.EQ.IMAXP1)THEN
          ROSOGI=ROSOGI-VGFL(M)*DTIME/DV(I)
         ELSE
          ROSOGI=ROSOGI1
         END IF
         SUMSG(M)=SUMSG(M)+ROSOGI*DV(I)
         SUMG(M)=SUMG(M)+ROGI*DV(I)
         ROSG(I,M)=MAX(ZERO,ROSOGI)
         ROG(I,M)=MAX(ZERO,ROGI)
         SUMPSG(M)=SUMPSG(M)+ROSG(I,M)*DV(I)
         SUMPG(M)=SUMPG(M)+ROG(I,M)*DV(I)
        END DO

       END DO
      END IF

! Calculate dust densities
      DO I=2,IMAXP1
       RODI=ROD0(I)+(XMFL(I)-XMFL(I-1))*DTIME/DV(I)
       SUMD=SUMD+RODI*DV(I)
       ROD(I)=MAX(ZERO,RODI)
       SUMPD=SUMPD+ROD(I)*DV(I)
      END DO

! Normalize composition (to compensate for negative densities - set to zero)

      IF(SUMPC.GT.ZERO) FACNC=MAX(ZERO,SUMC/SUMPC)
      IF(SUMPB.GT.ZERO) FACNB=MAX(ZERO,SUMB/SUMPB)
      IF(SUMPA.GT.ZERO) FACNA=MAX(ZERO,SUMA/SUMPA)
      IF(SUMPV.GT.ZERO) FACNV=MAX(ZERO,SUMV/SUMPV)
      IF(SUMPD.GT.ZERO) FACND=MAX(ZERO,SUMD/SUMPD)
      FACNMIN=MIN(FACNC,FACNB,FACNA,FACNV,FACND)

      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        IF(SUMPG(M).GT.ZERO) FACNG(M)=MAX(ZERO,SUMG(M)/SUMPG(M))
        IF(SUMPSG(M).GT.ZERO) FACNSG(M)=MAX(ZERO,SUMSG(M)/SUMPSG(M))
        FACNMIN=MIN(FACNMIN,FACNG(M),FACNSG(M))
       END DO
      END IF

      IAC=1
      DO I=2,IMAXP1

       ROC(I)=ROC(I)*FACNC
       ROB(I)=ROB(I)*FACNB
       ROA(I)=ROA(I)*FACNA
       ROV(I)=ROV(I)*FACNV
       ROD(I)=ROD(I)*FACND
       RHOI=ROA(I)+ROB(I)+ROC(I)+ROV(I)+ROD(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         ROSG(I,M)=ROSG(I,M)*FACNSG(M)
         ROG(I,M)=ROG(I,M)*FACNG(M)
         RHOI=RHOI+ROSG(I,M)+ROG(I,M)
        END DO
       END IF

       XA(I)=ROA(I)/RHOI
       XB(I)=ROB(I)/RHOI
       XC(I)=ROC(I)/RHOI
       XV(I)=ROV(I)/RHOI
       XD(I)=ROD(I)/RHOI
       ROICE(I)=ROA(I)+ROB(I)+ROC(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         XSG(I,M)=ROSG(I,M)/RHOI
         XG(I,M)=ROG(I,M)/RHOI
         ROICE(I)=ROICE(I)+ROSG(I,M)
         ROSGI(M)=ROSG(I,M)
        END DO
       END IF

       DELMT=MAX(DELMT,ABS(RHOI-RHO(I)))
       IF((XA(I-1).GT.XC(I-1).AND.XA(I).LT.XC(I)).OR.
     1    (XA(I-1).LT.XC(I-1).AND.XA(I).GT.XC(I))) IAC=I
       DM(I)=DV(I)*RHOI
       PO(I)=POR(ROA(I),ROB(I),ROC(I),ROD(I),RHOICE,RHODUS,
     1           MMAX,ROSGI,RHOSOL)
       RHO(I)=RHOI

      END DO

      RETURN
      END
