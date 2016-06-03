      PROGRAM  COMPROG

!  APRIL 20, 1992  (Includes KNUDSEN and POISEUILLE flows)
!  JUNE 3, 1992    (Dimension extended to IMX=300)
!  APRIL 27, 2000  (Includes vector of volatiles, besides H2O)

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*20  COMIN,COMOUT,COMPLOT,COMPROF,COMODI,COMODO,COMLOG
      CHARACTER*10  FSUFF
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION  A(IMX),B(IMX),C(IMX),D(IMX),XT(IMX)
      DIMENSION  HV(IMX),HVG(IMX,LG),HG(IMX,LG)
      DIMENSION  DELU(IMX)
      DIMENSION  FLJV(IMX),FLJG(IMX,LG)
      DIMENSION  DFLJVL(IMX),DFLJGL(IMX,LG)
      DIMENSION  DFLJVR(IMX),DFLJGR(IMX,LG)
      INTEGER H,FN


!  Description of some parameters:
!  -------------------------------

!  Symbol "V" stands for H2O vapor.
!  Symbol "G" stands for volatile, other than H2O ("SG" for solid).
!  ROICE includes all ices (A+B+C+SUM(SG)).

!  NDEX messages:
!  0 - Normal.
!  1 - THETA < 1.
!  2 - No convergence for T.
!  3 - Negative abundance.
!  4 - No convergence for RHO.

!  Iteration chains:
!  KN=0 - Correct composition at each iteration.
!  KN>0 - Correct composition when T has converged.
!  KN>1 - Allow THETA LT 1 KN times per iteration.

!  Orbital angles (for a rotating comet, PROT > 0.0):
!  FI - Eccentric anomaly, measured from aphelion.
!  TANOM - True anomaly.
!  DECLIN - Declination angle (of spin axis); default 0.0. --> comin input.
!  ALAT - Cometocentric latitude angle; default 0.0. --> comin input.
!  JANGLE - The angle between the projection of the spin axis and the comet-sun line in perihelion, on the orbit plane; default 0.0. --> comin input.

!  Mapping:
!  --------
!  Calls: CONSTANTS, INPUT, ENDSTP, MODEL, SOURCE, FLUX, MFLUX, SGTSL, COMPOS, PRINT.
!  Called by: None.


      CALL GETARG(1,comin)
      CALL GETARG(2,comout)
      CALL GETARG(3,complot)
      CALL GETARG(4,comprof)
      CALL GETARG(5,comodi)
      CALL GETARG(6,comodo)
      CALL GETARG(7,comlog)
      
      OPEN(1,file=comin)
      OPEN(2,file=comout)
      OPEN(3,file=complot)
      OPEN(4,file=comprof)
      OPEN(7,file=comodi)
      OPEN(8,file=comodo)
      OPEN(17,file=comlog)
      
      CALL CONSTANTS
      
      CALL INPUT
      
!  Determining the number of hour angle elements (from comin file).
!  Opening the relevant number of h-angle files to output (20 to HANGN).

      IF(RUNM.GT.0)THEN

       IF(HANGN.GT.180)THEN
        WRITE(2,*) 'ERROR: input HANGN > 180. See comin file. Exiting.'
        WRITE(17,*) 'ERROR: input HANGN > 180. See comin file. Exiting.'
        STOP
       END IF

       IF(MOD(360,HANGN).NE.0)THEN
        DO WHILE(MOD(360,HANGN).NE.0)
         HANGN=HANGN+1
        END DO
       END IF

       DHANG=360.0D0/HANGN
       HANG=ZERO
       FN=19
       DO H=1,HANGN
        IF(INT(HANG).LT.10)THEN 
         WRITE(FSUFF,'(I1)') INT(HANG)
         OPEN(FN+H,file='hang_00'//FSUFF)
        ELSEIF((INT(HANG).GE.10).AND.(INT(HANG).LT.100))THEN 
         WRITE(FSUFF,'(I2)') INT(HANG)
         OPEN(FN+H,file='hang_0'//FSUFF)
        ELSEIF(INT(HANG).GE.100)THEN 
         WRITE(FSUFF,'(I3)') INT(HANG)
         OPEN(FN+H,file='hang_'//FSUFF)
        END IF
        HANG=HANG+DHANG
       END DO

      ELSE
       ALAT=0.0
       HANGN=1
       OPEN(20,file='hang_000')

      END IF
        
      NTIME=0
      NRED=0
      INDRED=0
      NIT=0
      MIT=0

!  Check spin period and that time step is smaller than the maximum allowed
      IF(ITMODE.GT.0)THEN
       IF(PROT.GT.ZERO)THEN
        FIXSTEP=PROT/DFLOAT(NTDAY)
        IF(FIXSTEP.GT.DTMAX)THEN
         WRITE(2,*)'ERROR: time-step greater than maximum ', 
     1             '(FIXSTEP > DTMAX). Exiting.'
         WRITE(17,*)NTIME
         WRITE(17,*)'ERROR: time-step greater than maximum ',
     1              '(FIXSTEP > DTMAX). Exiting.'
         STOP
        END IF
       ELSE
        FIXSTEP=DTMAX
       END IF
      END IF

      CALL ENDSTP(IPRINT,INDRED,ITIMEA,ITIMER,NMERGE,NRED,
     1            NIT,MIT,DIST0)

      FLJV(1)=ZERO
      DFLJVL(1)=ZERO
      DO M=1,LG
       FLJG(1,M)=ZERO
       DFLJGL(1,M)=ZERO
      END DO

!  Begin time steps

      DO  1000 NTIM=1,NTIMAX

      NTIME=NTIM
      NRED=0

      DO I=2,IMAXP1
       TI=T0(I)
       XM0(I)=ONE/(ONE+EXP(EXMELT*(ONE-TI/TMELT))) 
       HV(I)=HSUB+TI*(CV-BETA-ALFA*TI)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         HVG(I,M)=HSUBG(M)+TI*(CG(M)-BETA-ALFA*TI)
         HG(I,M)=TI*(CG(M)-BETA-ALFA*TI)
        END DO
       END IF
      END DO

 100  MIT=0
      DELMT =ZERO
      DELMAX=ZERO

!  Determine time step

      IF(ITMODE.EQ.0)THEN
       FI=FI+DFI
       IF(ECCEN.EQ.ZERO)THEN
        DTIME=PSICON*DFI
        TIME=TIME+DTIME
       ELSE
        TIME=PSICON*(FI-ECCEN*SIN(FI)-PI)
        DTIME=TIME-TIME0
        IF(DTIME.LE.ZERO) DTIME=DFI*(PSICON*(ONE-ECCEN*COS(FI)))
       END IF
      ELSE
       DTIME=FIXSTEP
       TIME=TIME+DTIME
       DFI=DTIME/(PSICON*(ONE-ECCEN*COS(FI)))
       FI=FI+DFI
      END IF

      DIST=AXIS*(ONE-ECCEN*COS(FI))

      IF(PROT.GT.ZERO)THEN
       TANOM=ASIN(AXIS*SQRT(ONE-ECCEN**2)/DIST*SIN(FI-PI))
       PHROT=TWO*PI*MOD(TIME/PROT,ONE)
!       HANGLE=PHROT+HOUR
       HANGLE=PHROT
       OTHTA=TWO*ATAN(SQRT((ONE+ECCEN)/(ONE-ECCEN))*TAN(FI/TWO))
       SSTHTA=ACOS(SIN(DECLIN)*COS(JANGLE+OTHTA))
       IF(SSTHTA.LT.TWO*PI) SSTHTA=SSTHTA+TWO*PI
       TEMPX=SIN(PHROT)*SIN(OTHTA-JANGLE)
     1       +COS(PHROT)*COS(OTHTA+JANGLE)*COS(DECLIN)
       TEMPY=COS(PHROT)*SIN(OTHTA-JANGLE)
     1       -SIN(PHROT)*COS(OTHTA+JANGLE)*COS(DECLIN)
!       SSFI=ATAN(TEMPY/TEMPX)
       SSFI=ZERO
       IF((TEMPX.LT.ZERO).AND.(TEMPY.GT.ZERO)) SSFI=PI-SSFI
       IF((TEMPX.LT.ZERO).AND.(TEMPY.LT.ZERO)) SSFI=PI+SSFI
       IF((TEMPX.GT.ZERO).AND.(TEMPY.LT.ZERO)) SSFI=TWO*PI-SSFI
       CFI=HANGLE
       CTHTA=ABS(PI/TWO-ALAT)
       COSXI=SIN(SSTHTA)*SIN(CTHTA)*COS(CFI)+
     1       COS(SSTHTA)*COS(CTHTA)
        IF(COSXI.LT.ZERO) COSXI=ZERO
        SFLUX=MAX(ZERO,Q/(DIST*DIST)*COSXI)
      ELSE
        SFLUX=Q/(DIST*DIST)
      END IF

      IF(DTIME.GT.DTMAX)THEN
       FI=FI-DFI
       TIME=TIME0
       DFI=DFI/FACTIM
       GO TO 100 
      END IF

!  Start iterations for T(I)-NIT and RHO(I)-MIT

!  Initialization of parameters

  101 MIT=MIT+1
      DO I=2,IMAXP1
       XTP(I)=ZERO
       XT(I)=ZERO
      END DO
      NDEX=0
      NIT=0
      ACON=ZERO

!  New iteration --> Solve the equations throughout all the shells

    3 NIT=NIT+1
      ACONP=ACON

      CALL SOURCE
      CALL FLUX
      CALL MFLUX

      DO I=2,IMAXP1

       TI=T(I)
       T0I=T0(I)
       UCON=ROV(I)*CV+ROD(I)*CD
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         UCON=UCON+ROG(I,M)*CG(M)
        END DO
       END IF

       IF(TI.LT.TML)THEN
        DELU(I)=(ROICE(I)*(ALFA*(TI+T0I)+BETA)+UCON)*(TI-T0I)
        DUDT(I)=ROICE(I)*(TWO*ALFA*TI+BETA)+UCON
       ELSE
        XM=ONE/(ONE+EXP(EXMELT*(ONE-TI/TMELT)))
        UICE=(ONE-XM)*(ALFA*TI*TI+BETA*TI)+XM*(HMELT+CW*TI)
        UICE0=(ONE-XM0(I))*(ALFA*T0I*T0I+BETA*T0I)+
     1        XM0(I)*(HMELT+CW*T0I)
        DELU(I)=ROICE(I)*(UICE-UICE0)+UCON*(TI-T0I)
        DUDT(I)=ROICE(I)*((ONE-XM)*(TWO*ALFA*TI+BETA)+XM*CW+
     1    XM*(ONE-XM)*EXMELT/TMELT*(HMELT+CW*TI-ALFA*TI*TI-BETA*TI))+
     2    UCON
       END IF

       IF(I.LT.IMAXP1)THEN
        TIP=T(I+1)
        IF(VAFL(I).GT.ZERO)THEN
         FLJV(I)=VAFL(I)*CV*TIP
         DFLJVR(I)=VAFL(I)*CV+DVFDTR(I)*CV*TIP
         DFLJVL(I)=DVFDTL(I)*CV*TIP
        ELSE
         FLJV(I)=VAFL(I)*CV*TI
         DFLJVR(I)=DVFDTR(I)*CV*TI
         DFLJVL(I)=VAFL(I)*CV+DVFDTL(I)*CV*TI
        END IF
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          IF(GAFL(I,M).GT.ZERO)THEN
           FLJG(I,M)=GAFL(I,M)*CG(M)*TIP
           DFLJGR(I,M)=GAFL(I,M)*CG(M)+DGFDTR(I,M)*CG(M)*TIP
           DFLJGL(I,M)=DGFDTL(I,M)*CG(M)*TIP
          ELSE
           FLJG(I,M)=GAFL(I,M)*CG(M)*TI
           DFLJGR(I,M)=DGFDTR(I,M)*CG(M)*TI
           DFLJGL(I,M)=GAFL(I,M)*CG(M)+DGFDTL(I,M)*CG(M)*TI
          END IF
         END DO
        END IF
       ELSE
        FLJV(I)=VAFL(I)*CV*TI
        DFLJVR(I)=ZERO
        DFLJVL(I)=VAFL(I)*CV+DVFDTL(I)*CV*TI
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          FLJG(I,M)=GAFL(I,M)*CG(M)*TI
          DFLJGR(I,M)=ZERO
          DFLJGL(I,M)=GAFL(I,M)*CG(M)+DGFDTL(I,M)*CG(M)*TI
         END DO
        END IF
       END IF
      END DO

      ACON=ZERO
      DO I=2,IMAXP1

       IM=I-1
       DVI=DV(I)
       DVOT=DVI/DTIME
       AIM=DFLR(I)+DFLJVR(I)
       B(IM)=-DFLL(I)+DFLR(IM)-DFLJVL(I)+DFLJVR(IM)
     1       +DUDT(I)*DVOT+(HV(I)*DQVDT(I)-DSDT(I))*DVI
       CIM=-DFLL(IM)-DFLJVL(IM)
       D(IM)=FL(I)-FL(IM)+FLJV(I)-FLJV(IM)
     1      -DELU(I)*DVOT+(S(I)-HV(I)*QV(I))*DVI
       ELHS(I)=(FL(I)-FL(IM))/DVI
       ERHS(I)=(FLJV(I)-FLJV(IM))/DVI
       SVG(I)=-HV(I)*QV(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         AIM=AIM+DFLJGR(I,M)
         B(IM)=B(IM)-DFLJGL(I,M)+DFLJGR(IM,M)
     1        +(HVG(I,M)*DQVGDT(I,M)+HG(I,M)*DQGDT(I,M))*DVI
         CIM=CIM-DFLJGL(IM,M)
         D(IM)=D(IM)+FLJG(I,M)-FLJG(IM,M)
     1        -(HVG(I,M)*QVG(I,M)+HG(I,M)*QG(I,M))*DVI
         ERHS(I)=ERHS(I)+(FLJG(I,M)-FLJG(IM,M))/DVI
         SVG(I)=SVG(I)-(HVG(I,M)*QVG(I,M)+HG(I,M)*QG(I,M))
        END DO
       END IF
       A(IM)=-AIM
       C(IM)=-CIM
       DUDT(I)=DELU(I)/DTIME
       RATFL(I)=ZERO
       IF(ELHS(I).NE.ZERO) RATFL(I)=ERHS(I)/ELHS(I)
       ACON=ACON+ABS(D(IM))

      END DO

      CALL SGTSL(IMAXP1-1,C,B,A,D,INFO)
      IF(INFO.NE.0)THEN
       WRITE(2,*)'ERROR: Matrix inversion (MAIN). ',
     1           'INFO element of the diagonal=0 (INFO=',INFO,'). ', 
     2           'Exiting.'
       WRITE(17,*)NTIME
       WRITE(17,*)'ERROR: Matrix inversion (MAIN). ',
     1            'INFO element of the diagonal=0 (INFO=',INFO,'). ', 
     2            'Exiting.'
       STOP
      END IF

!  Convergence check

      ITHETA=0
      THETA=ONE
  611 DO I=2,IMAXP1
       XT(I)=D(I-1)
       DELT=XT(I)*THETA
       DELX1=ABS(DELT)
       IF(DELX1/T0(I).GT.2.D1*EPSTIM)THEN
        THETA=THETA*0.5D0
        ITHETA=ITHETA+1
        NDEX=1
        IF(ITHETA.GE.KN)THEN
         GO TO 7
        ELSE
         GO TO 611
        END IF
       END IF
      END DO

      KCONV=0
      DELXMA=ZERO
      DELMAX1=ZERO
      IDELM=1
      DO I=2,IMAXP1
       DELT=XT(I)*THETA
       DELX1=ABS(DELT)
       DELXMA=MAX(DELXMA,DELX1/T(I))
       TI=T(I)+DELT
       IF(TI.LE.ZERO)THEN
        TI=T(I)*EXP(DELT/T(I))
        IF(TI.LT.TIMIN)THEN
         DELTT=XTP(I)+DELT
         TI=T0(I)*EXP(DELTT/T0(I))
        END IF
       END IF
       T(I)=TI
       IF(DELX1.GT.EPSX1*TI) KCONV=KCONV+1
       XTP(I)=XTP(I)+DELT
       DELMAX=MAX(DELMAX,ABS((TI-T0(I))/T0(I)))
       DELMP=DELMAX1
       DELMAX1=MAX(DELMAX1,ABS(XTP(I)/T0(I)))
       IF(DELMAX1.NE.DELMP) IDELM=I
      END DO

      IF(KN.EQ.0)THEN
       DELMTT=DELMT
       CALL COMPOS
       IF(NDEX.EQ.3) GO TO 7
      END IF
      IF(IPRINT.NE.0.OR.MOD(NTIME,NPRINT).EQ.0)THEN
       WRITE(2,'(2X,6I5,1P5E15.3)')NTIME,MIT,NIT,KCONV,IDELM,ITIMEA,
     1                             THETA,DELMT,DELMAX1,DELXMA,ACON
      END IF

      IF(KCONV.EQ.0)THEN
!  Another iteration for T or repeat time step
       IF(NIT.LT.KMAX) GO TO 3
       IF(NIT.LT.5*KMAX.AND.ACON.LT.ACONP) GO TO 3
       NDEX=2
      ELSE
       IF(KN.GT.0)THEN
!  Correct composition and re-iterate
        DELMTT=DELMT
        CALL COMPOS
        IF(NDEX.EQ.3) GO TO 7
       END IF
       IF(DELMT.LT.EPSX2) GO TO 9
       IF(DELMTT.GT.ZERO.AND.DELMT/DELMTT.LT.5.D-2) GO TO 9
       IF(NMERGE.GT.0) GO TO 9
       IF(KN.GT.0.AND.MIT.LT.KMAX) GO TO 101
       NDEX=4
      END IF

!  Reduce time step and repeat calculations

    7 TIME=TIME-DTIME
      FI=FI-DFI
      DFI=DFI*FAC
      NRED=NRED+1
      INDRED=NDEX
      ITIMER=ITIMER+1
      ITIMEA=1
!      IF(MOD(ITIMER,10).EQ.0) DTMAX=DTMAX/1.D1

      IF(NRED.LE.10)THEN
        DO I=2,IMAXP1
         T(I)=T0(I)
         RHO(I)=RHO0(I)
         ROICE(I)=ROICE0(I)
         ROD(I)=ROD0(I)
         ROA(I)=ROA0(I)
         ROB(I)=ROB0(I)
         ROC(I)=ROC0(I)
         ROV(I)=ROV0(I)
         IF(MMAX.GT.0)THEN
          DO M=1,MMAX
           ROG(I,M)=ROG0(I,M)
           ROSG(I,M)=ROSG0(I,M)
          END DO
         END IF
        END DO
        GO TO 100
      END IF

      WRITE(2,50)NTIME,NDEX,KCONV,NIT,MIT,DFI,DELMT,DELMTT,DELXMA,
     1           DELMAX,DTIME,ACON,ACONP
      WRITE(17,*)NTIME
      WRITE(17,50)NTIME,NDEX,KCONV,NIT,MIT,DFI,DELMT,DELMTT,DELXMA,
     1            DELMAX,DTIME,ACON,ACONP,IB
   50 FORMAT(1X,' TOO MANY REPETITIONS - NTIME',I5,' NDEX',I5,
     1      ' KCONV',I5,' NIT',I5,' MIT',I5/1X,'DFI',E12.5,' DELMT',
     2      E12.5,' DELMTT',E12.5,' DELXMA',E12.5,' DELMAX',E12.5/1X,
     3      ' DTIME',E12.5,' ACON',E12.5,' ACONP',E12.5)

      CALL MODEL
      CALL PRINT

      STOP

!  End of time step

!  Calculate total sublimation energy for time step

    9 QVTOT=ZERO
      DO I=2,IMAXP1
       QVTOT=QVTOT+QV(I)*DV(I)*HV(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         QVTOT=QVTOT+QVG(I,M)*DV(I)*HVG(I,M)
        END DO
       END IF
      END DO

      CALL ENDSTP(IPRINT,INDRED,ITIMEA,ITIMER,NMERGE,NRED,
     1            NIT,MIT,DIST0)

!  NDEX=9 (set in ENDSTP) if TIME exceeds TIMAX (set in CONSTANTS)

      IF(NDEX.EQ.9) STOP

!  End of time steps loop

 1000 CONTINUE 

      STOP

      END
