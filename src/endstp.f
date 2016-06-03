!***************************************************************
      SUBROUTINE ENDSTP(IPRINT,INDRED,ITIMEA,ITIMER,NMERGE,NRED,
     1                  NIT,MIT,DIST0)
!***************************************************************

!       Subroutine description:
!       -----------------------
!       Determines parameters for the next time step calculations.

!       Subroutine mapping:
!       -------------------
!       Calls: PRINT, PROFILE, SURFL, DIVIS, SMERGE, GRIDEF, HYDRO, MODEL, IMPACT.
!       Called by: MAIN.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0,FOUR=4.D0)
      DIMENSION  TBREAK(IMX),TDYN(IMX),FM(IMX)
      DIMENSION  RHOT(IMX),PT(IMX),DMT(IMX),DVT(IMX)
      DIMENSION  DGTOT(LG),FMG(LG),FMG0(LG)
      DIMENSION  SUMG(LG),VGIP(LG),DEPTHG(LG)
      SAVE


      IF(NTIME.EQ.0)THEN

!  Set initial global parameters

       UTOTP=ZERO
       UINPUT=ZERO
       VFOUT=ZERO
       FMH2O0=ZERO
       FMDUS0=ZERO
       SUMDUS=ZERO
       SUMVAP=ZERO
       DO M=1,LG
        FMG0(M)=ZERO
        SUMG(M)=ZERO
        DEPTHG(M)=ZERO
        IBG(M)=1
       END DO
       IF(NPRINT.LT.0)THEN
        IPRINT=1
        NPRINT=-NPRINT
       ELSE
        IPRINT=0
       END IF
       ITIMER=0
       ITIMEA=1
       NMERGE=0
       DIST0=DIST
       NORBPREV=1
       PI2=TWO*PI
       EROT=PI2/3.6D2
       FIROT=ZERO
       ILEG=-1
       IF(IOP.GT.0)THEN
        AXIS=AXISFIT(TIME0)*ONEAU
        ECCEN=ECCENFIT(TIME0)
        WRITE(2,*) 'New orbital parameters: AXIS, ECCEN = ', AXIS,ECCEN
        WRITE(17,*)NTIME
        WRITE(17,*) 'New orbital parameters: (a,e) = ', AXIS,ECCEN
        PSICON=(AXIS**(THRE/TWO))/SQRT(GM)
       END IF
     
       HOUR(1)=ZERO
       DO J=2,HANGN
        HOUR(J)=HOUR(J-1)+2.0D0*PI/REAL(HANGN)
       END DO
      
!  It is implied that calculations start on the inbound leg

       IF(IHYDRO.GT.0)THEN
!  Calculate coefficients of the equation of state
        DO I=1,IMAXP1
         TBREAK(I)=ZERO
         RHOT(I)=RHO(I)
         PT(I)=P(I)
         DMT(I)=DM(I)
         DVT(I)=DV(I)
        END DO
        RADIUS=R(IMAXP1)
        RO0=RO0FAC*RHO(2)
        PK=THRE/(TWO*PI4)*G*FMTOT0**2/RADIUS**4
        IEOS=IMAXP1
        WRITE(2,*)'EOS parameters set: R,PK,RO0 =  ',RADIUS,PK,RO0
        WRITE(17,*)'EOS parameters set: R,PK,RO0 =  ',RADIUS,PK,RO0
        CALL SETEOS(PSOL,RHOT,PT,DMT,DVT,IEOS,G,PI4,PK,RO0,RADIUS
     1              EPSX2,EPSHY)
        IF(IEOS.EQ.0)THEN
         IHYDRO=0
         WRITE(2,*)'ERROR: No solution for EOS. Set IHYDRO=0.'
         WRITE(17,*)'ERROR: No solution for EOS. Set IHYDRO=0.'
        END IF
       END IF 

       RETURN

      END IF

!  If NTIME is not equal to 0:
!  Calculate global quantities
!  Calaculate the first depth in which each specie resides

      UTOT=ZERO
      FM(1)=ZERO
      FMDUS=ZERO
      FMH2O=ZERO
      EPSI=ZERO
      EPSD=ZERO
!  If there is more than 0 gas species:
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        FMG(M)=ZERO
        DEPTHG(M)=ZERO
       END DO
      END IF
      PMAX=ZERO
      DEPTH=ZERO
      TDYNMI=1.D20
      TETMAX=ZERO
      DEPTET=ZERO

!  Run on the number of grid shells (inner shells):
      IP=IMAXP1
      DO I=2,IP

       TI=T(I)
       RHOI=RHO(I)
       DVI=DV(I)
       DM(I)=RHOI*DVI
       DMI=DM(I)
       ZGAS=ZERO
       SUMX=ROV(I)/AIB
       CGI=ROV(I)*CV+ROD(I)*CD
       IF(MMAX.GT.0)THEN
        IR=IP+2-I
!  IR starts from the outer shell towards the core.
!  Go over the gas species:
        DO M=1,MMAX
         ZGIM=ZGI(I,M)
         ZGAS=ZGAS+ZGIM
         SUMX=SUMX+ROG(I,M)/AGB(M)
         CGI=CGI+ROG(I,M)*CG(M)
!  ROG is the mass fraction of the gas multiplied by the density.
         FMG(M)=FMG(M)+((ROA(I)+ROB(I))*ZGIM+ROG(I,M)+ROSG(I,M))*DVI
         IF(ROSG(IR,M).LT.PRECX)THEN
          DEPTHG(M)=(R(IP)-R(IR-1))/1.D2
          IBG(M)=IR-1
         END IF
        END DO
       END IF
       P(I)=TI*SUMX/PO(I)
       IF(TI.LT.TML)THEN
        UICE=ALFA*TI*TI+BETA*TI
       ELSE
        XM=ONE/(ONE+EXP(EXMELT*(ONE-TI/TMELT)))
        UICE=(ONE-XM)*(ALFA*TI*TI+BETA*TI)+XM*(HMELT+CW*TI)
       END IF
       U(I)=ROICE(I)*UICE+CGI*TI
       UTOT=UTOT+U(I)*DVI
       FM(I)=FM(I-1)+DMI
       FMDUS=FMDUS+ROD(I)*DVI
       FMH2O=FMH2O+((ROA(I)+ROB(I))*(ONE-ZGAS)+ROC(I)+ROV(I))*DVI
       IF(P(I).GT.PMAX) PMAX=P(I)
       PCRIT=PCR(I)
       IF(P(I).GT.PCRIT.AND.P(I-1).LT.PCRIT) DEPTH=(R(IP)-R(I-1))/1.D2
       TDYN(I)=TDYNMI
       DRAINV=ABS(VAFL(I)-VAFL(I-1))
       IF(DRAINV.GT.ZERO) TDYN(I)=ROV(I)*DV(I)/DRAINV
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         DRAING=ABS(GAFL(I,M)-GAFL(I-1,M))
         IF(DRAING.GT.ZERO) TDYN(I)=MAX(TDYN(I),ROG(I,M)*DV(I)/DRAING)
        END DO
       END IF
       P1=P(I)
       R1=R(I-1)
       R12=R(I)
       IF(I.LT.IP)THEN
        P2=P(I+1)
        R2=R(I+1)
       ELSE
        P2=ZERO
        R2=R(IP)
       END IF
!      A1=0.25D0*(P1*(R1+R12)**2-P2*(R12+R2)**2)/(R2-R1)
!      A2=PCR(I)*R12
       A1=(P1-P2)/(R2-R1)
       A2=4.D0*PCR(I)/R12
       IF(A1.GT.A2)THEN
        TBREAK(I)=TBREAK(I)+DTIME
        IF(TBREAK(I).GT.TDYN(I))THEN
         TETI=MIN(TETA*XPS,TET(I)*A1/A2)
         ALFP(I)=ALFP(I)+LOG(TETI/TET(I))/LOG(XPS)
         TET(I)=TETI
         PCR(I)=PCR(I)*SQRT(A1/A2)
         TBREAK(I)=ZERO
        END IF
       ELSE
        TBREAK(I)=ZERO
        TETI=MIN(TETA*XPS,TET(I)*
     1      (ONE+DTIME*1.26D5*(P(I)/PCRIT)**4*EXP(-9.1D11/(RGAS*T(I)))))
        ALFP(I)=ALFP(I)+LOG(TETI/TET(I))/LOG(XPS)
        TET(I)=TETI
       END IF
       TDYNMI=MIN(TDYNMI,TDYN(I))
       TETMAX=MAX(TETMAX,TET(I))
       IF(TETMAX.EQ.TET(I)) DEPTET=R(IP)-R(I)

!  End of running on all inner grid shells
      END DO

      FMTOT=FM(IP)
      VFIP=VAFL(IP)-VAPFL-VAPFLA
      SUMVAP=SUMVAP-VFIP*DTIME
      DVTOT=-VFIP/AMICE
      DELI1=FMH2O-FMH2O0
      DELI2=VFIP*DTIME
      IF(DELI1*DELI2.NE.ZERO) EPSI=0.5D0*(DELI1-DELI2)/(DELI1+DELI2)
      VDIP=XMFL(IP)
      SUMDUS=SUMDUS-VDIP*DTIME
      DDOT=-VDIP
      DELD1=FMDUS-FMDUS0
      DELD2=VDIP*DTIME
      IF(DELD1*DELD2.NE.ZERO) EPSD=0.5D0*(DELD1-DELD2)/(DELD1+DELD2)
      DELUIN=(FL(IP)+(CV*VAFL(IP)+CD*VDIP)*T(IP)+STOT-QVTOT)*DTIME
      VFOUT=VFOUT+(VFIP+VDIP)*DTIME
      EPSG=ZERO
      MCR=0
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        VGIP(M)=GAFL(IP,M)-VGFL(M)-VGFLA(M)
        SUMG(M)=SUMG(M)-VGIP(M)*DTIME
        DELUIN=DELUIN+CG(M)*GAFL(IP,M)*T(IP)*DTIME
        VFOUT=VFOUT+VGIP(M)*DTIME
        DGTOT(M)=-VGIP(M)/AMG(M)
        DELG1=FMG(M)-FMG0(M)
        DELG2=VGIP(M)*DTIME
        IF(DELG1*DELG2.NE.ZERO)THEN
         EPSG1=MAX(EPSG,ABS(0.5D0*(DELG1-DELG2)/(DELG1+DELG2)))
         IF(EPSG1.GT.EPSG) MCR=M
         EPSG=EPSG1
        END IF 
       END DO
      END IF
      UINPUT=UINPUT+DELUIN
      EPSEN=(UTOT-UTOTP-DELUIN)/(UTOT-UTOTP+DELUIN)*0.5D0
      DISTAU=DIST/ONEAU
      XORB=0.5D0*(ONE+FI/PI)
      NORB=INT(XORB)
      ILEG0=ILEG
      IF(DIST.LT.DIST0) ILEG=-1
      IF(DIST.GT.DIST0) ILEG=1
!      IF(ILEG.NE.ILEG0) CALL PRINT
      DIST0=DIST

!  Print time step data

! The following statements (if activated) are meant to record perihelion properties

!      PHDIST=AXIS*(ONE-ECCEN)
!      IF(ILEG.EQ.-1.AND.ABS(DIST-PHDIST)/PHDIST.LT.5.D-3)THEN
!       PHROTH=MOD(TIME/PROT,ONE)*360.D0
!       WRITE(3,66) NORB,TIME,DISTAU,PHROTH,T(IP),ROC(IP),ROC(IP-1),
!    1              DVTOT/FACR(IP),DGTOT(1)/FACR(IP),DGTOT(2)/FACR(IP),
!    2              DGTOT(3)/FACR(IP),DGTOT(4)/FACR(IP),
!    3              DGTOT(5)/FACR(IP),DDOT(FACR(IP)
!  66   FORMAT(1X,I5,1PE17.9,1PE14.6,1P11E10.3)
!       CALL PROFILE
!      END IF

      IF(IIMPAC.GT.0)THEN
        CALL IMPACT
       ELSE
        FIMPACT=ZERO
       END IF

      TIP=T(IMAXP1)
      IF(MOD(NTIME,NSTEP*10).EQ.0.OR.NTIME.EQ.1) WRITE(2,51)
  51  FORMAT(50H**************************************************/1X,
     1'     NT IB IN  M  N  AC  G1  G2    TIME    D(AU)/TC   T-SURF',4X,
     2'H2ODOT    CODOT     CO2DOT   DUSTDOT    DELMAX    SFLUX      XI 
     3    DTIME '//)

      IF(MOD(NTIME,NSTEP).EQ.0)THEN
       IF(ECCEN.GT.ZERO)THEN
        DISTC=DISTAU
       ELSE
        DISTC=T(2)
       END IF
       SFLPR=FL(IMAXP1)/FACR(IMAXP1)
       INDEP=NDEX
       IF(NRED.GT.0) INDEP=INDRED
       WRITE(2,52)NTIME,IB  ,INDEP,MIT,NIT,IAC,IBG(1),IBG(2),TIME,DISTC,
     1            TIP,DVTOT,DGTOT(1),DGTOT(2),DDOT,DELMAX,SFLPR,
     2            ACOS(COSXI)*180.D0/PI,DTIME
  52   FORMAT(I9,I2,3I3,3I4,1P11E10.2)
      END IF

      TMAX=ZERO
      DO I=2,IMAXP1
       TMAX=MAX(TMAX,T(I))
      END DO
      IF(ECCEN.GT.ZERO)THEN
       TEMP=TMAX
       DISTM=DISTAU
      ELSE
       TEMP=T(2)
       DISTM=TMAX
      END IF
      DRIAC=ZERO
      IF(R(IAC).GT.ZERO) DRIAC=(R(IMAXP1)-R(IAC))

      IF(PROT.GT.ZERO)THEN

       TANOM=ASIN(AXIS*SQRT(ONE-ECCEN**2)/DIST*SIN(FI-PI))
       PHROT=TWO*PI*MOD(TIME/PROT,ONE)
!       HANGLE=PHROT-TANOM
       HANGLE=PHROT
!       DELPD=ABS(TANOM-PHROT)
       DELPD=ABS(PHROT)
!       DELPN=ABS(TANOM-(PHROT-PI))
       DELPN=ABS(PHROT-PI)
       
       IF(IIMPAC.GT.0)THEN
        CALL IMPACT
       ELSE
        FIMPACT=ZERO
       END IF

!  In order to write the h-files, we need to know which 
!  part of the self-spin the comet is at. HANGLE denotes
!  the degrees by which the comet is turned, relative to
!  the starting point

       JJ=0
       DO J=1,HANGN-1
        IF(ABS(HANGLE-HOUR(J)).LT.PI/180.D0) JJ=19+J
       END DO
       IF((ONE-COS(HANGLE)).LT.1.5D-4) JJ=19+HANGN
       IF(JJ.GT.0)THEN
        WRITE(JJ,*)NORB,TIME,DISTAU,TIP,DVTOT,DDOT,
     1             DGTOT(1),DGTOT(2),DGTOT(3),DGTOT(4),DGTOT(5),
     2             DGTOT(6),DGTOT(7),DGTOT(8),DGTOT(9),DGTOT(10),
     3             DRIAC,R(IMAXP1),HANGLE/PI*180.D0
       END IF
      ELSE
       IF(MOD(NTIME,NSTEP).EQ.0.OR.DELMAX.GT.1.D-2*EPSTIM)THEN
        WRITE(3,60)NORB,TIME,DISTM,TEMP,TIP,DVTOT,
     1             DGTOT(1),DGTOT(2),DDOT,DEPTHG(1),DEPTHG(2),
     2             DRIAC,R(IMAXP1)
       END IF
  60   FORMAT(1X,I5,1PE15.7,1PE13.5,1P8E10.2,1PE10.3,1PE12.5,1PE10.3)
!  60   FORMAT(1X,I5,1PE15.7,1PE13.5,1P13E10.2,1PE10.3,1PE12.5,1PE10.3)
      END IF

      IF(MOD(NTIME,NPRINT).EQ.0)THEN
       CALL PRINT
       UGAIN=UTOT-UTOT0
       DELTAU=ABS(UGAIN-UINPUT)/(0.5D0*(ABS(UGAIN)+ABS(UINPUT)))
       WRITE(2,53)UGAIN,UINPUT,DELTAU,EPSEN
  53   FORMAT(1X,'ENERGY BALANCE- CHANGE IN ENERGY=',1PE11.3,
     1 ' INTEGRATED GAINS AND LOSSES=',1PE11.3,' CONSERVATION ACCURACY='
     2 ,1P2E11.3)
       FMGAIN=FMTOT-FMTOT0
       DELTAM=ABS(FMGAIN-VFOUT)/(0.5D0*(ABS(FMGAIN)+ABS(VFOUT)))
       WRITE(2,54)FMTOT,FMGAIN,VFOUT,DELTAM,EPSI,EPSD,EPSG,MCR
  54   FORMAT(1X,'MASS BALANCE- MASS=',1PE12.5,' MASS CHANGE=',
     1      1PE11.4,' LOST=',1PE11.4,' ACCURACY=',1PE10.3,
     2      1P3E11.3,I3)
       CALL PROFILE
      END IF


! Set values for the next time step

      IF(TIME.GT.TIMAX)THEN
       WRITE(2,*)'ERROR: time greater than maximum (TIME > TIMAX).'
       WRITE(17,*) NTIME
       WRITE(17,*)'ERROR: time greater than maximum (TIME > TIMAX). '
       NDEX=9
        STOP
      END IF

      IF(DELMAX.LT.EPSTIM/THRE) DFI=MIN(DFIMAX,DFI*FACTIM)
      IF(DELMAX.GT.THRE*EPSTIM.AND.NMERGE.EQ.0) DFI=DFI/FACTIM
      IF(DELMAX.LT.1.D-2) ITIMEA=ITIMEA+1
!      IF(MOD(ITIMEA,10).EQ.0) DTMAX=MIN(DTMAX*1.D1,TIMAX)
      TIME0=TIME
      DIST0=DIST
      IF(NORB.NE.NORBPREV)THEN
       NORBPREV=NORB
       IF(IOP.GT.0)THEN
        AXIS=AXISFIT(TIME0)*ONEAU
        ECCEN=ECCENFIT(TIME0)
        WRITE(2,*) 'New orbital parameters: AXIS, ECCEN =  ',AXIS,ECCEN
        WRITE(17,*)NTIME
        WRITE(17,*) 'New orbital parameters: (a,e) =  ',AXIS,ECCEN
        PSICON=(AXIS**(THRE/TWO))/SQRT(GM)
       END IF
      END IF

      NMERGE=0
      I=IMAXP1
      IM=I-1
      DMI=DM(I)
      DMM=DM(IM)
      IF(DMI/DMM.LT.0.4D0)THEN
       IF(T(I)/T(IM).LT.1.25D0.OR.T(I)/T(IM).GT.0.8D0)THEN
        MODIV=0
!  Only DIVIS is applicable here, since the last shell is dealt with.
!        CALL DIVIS(MODIV,I)
        NMERGE=1
       END IF
      END IF

      IP1=IMAXP1
      DO I=3,IP1-1

       IF(I.GT.IMAXP1-1) GO TO 567

!  This question is necessary since SMERGE reduces the number of shells
!  Re-divide the grid shells in case there are too big differences
!  between consecutive shells.

       IM=I-1
       IP=I+1
       DMI=DM(I)
       DMM=DM(IM)
       DMP=DM(IP)
!       IF(DMI/DMP.LT.0.4D0.OR.DMI/DMM.LT.0.4D0)THEN
       IF(DMI/DMP.LT.0.1D0.OR.DMI/DMM.LT.0.1D0)THEN
        IF(DMP.GT.DMM)THEN
         IF(T(I)/T(IP).LT.1.25D0.OR.T(I)/T(IP).GT.0.8D0)THEN
          MODIV=0
          IF(IZONE1.EQ.0)THEN
           WRITE(2,*) 'SMERGE call. I,DM,I-1,DM,I+1,DM,T(I)/T(I+1) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(I)/T(IP)
           WRITE(17,*) NTIME
           WRITE(17,*) 'SMERGE call. I,DM,I-1,DM,I+1,DM,T(I)/T(I+1) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(I)/T(IP)
           CALL SMERGE(I)
           NMERGE=1
          ELSE
           WRITE(2,*) 'DIVIS call. I,DM,I-1,DM,I+1,DM,T(I)/T(I+1) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(I)/T(IP)
           WRITE(17,*) NTIME
           WRITE(17,*) 'DIVIS call. I,DM,I-1,DM,I+1,DM,T(I)/T(I+1) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(I)/T(IP)
           CALL DIVIS(MODIV,IP)
           NMERGE=1
          END IF
         END IF
        ELSE
         IF(T(IM)/T(I).LT.1.25D0.OR.T(IM)/T(I).GT.0.8D0)THEN
          MODIV=0
          IF(IZONE1.EQ.0)THEN
           WRITE(2,*) 'SMERGE call. I,DM,I-1,DM,I+1,DM,T(I-1)/T(I) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(IM)/T(I)
           WRITE(17,*) NTIME
           WRITE(17,*) 'SMERGE call. I,DM,I-1,DM,I+1,DM,T(I-1)/T(I) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(IM)/T(I)
           CALL SMERGE(I)
           NMERGE=1
          ELSE
           WRITE(2,*) 'DIVIS call. I,DM,I-1,DM,I+1,DM,T(I-1)/T(I) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(IM)/T(I)
           WRITE(17,*) NTIME
           WRITE(17,*) 'DIVIS call. I,DM,I-1,DM,I+1,DM,T(I-1)/T(I) = ',
     1                 I,DMI,IM,DMM,IP,DMP,T(IM)/T(I)
           CALL DIVIS(MODIV,I)
           NMERGE=1
          END IF
         END IF
        END IF
       END IF
      END DO

!  Calculate erosion of nucleus

!  GRIDEF should be used only for a heavily eroding nucleus (e.g., JFCs, sun-grazers, etc.)

 567  IP=IMAXP1
      IF(IZONE2.EQ.1)THEN
       IF(DV(IP).LT.0.90D0*DVS.AND.XC(IP).GT.PRECX)THEN
        DELTAV=DVS-DV(IP)
        WRITE(2,*) 'GRIDEF call. NTIME,IP,DV,DV(IMAXP1),DELTAV = ',
     1              NTIME,IP,DV(IP),DVS,DELTAV
        WRITE(17,*)NTIME
        WRITE(17,*) 'GRIDEF call. IP,DV,DV(IMAXP1),DELTAV = ',
     1              IP,DV(IP),DVS,DELTAV
        CALL GRIDEF(DELTAV)
        NMERGE=1
       END IF
      ELSE
       IF(DM(IP).LT.0.4D0*DM(IP-1))THEN
        WRITE(2,*) 'SMERGE call. NTIME,IP,DM = ',NTIME,IP,DM(IP)
        WRITE(17,*)NTIME
        WRITE(17,*) 'SMERGE call. IP,DM = ',IP,DM(IP)
        CALL SMERGE(IP)
        NMERGE=1
       END IF
      END IF

!  Calculate collapse (PO>PORCOL) or expansion (PO<POREXP).
!  These values are set in the COMIN file.

      IF(IHYDRO.EQ.0)THEN

       VI=ZERO
       DO I=2,IP
        POI=PO(I)
        FACTOR=ONE
        IF(POI.GT.PORCOL)THEN
         PO(I)=POI/(TWO-POI)
         FACTOR=TWO/(TWO-POI)
        END IF
        IF(POI.LT.POREXP)THEN
         PO(I)=TWO*POI/(ONE+POI)
         FACTOR=ONE/(ONE+POI)
        END IF
        IF(FACTOR.NE.ONE)THEN
         DV(I)=DV(I)/FACTOR
         RHO(I)=FACTOR*RHO(I)
         ROICE(I)=FACTOR*ROICE(I)
         ROA(I)=FACTOR*ROA(I)
         ROB(I)=FACTOR*ROB(I)
         ROC(I)=FACTOR*ROC(I)
         ROV(I)=FACTOR*ROV(I)
         ROD(I)=FACTOR*ROD(I)
         IF(MMAX.GT.0)THEN
          DO M=1,MMAX
           ROG(I,M)=FACTOR*ROG(I,M)
           ROSG(I,M)=FACTOR*ROSG(I,M)
          END DO
         END IF
        END IF
        VI=VI+DV(I)
        R(I)=(VI/PI43)**(ONE/THRE)
        FACR(I)=PI4*R(I)**2
        FSQ(I)=FACR(I)*FACR(I)
       END DO

      ELSE

!  Calculate hydrostatic structure (every NPRINT time steps), if IHYDRO > 0 

       IF(NTIME.EQ.1.OR.MOD(NTIME,NPRINT).EQ.0)THEN
        CALL HYDRO(PSOL,RHO,P,DM,DV,IMAXP1,G,PI4,EPSX2,EPSHY,NTIME)
        VI=ZERO
        DO I=2,IMAXP1
         VI=VI+DV(I)
         R(I)=(VI/PI43)**(ONE/THRE)
         FACR(I)=PI4*R(I)**2
         FSQ(I)=FACR(I)*FACR(I)
         ROA(I)=XA(I)*RHO(I) 
         ROB(I)=XB(I)*RHO(I)
         ROC(I)=XC(I)*RHO(I) 
         ROD(I)=XD(I)*RHO(I) 
         ROV(I)=XV(I)*RHO(I) 
         ROICE(I)=ROA(I)+ROB(I)+ROC(I)
         IF(MMAX.GT.0)THEN
          DO M=1,MMAX
           ROG(I,M)=XG(I,M)*RHO(I)
           ROSG(I,M)=XSG(I,M)*RHO(I)
           ROICE(I)=ROICE(I)+ROSG(I,M)
          END DO
         END IF
        END DO
       END IF

      END IF

      DO I=2,IMAXP1

       T0(I)=T(I)
       U0(I)=U(I)
       DM0(I)=DM(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         ROSGI(M)=ROSG0(I,M)
        END DO
       END IF
       POI0=POR(ROA0(I),ROB0(I),ROC0(I),ROD0(I),RHOICE,RHODUS,
     1          MMAX,ROSGI,RHOSOL)
       RHO0(I)=RHO(I)
       ROICE0(I)=ROICE(I)
       ROA0(I)=ROA(I)
       ROB0(I)=ROB(I)
       ROC0(I)=ROC(I)
       ROV0(I)=ROV(I)
       ROD0(I)=ROD(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         ROG0(I,M)=ROG(I,M)
         ROSG0(I,M)=ROSG(I,M)
         ROSGI(M)=ROSG(I,M)
        END DO
       END IF
       POI=POR(ROA(I),ROB(I),ROC(I),ROD(I),RHOICE,RHODUS,
     1         MMAX,ROSGI,RHOSOL)
       SVCORR=ONE
       ALFPS=ALFP(I)
       IF(ALFPS.GT.ZERO)THEN
        SVCORR=(ALFPS-THRE)**2/((ALFPS-TWO)*(ALFPS-FOUR))
     1         *(XPS**(ALFPS-TWO)-ONE)*(XPS**(ALFPS-FOUR)-ONE)
     2         /(XPS**(ALFPS-THRE)-ONE)**2
       END IF

!       TET(I)=TET(I)*SQRT(POI/POI0)
!       TET(I)=TET(I)*(ONE+(POI-POI0)/(TWO*POI0*SVCORR))
       TET(I)=TET(I)*EXP((POI-POI0)/(TWO*POI0*SVCORR))

       IF(LMAX.GT.0)THEN
        DO L=1,LMAX
         XMRAD0(I,L)=XMRAD(I,L)
        END DO
       END IF

      END DO

      FMH2O0=FMH2O
      FMDUS0=FMDUS
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        FMG0(M)=FMG(M)
       END DO
      END IF
      UTOTP=UTOT
      IF(NTIME.EQ.NTIMAX.OR.MOD(NTIME,NSAVE).EQ.0)THEN
       CALL MODEL
      END IF
      IF(TIME.GT.TIMEFIN.AND.TIMEFIN.GT.0.0D0)THEN
       CALL MODEL
       WRITE(2,*) 'End of Run (TIME > TIMEFIN) --> t = ',TIME
       WRITE(17,*)NTIME
       WRITE(17,*) 'End of Run (TIME > TIMEFIN) --> t = ',TIME
       STOP
      END IF
      IF(NORB.GT.ORBFIN.AND.ORBFIN.GT.0)THEN
       CALL MODEL
       WRITE(2,*) 'End of Run (NORB > ORBFIN) --> Orbit = ',NORB
       WRITE(17,*)NTIME
       WRITE(17,*) 'End of Run (NORB > ORBFIN) --> Orbit = ',NORB
       STOP
      END IF
!      IF(JMAX.GT.0.AND.IB.GT.0)THEN
!       CALL SURFL(0.D0,0.D0,JMAX,EMIS,SIG,SFLUX,PI,0.D0,0.D0,0.D0,
!     1             50  ,HSUB,QTOT,PZERO,EXPT,AIB,TETAB,0.D0,
!     2             0.D0,0.D0,NTIME,-1    ,IB,EPSX1,EPSX2,
!     3             PQV,PFL,PXMFL,PAT,VAFLIP,DVFIPL,FSUR)
!      END IF
      
      RETURN
      END
