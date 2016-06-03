!*********************
      SUBROUTINE INPUT
!*********************

!       Subroutine description:
!       -----------------------
!       This subroutine initializes arrays and matrices, reads the data
!       file (comin) and builds an initial model if NMODEL is negative.
!       Otherwise - reads a model from the COMODI file.
!       After the model is built - prints data to output files.

!       Subroutine mapping:
!       -------------------
!       Calls: PRINT.
!       Called by: Main.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION XRADIN(LR),XSGIN(LG)
      SAVE


!  Initialize all physical variables of the nucleus model

      DO I=1,IMX
       ALFP(I)=ZERO
       DM(I)=ZERO
       DM0(I)=ZERO
       DV(I)=ZERO
       FL(I)=ZERO
       XMFL(I)=ZERO
       RHO(I)=ZERO
       RHO0(I)=ZERO
       R(I)=ZERO
       FACR(I)=ZERO
       FSQ(I)=ZERO
       P(I)=ZERO
       PSOL(I)=ZERO
       QV(I)=ZERO
       S(I)=ZERO
       T(I)=ZERO
       T0(I)=ZERO
       U(I)=ZERO
       U0(I)=ZERO
       XA(I)=ZERO
       ROA0(I)=ZERO
       ROA(I)=ZERO
       XA(I)=ZERO
       ROA0(I)=ZERO
       ROA(I)=ZERO
       XB(I)=ZERO
       ROB0(I)=ZERO
       ROB(I)=ZERO
       XC(I)=ZERO
       ROC0(I)=ZERO
       ROC(I)=ZERO
       XV(I)=ZERO
       ROV0(I)=ZERO
       ROV(I)=ZERO
       XD(I)=ZERO
       ROD0(I)=ZERO
       ROD(I)=ZERO
       ROICE(I)=ZERO
       ROICE0(I)=ZERO
       XTP(I)=ZERO
       TET(I)=ZERO
       ZKN(I)=ONE
       DO M=1,LG
        GAFL(I,M)=ZERO
        QVG(I,M)=ZERO
        QG(I,M)=ZERO
        XG(I,M)=ZERO
        ROG0(I,M)=ZERO
        ROG(I,M)=ZERO
        XSG(I,M)=ZERO
        ROSG0(I,M)=ZERO
        ROSG(I,M)=ZERO
        ZGI(I,M)=ZERO
       END DO
       DO L=1,LR
        XMRAD(I,L)=ZERO
        XRAD(I,L)=ZERO
        XMRAD0(I,L)=ZERO
       END DO
      END DO
      DO M=1,LG
       ROSGI(M)=ZERO
      END DO
      ACOVAM=ZERO
      DFIMOD=ZERO
      TNS=ZERO
      IAC=1
      
!  Read data from COMIN file
 
 800  FORMAT(BZ,10X,12A4,I7)
 801  FORMAT(BZ,15X,I10)
 802  FORMAT(BZ,15X,F20.1)
 803  FORMAT(A36)
      
      READ(1,800) JOBNAM,IJOB
      READ(1,803) JUNK
      READ(1,801) NMODEL,RUNM
      READ(1,803) JUNK
!  Read-in output parameters:
      READ(1,801) HANGN
      READ(1,801) NPRINT,NSAVE,NSTEP,ITABLE
      READ(1,803) JUNK
!  Read-in flags:
      READ(1,801) KP,IOP,ICORR
      READ(1,801) IHYDRO,EOSF
      READ(1,801) IZONE1,IZONE2
      READ(1,801) IIMPAC
      READ(1,803) JUNK
!  Read-in temporal parameters:
      READ(1,801) NTIMAX
      READ(1,802) TIME0,TIMEFIN
      READ(1,801) ORBFIN
      READ(1,801) ITMODE,NTDAY
      READ(1,803) JUNK
!  Read-in numerics parameters:
      READ(1,801) KMAX,KN
      READ(1,802) DFI,DFIMAX,DTMAX
      READ(1,802) EPSTIM,EPSX1,EPSX2,EPSHY
      READ(1,802) FAC,FACTIM
      READ(1,802) PAT,PRECX
      READ(1,801) IMAXP1
      READ(1,801) ISL,JMAX
      READ(1,802) DELR,QDM
      READ(1,803) JUNK
!  Read-in physics parameters:
      READ(1,802) AXIS,ECCEN
      READ(1,802) ALB,EMIS
      READ(1,802) PROT
      READ(1,802) ALAT,DECLIN,JANGLE
      READ(1,802) RADIUS,RO
      READ(1,802) TT
      READ(1,802) TETA,TETAB
      READ(1,802) ALFPS,XPS
      READ(1,802) XPCR,TZERO
      READ(1,802) ADMAX,BDUST
      READ(1,802) RO0FAC
      READ(1,802) EFRAD,HERTZF
      READ(1,802) PORCOL,POREXP
      READ(1,803) JUNK
!  Read-in composition parameters:
      READ(1,801) MMAX
      READ(1,802) (ZG(M),M=1,LG)
      READ(1,802) XAI,XCI,XDI
      READ(1,802) ZAC,CLIMIT
      READ(1,802) (XSGIN(M),M=1,LG)
      READ(1,801) LMAX 
      READ(1,802) (XRADIN(L),L=1,LR)
      READ(1,803) JUNK
!  Read-in optional parameters:
      READ(1,803) JUNK
!  Read-in impact-related parameters:
      READ(1,801) NDI
      READ(1,802) DDI
      READ(1,802) THICK,CRDIAM,DELTDI
      READ(1,802) KENRAT,ENPART
      READ(1,803) JUNK

      PSICON=(AXIS**(THRE/TWO))/SQRT(GM)
      ALAT=ALAT*PI/180.0D0
      DECLIN=DECLIN*PI/180.0D0
      JANGLE=JANGLE*PI/180.0D0

      IF(NMODEL.LT.0)THEN

!  Build initial model (isothermal, homogeneous) for given R,T,RO,XA,XC
!  and mass shells as a geometric series increasing inwards (with factor QDM).
!  Below (IMAXP1-100) mass shells are of equal thickness - fine zoned.

       NMOD1=0

!  Normally, an evolution run starts at aphelion (FI=PI), otherwise the
!  initial FI must be set to a different value 

       TIME=ZERO
       FI=PI

       RADCOR=RADIUS
       TOTVOL=PI43*RADIUS**3
       CORVOL=TOTVOL

!  If IMAXP1 > ISL then divide the nucleus in two parts:
!  The outer shells are of constant dR and the inner shells are a geometric 
!  series in dV. Otherwise, divide the comet only into shells that are a 
!  geometric series in dV.
!  Note: DV(I) holds the volume between shell number I and shell number I-1.

       IF(IMAXP1.GT.ISL)THEN
        DO I=IMAXP1,ISL+1,-1
         RADCOR=RADCOR-DELR
         CORV=PI43*RADCOR**3
         DV(I)=CORVOL-CORV
         CORVOL=CORV
        END DO
        IP=ISL
       ELSE
        IP=IMAXP1
       END IF
       DV(IP)=CORVOL*(QDM-ONE)/(QDM**(IP-1)-ONE)
       DO I=IP-1,2,-1
        DV(I)=QDM*DV(I+1)
       END DO

!  For the entire comet, calculate the dM array and initialize the
!  parameters of each shell.

       DO I=IMAXP1,2,-1

        DM(I)=DV(I)*RO
        T(I)=TT
        RHO(I)=RO
        TET(I)=ABS(TETA)
        XD(I)=XDI
 
        IF(ZAC.LT.ONE)THEN
!  ZAC = non-crystallizing fraction of amorphous ice
         ICLAT=0
         XA(I)=XAI*(ONE-ZAC)
         XB(I)=XAI*ZAC
         XC(I)=XCI
        ELSE
!  ZAC = clathrate (crystalline) ice mass fraction	
         ICLAT=1
         XB(I)=XAI+XCI
         XA(I)=ZERO
         XC(I)=ZERO
        END IF

!  Define trapped gas species in amorphous ice.
        IF(MMAX.GT.0)THEN
         ZGAS=ZERO
         DO M=1,MMAX
          XSG(I,M)=XSGIN(M)
          ZGI(I,M)=ZG(M)
          ZGAS=ZGAS+ZG(M)
         END DO
         IF(ZAC.EQ.ONE.AND.ZGAS.GT.CLIMIT)THEN
          WRITE(2,*)'ERROR: too much gas in clathrate',I,ZGAS,CLIMIT
          STOP
         END IF
        END IF

        IF(LMAX.GT.0)THEN
         DO L=1,LMAX
          XMRAD(I,L)=DM(I)*XRADIN(L)
         END DO
        END IF

       END DO

       PCRIT=10.0D0**XPCR
       DO I=1,IMAXP1
        PCR(I)=PCRIT
       END DO

       PSICON=(AXIS**(THRE/TWO))/SQRT(GM)
       IF(ECCEN.GT.ZERO) TIME=PSICON*(FI-ECCEN*SIN(FI)-PI)

      ELSE

!  Read model from COMODI file
!  The file contains parameters for each shell.

       READ(7,800)JOBNAM,IJOB
       READ(7,*)NMOD1
       READ(7,*)IMAXP1,JMAX,MMAX,LMAX,ICLAT
       READ(7,*)TIME
       TIME_OLD=TIME
       READ(7,*)FI,DFIMOD
       READ(7,*)AXIS_OLD,ECCEN_OLD

!  Check if the orbital elements in COMIN changed. 
!  If so, adjust orbital position through time (a very rough approximation).

       IF((AXIS_OLD.NE.AXIS).OR.(ECCEN_OLD.NE.ECCEN))THEN
        PSICON=(AXIS**(THRE/TWO))/SQRT(GM)
        MA=TIME/PSICON
        FI=MA+PI
        TIME=PSICON*(FI-ECCEN*SIN(FI)-PI)
        WRITE(2,*)'a and e changed from last run:'
        WRITE(2,*)'Old = ',AXIS_OLD/ONEAU,ECCEN_OLD
        WRITE(2,*)'New = ',AXIS/ONEAU,ECCEN
        WRITE(17,*)NTIME
        WRITE(17,*)'a and e changed from last run:'
        WRITE(17,*)'Old = ',AXIS_OLD/ONEAU,ECCEN_OLD
        WRITE(17,*)'New = ',AXIS/ONEAU,ECCEN
        WRITE(*,*)'Old: a,e,fi,t = ',AXIS_OLD/ONEAU,ECCEN_OLD,FI_OLD,
     1                             TIME_OLD
        WRITE(*,*)'New: a,e,fi,t = ',AXIS/ONEAU,ECCEN,FI,TIME
       END IF

       DO I=2,IMAXP1
        READ(7,*)K,DM(I),RHO(I),T(I),TET(I),XA(I),XB(I),XC(I),
     1           (XSG(I,M),M=1,LG),XD(I),(XG(I,M),M=1,LG),XV(I),
     2           (ZGI(I,M),M=1,LG),PCR(I),PSOL(I),ALFP(I)
        IF(LMAX.GT.0)THEN
         READ(7,*)K,(XMRAD(I,L),L=1,LMAX)
        END IF
       END DO

      END IF


      IJOB=IJOB+1

!  Write to COMODO file.
      WRITE(8,800) JOBNAM,IJOB

!  Prepare and print initial parameters

      IB=0
      VI=ZERO
      DO I=2,IMAXP1
       DMI=DM(I)
       DVI=DMI/RHO(I)
       VI=VI+DVI
       R(I)=(VI/PI43)**(ONE/THRE)
       FACR(I)=PI4*R(I)**2
       FSQ(I)=FACR(I)*FACR(I)
       DV(I)=DVI
       SUMX=XA(I)+XB(I)+XC(I)+XV(I)+XD(I)
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         SUMX=SUMX+XSG(I,M)+XG(I,M)
        END DO
       END IF
       IF(ABS(SUMX-ONE).GT.1.D-6)THEN
        WRITE(2,*)'ERROR: Unbalanced composition. I,SUM,XA,XB,XC,XD = ',
     1             I,SUMX,XA(I),XB(I),XC(I),XD(I)
        WRITE(17,*)NTIME
        WRITE(17,*)'ERROR: Unbalanced composition. I,SUM,XA,XB,XC,XD =',
     1              I,SUMX,XA(I),XB(I),XC(I),XD(I)
       END IF
      END DO

!  Spin period may be supplied as fraction of orbital period (negative PROT)
      IF(PROT.LT.ZERO)THEN
       PCYC=(AXIS/ONEAU)**(THRE/TWO)*YEAR
       PROT=-PROT*PCYC
      END IF

      IF(PROT.GT.ZERO)THEN
!  Solar flux (to be multiplied by local zenith angle)
       Q=(ONE-ALB)*SLM
      ELSE
!  Fast-rotator approximation (equally heated surface by averaged flux).
       Q=(ONE-ALB)*SLM*0.25D0
      END IF

      RADIUS=R(IMAXP1)
      UTOT0=ZERO
      FMTOT0=ZERO

!  Run on all shells

      DO I=2,IMAXP1

       ETA(I)=ZERO
       DMI=DM(I)
       DM0(I)=DMI
       FMTOT0=FMTOT0+DMI
       RHO0(I)=RHO(I)
       TI=T(I)
       T0(I)=TI
       ROICE(I)=(XA(I)+XB(I)+XC(I))*RHO(I)
       ROD(I)=XD(I)*RHO(I)
       ROD0(I)=ROD(I)
       ROA(I)=XA(I)*RHO(I)
       ROA0(I)=ROA(I)
       ROB(I)=XB(I)*RHO(I)
       ROB0(I)=ROB(I)
       ROC(I)=XC(I)*RHO(I)
       ROC0(I)=ROC(I)
       ROV(I)=XV(I)*RHO(I)
       ROV0(I)=ROV(I)
       U0(I)=(ROV0(I)*CV+ROD0(I)*CD)*TI
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         ROSG(I,M)=XSG(I,M)*RHO(I)
         ROSG0(I,M)=ROSG(I,M)
         ROICE(I)=ROICE(I)+ROSG(I,M)
         ROG(I,M)=XG(I,M)*RHO(I)
         ROG0(I,M)=ROG(I,M)
         U0(I)=U0(I)+ROG0(I,M)*CG(M)*TI
         ROSGI(M)=ROSG(I,M)
        END DO
       END IF
       ROICE0(I)=ROICE(I)
       PO(I)=POR(ROA(I),ROB(I),ROC(I),ROD(I),RHOICE,RHODUS,
     1           MMAX,ROSGI,RHOSOL)
       IF(TI.LT.TML)THEN
        UICE=ALFA*TI*TI+BETA*TI
       ELSE
        XM=ONE/(ONE+EXP(EXMELT*(ONE-TI/TMELT)))
        UICE=(ONE-XM)*(ALFA*TI*TI+BETA*TI)+XM*(HMELT+CW*TI)
       END IF
       U0(I)=U0(I)+ROICE0(I)*UICE
       UTOT0=UTOT0+U0(I)*DV(I)
       U(I)=U0(I)
       IF(LMAX.GT.0)THEN
        DO L=1,LMAX
         XMRAD0(I,L)=XMRAD(I,L)
        END DO
       END IF

      END DO

!  Write to COMOUT file

      WRITE(2,810) JOBNAM,IJOB
 810  FORMAT(///1H1,10X,12A4,I7//)
      IF(TIME0.LT.ZERO) TIME0=TIME
      IF(DFIMOD.GT.ZERO) DFI=DFIMOD
      DFI=TWO*PI/DFI
      DFIMAX=TWO*PI/DFIMAX
      PROTD=PROT/(2.4D1*3.6D3)
      DVS=DV(IMAXP1)

      WRITE(2,811)IMAXP1,JMAX,MMAX,LMAX,KMAX,NPRINT,NSAVE,NSTEP,ITABLE,
     1            NMODEL,NTIMAX,ISL,KN,KP
 811  FORMAT(//8H IMAXP1=,I10,8H JMAX  =,I10,8H MMAX  =,I10,
     1         8H LMAX  =,I10,8H KMAX  =,I10/8H NPRINT=,I10,
     2         8H NSAVE =,I10,8H NSTEP =,I10,8H ITABLE=,I10,
     3         8H NMODEL=,I10,8H NTIMAX=,I10,
     4         8H ISL   =,I10,8H KN    =,I10,8H KP    =,I10/)
      WRITE(2,812)TIME0,DFI,DFIMAX,DTMAX,EPSTIM,EPSX1,EPSX2,FACTIM,
     1            ALB,AXIS,ECCEN,EMIS,FMTOT0,PRECX,TETA
 812  FORMAT(//8H TIME0 =,1PE14.6,8H DFI   =,1PE14.6,8H DFIMAX=,1PE14.6,
     1         8H DTMAX =,1PE14.6,8H EPSTIM=,1PE14.6/8H EPSX1 =,1PE14.6,
     2         8H EPSX2 =,1PE14.6,8H FACTIM=,1PE14.6,8H ALB   =,1PE14.6,
     3         8H AXIS  =,1PE14.6/8H ECCEN =,1PE14.6,8H EMIS  =,1PE14.6,
     4         8H FMTOT =,1PE14.6,8H PRECX =,1PE14.6,8H TETA  =,1PE14.6)
      WRITE(2,813)TETAB,PAT,FAC,XPCR,QDM,ZAC,(ZG(M),M=1,LG),TT,RO,
     1            RADIUS,PROTD,TZERO,ADMAX,BDUST
 813  FORMAT(  8H TETAB =,1PE14.6,8H PAT   =,1PE14.6,8H FAC   =,1PE14.6,
     1         8H EXPCR =,1PE14.6,8H QDM   =,1PE14.6/8H ZAC   =,1PE14.6,
     2         8H Z(CO) =,1PE14.6,8H Z(CO2)=,1PE14.6,8H Z(3)  =,1PE14.6,
     3         8H Z(4)  =,1PE14.6/8H Z(5)  =,1PE14.6,8H TINIT =,1PE14.6,
     4         8H ROINIT=,1PE14.6,8H RADIUS=,1PE14.6,8H PROTD =,1PE14.6/
     5         8H TZERO =,1PE14.6,8H ADMAX =,1PE14.6,8H BDUST =,1PE14.6)
      WRITE(2,814)NMOD1,IJOB
 814  FORMAT(//1X,' INPUT MODEL NUMBER=',I5,'  JOB NUMBER=',I5/)
      
      DIST=AXIS*(ONE-ECCEN*COS(FI))
      XORB=0.5D0*(ONE+FI/PI)
      NORB=INT(XORB)

      CALL PRINT

      RETURN
      END
