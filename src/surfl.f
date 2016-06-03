!*****************************************************************
      SUBROUTINE SURFL(TIP,PP,JM0 ,EMIS,SIG,SFLUX,PI,TNS,XIN,DXIN,
     1                KMAX,HSUB,QTOT,PZERO,EXPT,AIB,TETA,DRMAX,
     2                FIN,DFIN,NTIME,NPRINT,IB,EPSX1,EPSX2,
     3                PQV,PFL,PXMFL,PAT,VAFLIP,DVFIPL,FSUR)
!*****************************************************************

!       Subroutine description:
!       -----------------------
!       Calculate fluxes (heat and mass) through a thin boundary layer, 
!       assuming steady state and plane parallel geometry.
!       Input: pressure and temperature (PP, TIP) (at contact surface).
!       Output: heat and mass fluxes (FIN, DFIN, XIN, DXIN).

!       Subroutine mapping:
!       -------------------
!       Calls: PRINT, TABLE.
!       Called by: PRINT, ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      PARAMETER (IS=64)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION  DR(IS),P(IS),T(IS),FL(IS),DFLL(IS),DFLR(IS)
      DIMENSION  A(2,2),B(2,2),C(2,2),D(2),XT(2,IS),XTP(2,IS)
      DIMENSION  H(2,2),E(2,2),EIM(2,2),GW(2),G(2,2,IS),W(2,IS)
      DIMENSION  QV(IS),DQVDT(IS),DQVDP(IS),XMFL(IS)
      DIMENSION  DXMDTL(IS),DXMDTR(IS),DXMDPL(IS),DXMDPR(IS)
      DIMENSION  PRIMAT(IMX,LR),AM(LR),INDN(LR)
      SAVE
      DATA  KIN/0/

      IF(KIN.EQ.0)THEN
       IF(TETA.EQ.ZERO)THEN
        JMAX=0
        IB=0
        RETURN
       END IF
       QV(1)=ZERO
       DQVDT(1)=ZERO
       DQVDP(1)=ZERO
       XT(1,1)=ZERO
       XT(2,1)=ZERO
       XTP(1,1)=ZERO
       XTP(2,1)=ZERO
       JMAX=0
       KIN=1
      END IF

      IF(NPRINT.LE.0)THEN
       IF(NPRINT.LT.0)THEN
        DO J=2,JMAX
!        WRITE(8,*) DR(J),T(J),P(J)
        END DO
        RETURN
       ELSE
        IF(KIN.GT.1.AND.IB.GT.0)THEN
         WRITE(2,52)
  52     FORMAT(/6H  SURF,12H    DRI     ,12H TEMPERATURE,12H   P-VAP
     1       ,12H   PRESSURE ,12H    FLUX    ,12H   VAP FLUX ,
     2        12H   EDOT     ,12H  VELOCITY  ,12H   V-SOUND  /)
         DO J=1,JP1
          PRIMAT(J,1)=DR(J)
          PRIMAT(J,2)=T(J)
          PRIMAT(J,3)=ZERO
          IF(T(J).GT.ZERO) PRIMAT(J,3)=PZERO*EXP(EXPT/T(J))
          PRIMAT(J,4)=P(J)
          PRIMAT(J,5)=FL(J)
          PRIMAT(J,6)=XMFL(J)
          PRIMAT(J,7)=QV(J)
          PRIMAT(J,8)=ZERO
          IF(XMFL(J).GT.ZERO.AND.P(J+1).GT.ZERO)
     1       PRIMAT(J,8)=XMFL(J)*T(J+1)/P(J+1)/AIB
          IF(XMFL(J).LE.ZERO.AND.P(J).GT.ZERO)
     1       PRIMAT(J,8)=XMFL(J)*T(J)/P(J)/AIB
          PRIMAT(J,9)=SQRT(4.D0/THRE*T(J)/AIB)
         END DO
         KD=-9
         IJ=1
         CALL TABLE(KD,1,JP1,IJ,PRIMAT,AM,INDN)
         WRITE(2,53)FL(1),FL(JP1),XMFL(1),XMFL(JP1),PREC
  53     FORMAT(1X,'BOUNDARY CONVERGENCE ACCURACY',1P4E11.3,
     1       10X,'PREC=',E12.4)
        END IF
        RETURN
       END IF
      END IF

!  Start iterations for T(J) and P(J)

      NDEX=0
      KIN=KIN+1
      PREC=ZERO
      T(1)=TIP
      P(1)=PP
      JMPREV=JMAX
!      JM=LOG(TWO*DRMAX/TETA+ONE)/LOG(TWO)
      JM=LOG(TWO*DRMAX/(1.D1*TETA)+ONE)/LOG(TWO)
!      JM=LOG(DRMAX/TETA+ONE)/LOG(TWO)
      IF(JM.LT.JM0)THEN
       JMAX=JM
      ELSE
       JMAX=JM0
      END IF

      JP1=JMAX+1
!      DR(JP1)=TETA/TWO
      DR(JP1)=1.D1*TETA/TWO
!      DR(JP1)=TETA
      DR(1)=ZERO
      DO J=JMAX,2,-1
       DR(J)=TWO*DR(J+1)
      END DO
      DO J=2,JP1
       XTP(1,J)=ZERO
       XTP(2,J)=ZERO
       XT(1,J)=ZERO
       XT(2,J)=ZERO
       IF(JMAX.NE.JMPREV)THEN
        T(J)=TIP
        P(J)=ZERO
       END IF
      END DO

      NIT=0
      ACON=ZERO
  3   NIT=NIT+1
      ACONP=ACON
      DO J=1,2
       W(J,1)=ZERO
       DO  I=1,2
        H(I,J)=ZERO
       END DO
      END DO

      IQV=0
      TJP=T(JP1)
      PJP=P(JP1)
      DRJP=DR(JP1)*0.5D0
      FL(JP1)=SFLUX-EMIS*SIG*TJP**4
      DFLR(JP1)=ZERO
      DFLL(JP1)=-4.D0*EMIS*SIG*TJP**3
      XQV=PQV*SQRT(AIB/(TWO*PI*TJP))
      DXQVDT=-0.5D0*XQV/TJP
      QV(JP1)=XQV*(PZERO*EXP(EXPT/TJP)-PJP)
      IF(QV(JP1).NE.ZERO) IQV=IQV+1
      DQVDT(JP1)=PZERO*EXP(EXPT/TJP)*(DXQVDT-XQV*EXPT/TJP**2)-DXQVDT*PJP
      DQVDP(JP1)=-XQV
      TERMP=PJP/SQRT(TJP)
      XMFL(JP1)=-PXMFL*TERMP/DRJP
      DTRMPT=-0.5D0*TERMP/TJP
      DTRMPP=ONE/SQRT(TJP)
      DXMDTR(JP1)=ZERO
      DXMDPR(JP1)=ZERO
      DXMDTL(JP1)=-0.5D0*XMFL(JP1)/TJP
      DXMDPL(JP1)=-PXMFL/DRJP/SQRT(TJP)

      DO J=JMAX,2,-1

       DRJ=DR(J)*0.5D0
       DRPDR=DRJP+DRJ
       TJ=T(J)
       PJ=P(J)

       XFL=PFL/DRPDR
       FL(J)=XFL*(TJP-TJ)
       DFLL(J)=-XFL
       DFLR(J)=XFL

!  Change introduced for ISSI calculations:

!       XFL=1.D-2*5.67D7*(DRJ/TJ+DRJP/TJP)/DRPDR**2
!       DXFLL=-1.D-2*5.67D7*DRJ/TJ**2/DRPDR**2
!       DXFLR=-1.D-2*5.67D7*DRJP/TJP**2/DRPDR**2
!       FL(J)=XFL*(TJP-TJ)
!       DFLL(J)=-XFL+DXFLL*(TJP-TJ)
!       DFLR(J)=XFL+DXFLR*(TJP-TJ)

       XQV=PQV*SQRT(AIB/(TWO*PI*TJ))
       DXQVDT=-0.5D0*XQV/TJ
       QV(J)=XQV*(PZERO*EXP(EXPT/TJ)-PJ)
       IF(QV(J).NE.ZERO) IQV=IQV+1
       DQVDT(J)=QV(J)/XQV*DXQVDT-XQV*PZERO*EXP(EXPT/TJ)*EXPT/TJ**2
       DQVDT(J)=PZERO*EXP(EXPT/TJ)*(DXQVDT-XQV*EXPT/TJ**2)-DXQVDT*PJ
       DQVDP(J)=-XQV
       XFL=PXMFL/DRPDR
       TERM=PJ/SQRT(TJ)
       DTRMT=-0.5D0*TERM/TJ
       DTRMP=ONE/SQRT(TJ)
       XMFL(J)=XFL*(TERMP-TERM)
       DXMDTL(J)=-XFL*DTRMT
       DXMDTR(J)=XFL*DTRMPT
       DXMDPL(J)=-XFL*DTRMP
       DXMDPR(J)=XFL*DTRMPP
       TJP=TJ
       PJP=PJ
       DRJP=DRJ
       TERMP=TERM
       DTRMPT=DTRMT
       DTRMPP=DTRMP

      END DO

      IF(IQV.EQ.0)THEN
       FIN=FL(JP1)
       DFIN=DFLL(JP1)
       FSUR=FL(JP1)
       QTOT=ZERO
       XIN=ZERO
       DXIN=ZERO
       RETURN
      END IF
      DRJ=0.5D0*DR(2)

!  Change introduced for ISSI calculations:

!      PFL=1.D-2*5.67D7/T(2)
!      DFLR(1)=PFL/DRJ*(ONE-(T(2)-TIP)/T(2))

      FL(1)=PFL/DRJ*(T(2)-TIP)
      FIN=FL(1)
      DFIN=-PFL/DRJ
      DFLL(1)=ZERO
      DFLR(1)=PFL/DRJ
      XMFL(1)=PXMFL*(P(2)/SQRT(T(2))-PP/SQRT(TIP))/DRJ
      XIN=XMFL(1)
      DXIN=0.5D0*PXMFL*PP/(TIP*SQRT(TIP))/DRJ
      DXMDTL(1)=ZERO
      DXMDPL(1)=ZERO
      DXMDTR(1)=-0.5D0*PXMFL*P(2)/(T(2)*SQRT(T(2)))/DRJ
      DXMDPR(1)=PXMFL/SQRT(T(2))/DRJ
      ACON=ZERO

      DO I=2,JP1

       DRI=DR(I)
       A(1,1)=DFLR(I)
       B(1,1)=-DFLL(I)+DFLR(I-1)+DRI*HSUB*DQVDT(I)
       C(1,1)=-DFLL(I-1)
       A(1,2)=ZERO
       B(1,2)=DRI*HSUB*DQVDP(I)
       C(1,2)=ZERO
       D(1)=FL(I)-FL(I-1)-DRI*HSUB*QV(I)
       A(2,1)=DXMDTR(I)
       B(2,1)=-DXMDTL(I)+DXMDTR(I-1)-DRI*DQVDT(I)
       C(2,1)=-DXMDTL(I-1)
       A(2,2)=DXMDPR(I)
       B(2,2)=-DXMDPL(I)+DXMDPR(I-1)-DRI*DQVDP(I)
       C(2,2)=-DXMDPL(I-1)
       D(2)=XMFL(I)-XMFL(I-1)+DRI*QV(I)
       ACON=ACON+SQRT(D(1)**2+D(2)**2)

       DO L=1,2
        DO J=1,2
         E(L,J)=B(L,J)
         DO K=1,2
          E(L,J)=E(L,J)-C(L,K)*H(K,J)
         END DO
        END DO
       END DO

       DET=E(1,1)*E(2,2)-E(1,2)*E(2,1)
       IF(DET.EQ.ZERO)THEN 
        WRITE(2,*) I,J,K,L,NIT,DET,ACON,ACONP
        WRITE(17,*)NTIME
        WRITE(17,*)I,J,K,L,NIT,DET,ACON,ACONP
       END IF
       EIM(1,1)=E(2,2)/DET
       EIM(1,2)=-E(1,2)/DET
       EIM(2,1)=-E(2,1)/DET
       EIM(2,2)=E(1,1)/DET

       DO L=1,2
        GW(L)=D(L)
        DO K=1,2
         GW(L)=GW(L)+C(L,K)*W(K,I-1)
        END DO
       END DO

       DO L=1,2
        DO J=1,2
         H(L,J)=ZERO
         DO  K=1,2
          H(L,J)=H(L,J)+EIM(L,K)*A(K,J)
         END DO
        END DO
       END DO

       DO L=1,2
        W(L,I)=ZERO
        DO K=1,2
         W(L,I)=W(L,I)+EIM(L,K)*GW(K)
         G(L,K,I)=H(L,K)
        END DO
       END DO

      END DO

      XT(1,JP1)=W(1,JP1)
      XT(2,JP1)=W(2,JP1)
      DO I=JMAX,2,-1
       DO L=1,2
        XT(L,I)=W(L,I)
        DO K=1,2
         XT(L,I)=XT(L,I)+G(L,K,I)*XT(K,I+1)
        END DO
       END DO
      END DO

!  Check convergence

      KCONV=0
      DTMAX=ZERO
      DPMAX=ZERO

      DO J=2,JP1

       DELT=XT(1,J)
       DELX1=ABS(DELT)
       DELP=XT(2,J)
       DELX2=ABS(DELP)
       TJ=T(J)+DELT
       IF(TJ.LE.ZERO) TJ=T(J)*EXP(DELT/T(J))
       T(J)=TJ
       IF(DELX1.GT.EPSX1*TJ)THEN
        KCONV=KCONV+1
        NDEX=2
       END IF
       PJ=P(J)+DELP
       IF(PJ.LT.ZERO.AND.P(J).EQ.ZERO) PJ=ZERO
       IF(PJ.LT.ZERO.AND.P(J).GT.ZERO) PJ=P(J)*EXP(DELP/P(J))
       P(J)=PJ
       IF(PJ.GT.ZERO.AND.DELX2.GT.EPSX2*PJ)THEN
        KCONV=KCONV+1
        IF(NDEX.EQ.0)THEN
         NDEX=3
        ELSE
         NDEX=5
        END IF
       END IF
       XTP(1,J)=XTP(1,J)+DELT
       XTP(2,J)=XTP(2,J)+DELP
       IF(TJ.GT.ZERO) DTMAX=MAX(DTMAX,ABS(XTP(1,J)/TJ))
       IF(PJ.GT.ZERO) DPMAX=MAX(DPMAX,ABS(XTP(2,J)/PJ))
      END DO

      IF(KCONV.GT.0)THEN
       IF(NIT.LT.KMAX) GO TO 3
       IF(NIT.LT.KMAX+10.AND.ACON.LT.ACONP) GO TO 3
       WRITE(2,50)NDEX,KCONV,NIT,ACON,ACONP,TIP,PP,T(JP1),P(JP1),
     1            FL(JP1),FL(1),DTMAX,DPMAX
       WRITE(17,*)NTIME
       WRITE(17,50)NDEX,KCONV,NIT,ACON,ACONP,TIP,PP,T(JP1),P(JP1),
     1             FL(JP1),FL(1),DTMAX,DPMAX
  50   FORMAT(1X,' NO SOLUTION FOR SURFACE LAYER',3I5/1P10E13.5)
       CALL PRINT
       STOP
      ELSE
       QTOT=ZERO
       DO J=2,JP1
        QTOT=QTOT+QV(J)*DR(J)
       END DO
       JS=ONE/LOG(TWO)+ONE
!       TNS=T(JP1-JS)
!       FSUR=FL(JP1-JS)
       TNS=T(JP1-1)
       FSUR=FL(JP1-1)
       PREC=(FL(1)-FL(JP1)+QTOT*HSUB)/FL(1)
      END IF

      RETURN
      END
