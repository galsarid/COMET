!**************************
      SUBROUTINE SMERGE(II)
!**************************

!       Subroutine description:
!       -----------------------
!       Merge two shells into one. We do it in case one shell is much
!       Thinner than the next one.

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      SAVE

      IP=II
      I=IP-1
      DMI=DM(I)
      DMIP=DM(IP)
      DMPDM=DMI+DMIP
      DVI=DV(I)
      DVIP=DV(IP)
      DVPDV=DVI+DVIP
      RHO(I)=DMPDM/DVPDV
      U(I)=(U(I)*DVI+U(IP)*DVIP)/DVPDV
      ROICE(I)=(ROICE(I)*DVI+ROICE(IP)*DVIP)/DVPDV
      ROAI=ROA(I)*DVI+ROA(IP)*DVIP
      ROB(I)=(ROB(I)*DVI+ROB(IP)*DVIP)/DVPDV
      ROC(I)=(ROC(I)*DVI+ROC(IP)*DVIP)/DVPDV
      ROV(I)=(ROV(I)*DVI+ROV(IP)*DVIP)/DVPDV
      ROD(I)=(ROD(I)*DVI+ROD(IP)*DVIP)/DVPDV
      UCON=ROICE(I)*BETA+ROV(I)*CV+ROD(I)*CD
      SUMX=ROV(I)/AIB
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        ROSG(I,M)=(ROSG(I,M)*DVI+ROSG(IP,M)*DVIP)/DVPDV
        ROG(I,M)=(ROG(I,M)*DVI+ROG(IP,M)*DVIP)/DVPDV
        ROSGI(M)=ROSG(I,M)
        UCON=UCON+ROG(I,M)*CG(M)
        SUMX=SUMX+ROG(I,M)/AGB(M)
        IF(ROAI.GT.0.)THEN
         ZGI(I,M)=(ZGI(I,M)*ROA(I)*DVI+ZGI(IP,M)*ROA(IP)*DVIP)/ROAI
        ELSE
         ZGI(I,M)=ZERO
        END IF
       END DO
      END IF
      ROA(I)=ROAI/DVPDV
      TET(I)=(TET(I)*DMI+TET(IP)*DMIP)/DMPDM
      TI=T(I)
      UI=U(I)
      IF(TI.GT.TML)THEN
       T(I)=RTBIS(TI,UI,ROICE(I),UCON,ALFA,BETA,EXMELT,HMELT,TMELT,CW)
      ELSE
       IF(ROICE(I).GT.EPSX2)THEN
        ALF=ROICE(I)*ALFA
        T(I)=(SQRT(UCON*UCON+4.D0*ALF*UI)-UCON)/(TWO*ALF)
       ELSE
        T(I)=UI/UCON
       END IF
      END IF

      VI=PI43*R(I-1)**3+DVPDV
      R(I)=(VI/PI43)**(ONE/THRE)
      PO(I)=POR(ROA(I),ROB(I),ROC(I),ROD(I),RHOICE,RHODUS,
     1          MMAX,ROSGI,RHOSOL)
      P(I)=T(I)*SUMX/PO(I)
      DM(I)=DMPDM
      DV(I)=DVPDV
      IF(LMAX.GT.0)THEN
       DO L=1,LMAX
        XMRAD(I,L)=XMRAD(I,L)+XMRAD(IP,L)
       END DO
      END IF

      WRITE(2,*)'MERGE SHELLS: NTIME,I-1,I,IMAXP1,DM,DV,U,T,ROICE,P.'
      WRITE(2,11)NTIME,I,IP,IMAXP1,DM(I),DV(I),U(I),T(I),ROICE(I),P(I)
      WRITE(17,*)NTIME
      WRITE(17,*)'MERGE SHELLS: NTIME,I-1,I,IMAXP1,DM,DV,U,T,ROICE,P.'
      WRITE(17,11)NTIME,I,IP,IMAXP1,DM(I),DV(I),U(I),T(I),ROICE(I),P(I)
  11  FORMAT(1X,'MERGE SHELLS',I8,3I5,1P6E12.5)

      IF(IP.LT.IMAXP1)THEN
       IL=IMAXP1-1
       DO I=IP,IL
        INN=I+1
        U(I)=U(INN)
        ROICE(I)=ROICE(INN)
        ROA(I)=ROA(INN)
        ROB(I)=ROB(INN)
        ROC(I)=ROC(INN)
        ROD(I)=ROD(INN)
        ROV(I)=ROV(INN)
        TET(I)=TET(INN)
        T(I)=T(INN)
        RHO(I)=RHO(INN)
        PO(I)=PO(INN)
        P(I)=P(INN)
        DM(I)=DM(INN)
        DV(I)=DV(INN)
        VI=VI+DV(I)
        R(I)=(VI/PI43)**(ONE/THRE)
        FACR(I)=PI4*R(I)**2
        FSQ(I)=FACR(I)*FACR(I)
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          ZGI(I,M)=ZGI(INN,M)
          ROSG(I,M)=ROSG(INN,M)
          ROG(I,M)=ROG(INN,M)
         END DO
        END IF
        IF(LMAX.GT.0)THEN
         DO L=1,LMAX
          XMRAD(I,L)=XMRAD(INN,L)
         END DO
        END IF
       END DO
      END IF

      IP=IMAXP1
      U(IP)=ZERO
      ROICE(IP)=ZERO
      ROA(IP)=ZERO
      ROB(IP)=ZERO
      ROC(IP)=ZERO
      ROD(IP)=ZERO
      ROV(IP)=ZERO
      TET(IP)=ZERO
      T(IP)=ZERO
      RHO(IP)=ZERO
      PO(IP)=ZERO
      P(IP)=ZERO
      DM(IP)=ZERO
      DV(IP)=ZERO
      R(IP)=ZERO
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        ZGI(IP,M)=ZERO
        ROG(IP,M)=ZERO
        ROSG(IP,M)=ZERO
       END DO
      END IF
      IF(LMAX.GT.0)THEN
       DO L=1,LMAX
        XMRAD(IP,L)=ZERO
       END DO
      END IF 
      IMAXP1=IMAXP1-1

      RETURN
      END
