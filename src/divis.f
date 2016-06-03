!*******************************
      SUBROUTINE DIVIS(MODIV,II)
!*******************************

!       Subroutine description:
!       -----------------------
!       MODIV=0 : Add two adjacent shells and divide the new shell into 2
!                 (nr. of grid points remains unchanged)
!       MODIV=1 : Divide shell in to 2 (nr. of grid points increases by 1)

!       Subroutine mapping:
!       -------------------
!       Calls: None. 
!       Called by: ENDSTP, FLUX.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      SAVE

      IF(MODIV.EQ.0)THEN

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
         IF(ROAI.GT.ZERO)THEN
          ZGI(I,M)=(ZGI(I,M)*ROA(I)*DVI+ZGI(IP,M)*ROA(IP)*DVIP)/ROAI
         ELSE
          ZGI(I,M)=ZERO
         ENDIF
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
     1            MMAX,ROSGI,RHOSOL)
       P(I)=T(I)*SUMX/PO(I)
       DM(I)=DMPDM/TWO
       DV(I)=DVPDV/TWO

      ELSE 
!  For MODIV=1

       IP=I+1
       IMAXP1=IMAXP1+1
       IF(I.LT.IMAXP1)THEN
        DO I=IMAXP1,IP+1,-1
         INN=I-1 
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
         R(I)=R(INN)
         FACR(I)=FACR(INN)
         FSQ(I)=FSQ(INN)
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
        I=IP-1
       END IF 
       DV(I)=DV(I)/TWO
       DM(I)=DM(I)/TWO

      END IF

      VI=PI43*R(I-1)**3+DV(I)
      R(I)=(VI/PI43)**(ONE/THRE)
      FACR(I)=PI4*R(I)**2
      FSQ(I)=FACR(I)*FACR(I)
      IF(LMAX.GT.0)THEN
       DO L=1,LMAX
        IF(MODIV.EQ.0)THEN 
         XMRAD(I,L)=(XMRAD(I,L)+XMRAD(IP,L))/TWO
        ELSE
         XMRAD(I,L)=XMRAD(I,L)/TWO
        END IF
       END DO
      END IF

      WRITE(2,*)'REZONE SHELLS: NTIME,I,I+1,IMAXP1+1,DM,DV,U,T,ROICE,P.'
      WRITE(2,11)NTIME,I,IP,IMAXP1,DM(I),DV(I),U(I),T(I),ROICE(I),P(I)
      WRITE(17,*)NTIME
      WRITE(17,*)'REZONE SHELLS: NTIME,I,I+1,IMAXP1+1,DM,DV,U,T,ROICE,P'
      WRITE(17,11)NTIME,I,IP,IMAXP1,DM(I),DV(I),U(I),T(I),ROICE(I),P(I)
  11  FORMAT(1X,'REZONE SHELLS',I8,3I4,1P6E12.5)

      INN=I
      I=I+1
      U(I)=U(INN)
      ROICE(I)=ROICE(INN)
      ROA(I)=ROA(INN)
      ROB(I)=ROB(INN)
      ROC(I)=ROC(INN)
      ROD(I)=ROD(INN)
      ROV(I)=ROV(INN)
      T(I)=T(INN)
      TET(I)=TET(INN)
      RHO(I)=RHO(INN)
      PO(I)=PO(INN)
      P(I)=P(INN)
      DM(I)=DM(INN)
      DV(I)=DV(INN)
      VI=PI43*R(INN)**3+DV(I)
      R(I)=(VI/PI43)**(ONE/THRE)
      FACR(I)=PI4*R(I)**2
      FSQ(I)=FACR(I)*FACR(I)
      IF(LMAX.GT.0)THEN
       DO L=1,LMAX
        XMRAD(I,L)=XMRAD(INN,L)
       END DO
      END IF
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        ZGI(I,M)=ZGI(INN,M)
        ROSG(I,M)=ROSG(INN,M)
        ROG(I,M)=ROG(INN,M)
       END DO
      END IF

      RETURN
      END
