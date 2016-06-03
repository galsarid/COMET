!*********************
      SUBROUTINE MFLUX
!*********************

!       Subroutine description:
!       -----------------------
!       Calculate mass fluxes (mass per unit time) - vapor and gas - 
!       and their derivatives for given T, RHO and fixed mass fractions.
!       Take into account Knudsen and Reinolds flow, depending on the
!       interaction between the material and the pores.

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: Main.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION FIM(IMX),SIGMA(IMX),COPG(LG),XGJ(LG),XGJP(LG)
      SAVE
 

!  KP=0 - Knudsen flow.
!  KP=1 - Interpolation between Knudsen and Poiseuille (perfect mixture).
!  KP=2 - Transition to Poiseuille at low Knudsen numbers (immiscible flows).


      CONSTK=3.2D1/THRE*SQRT(RGAS/(TWO*PI))
      CONSTP=THRE/3.2D1*SQRT(RGAS/(TWO*PI))
      CONST=(16.D0/THRE)**2/PI
      CONSTD=THRE/(3.2D1*PI)/(6.67D-8*FMTOT0*RHODUS)
      FIM(1)=ZERO
      XMFL(1)=ZERO
      VAFL(1)=ZERO
      VFL(1)=ZERO
      VFR(1)=ZERO
      DVFDTL(1)=ZERO
      DVFDTR(1)=ZERO
      COPV=DH2O**2/SQRT(AI)
      DO M=1,LG
       GAFL(1,M)=ZERO
       GFL(1,M)=ZERO
       GFR(1,M)=ZERO
       DGFDTL(1,M)=ZERO
       DGFDTR(1,M)=ZERO
       COPG(M)=DIAMG(M)**2/SQRT(AG(M))
      END DO
!  Normal setting: maximal dust ejection efficiency at the surface, regardless
!                  of interior dust flow (BDUST>0); otherwise, same efficiency
!                  everywhere (BDUST<0 in data), including BDUST=0.
      BDUSTS=ONE
      IF(BDUST.LT.ZERO)THEN
       BDUST=-BDUST
       BDUSTS=BDUST
      END IF
     

      IF(NTIME.GT.0)THEN

       DO I=IMAXP1,2,-1

        TETI=TET(I)
        PORI=PO(I)
        IF(PORI.LT.0.6D0)THEN
!  Capillary model for porous medium 
         FIMI=PORI/8.D0
        ELSE
!  Randomly packed spherical grains model (optional)
         FIMI=PORI**(THRE/TWO)/(ONE-PORI)**(ONE/3.D0)
        END IF
        IF(I.EQ.IMAXP1)THEN
         IF(TET(I).GT.ZERO)THEN
          TETIP=TET(I)
         ELSE
          TETIP=TETAB
         END IF
         PXMFL=FIMI*TETIP*CONSTK/RGAS*SQV
         FIMI=FIMI*TETI*CONSTK
         FIM(I)=FIMI
        ELSE
         FIMI=FIMI*TETI*CONSTK
         FIM(I)=(FIMP+FIMI)/TWO
        END IF
        FIMP=FIMI
!  Calculate Knudsen and Reynolds numbers
        ZKN(I)=ONE
        RE(I)=ZERO
        XGAS=ROV0(I)/AMICE*DH2O**2
        MCR=0
        FMAX=ABS(VAFL(I))
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          XGAS=XGAS+ROG0(I,M)/AMG(M)*DIAMG(M)**2
          IF(GAFL(I,M).GT.FMAX)THEN
           MCR=M
           FMAX=ABS(GAFL(I,M))
          END IF
         END DO
        END IF
        XGASAV=XGAS
        IF(I.LT.IMAXP1) XGASAV=0.5D0*(XGAS+XGASP)
        IF(XGASAV.GT.ZERO)THEN
         ZKN(I)=ONE/(TWO*PI*XGASAV*TET(I))
         IF(KP.EQ.1.AND.ZKN(I).LT.1.D-2) ZKN(I)=1.D-2
!         IF(KP.EQ.2.AND.ZKN(I).LT.1.D-2) ZKN(I)=1.D-2
         IF(KP.EQ.2.AND.ZKN(I).LT.1.D-5) ZKN(I)=1.D-5
        END IF
        IF(MCR.EQ.0)THEN
         REY=ROV(I)*SQRT(8.D0*T(I)/(PI*AIB))
         IF(REY.GT.ZERO) RE(I)=ABS(VAFL(I)/(FACR(I)*ZKN(I)))/REY
        ELSE
         REY=ROG(I,M)*SQRT(8.D0*T(I)/(PI*AGB(M)))
         IF(REY.GT.ZERO) RE(I)=ABS(GAFL(I,M)/(FACR(I)*ZKN(I)))/REY
        END IF
        XGASP=XGAS
        XGAS=ROV0(I)
        SIGMAI=ROV0(I)*COPV
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          XGAS=XGAS+ROG0(I,M)
          SIGMAI=SIGMAI+ROG0(I,M)*COPG(M)
         END DO
        END IF
        IF(XGAS.GT.ZERO)THEN
         SIGMA(I)=SIGMAI/XGAS
        ELSE
         SIGMA(I)=ZERO
        END IF
        IF(KP.EQ.2) FIM(I)=FIM(I)*(ONE+ONE/(ZKN(I)*CONST))

       END DO

      END IF

!  Calculate Knudsen flux and derivatives

      DVFDTR(IMAXP1)=ZERO
      VFR(IMAXP1)=ZERO
      ACOMAX=ZERO
      DVIP=ZERO
      XVJP=ZERO
      DO M=1,LG
       XGJP(M)=ZERO
       DGFDTR(IMAXP1,M)=ZERO
       GFR(IMAXP1,M)=ZERO
      END DO
      TERP=ZERO
      DTERPT=ZERO

      DO I=IMAXP1,2,-1

       DVI=DV(I)
       DVPDV=0.5D0*(DVI+DVIP)
       XFL=FSQ(I)*FIM(I)/DVPDV
       TI=T(I)
       TER=SQRT(TI)/PO(I)
       DTERT=0.5D0*TER/TI
       XVJ=ROV(I)/SQV
       VFR(I)=XFL*TERP/SQV
       VFL(I)=XFL*TER/SQV
       VAFL(I)=XFL*(XVJP*TERP-XVJ*TER)
       DVFDTL(I)=-XFL*XVJ*DTERT
       DVFDTR(I)=XFL*XVJP*DTERPT
       XVJP=XVJ
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         XGJ(M)=ROG(I,M)/SQG(M)
         GFR(I,M)=XFL*TERP/SQG(M)
         GFL(I,M)=XFL*TER/SQG(M)
         GAFL(I,M)=XFL*(XGJP(M)*TERP-XGJ(M)*TER)
         DGFDTL(I,M)=-XFL*XGJ(M)*DTERT
         DGFDTR(I,M)=XFL*XGJP(M)*DTERPT
         XGJP(M)=XGJ(M)
        END DO
       END IF
!  Calculate dust flux
       XMFL(I)=ZERO
       IF(ROD0(I).GT.ZERO)THEN
        RAMP=ABS(VAFL(I)*SQRT(8.D0*RGAS*T(I)/(PI*AI)))
        SUMF=VAFL(I)
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          RAMP=MAX(RAMP,ABS(GAFL(I,M)*SQRT(8.D0*RGAS*T(I)/(PI*AG(M)))))
          SUMF=SUMF+GAFL(I,M)
         END DO
        END IF
        ACRIT=MIN(TET(I),CONSTD*RAMP)
        ACOVAM=ACRIT/ADMAX
        IF(ACOVAM.GT.ONE) ACOVAM=ONE
!       IF(ACOVAM.LT.1.D-3) ACOVAM=ZERO 
        ETA(I)=ROD0(I)/RHO0(I)*SQRT(ACOVAM)*BDUST
        XMFL(I)=SUMF*ETA(I)
        ACOMAX=MAX(ACOMAX,ACOVAM)
       END IF
       DVIP=DVI
       TERP=TER
       DTERPT=DTERT

      END DO

      IF(KP.EQ.1)THEN

!   Calculate Poiseuille flux and derivatives
!                    and
!   interpolate for intermediate Knudsen numbers
       
       ROTOTP=ZERO
       DVIP=ZERO
       TIP=ZERO
       TERMP=0.5D0*PO(IMAXP1)*SIGMA(IMAXP1)/SQRT(T(IMAXP1))
       DTERMP=-0.5D0*TERMP/T(IMAXP1)
       DO I=IMAXP1,2,-1

        DVI=DV(I)
        DVPDV=DVIP+DVI
        TI=T(I)
        TERM=0.5D0*PO(I)*SIGMA(I)/SQRT(TI)
        DTERM=-0.5D0*TERM/TI
        XFL=0.5D0*FSQ(I)*CONSTP*(PI*TET(I))**2
        ROTOT=ROV0(I)/AMICE
        IF(MMAX.GT.0)THEN
         DO M=1,MMAX
          ROTOT=ROTOT+ROG0(I,M)/AMG(M)
         END DO
        END IF
        DPDV=(ROTOTP*TIP-ROTOT*TI)/DVPDV
        DPDVT=-ROTOT/DVPDV
        DPDVTP=ROTOTP/DVPDV
        AKN=ZKN(I)
        IF(AKN.LT.ONE)THEN
         A1=(ONE-AKN)/0.99D0
         A2=(AKN-1.D-2)/0.99D0
         ROAV=TWO*ROV(I)
         IF(I.LT.IMAXP1) ROAV=ROV(I+1)+ROV(I)
         VFRP=XFL*(TERM+TERMP)*DPDV
         VFLP=-VFRP
         IF(I.EQ.IMAXP1) VFLP=-TWO*VFRP
         VAFLP=ROAV*VFRP
         DVFDTRP=XFL*ROAV*((TERM+TERMP)*DPDVTP+DTERMP*DPDV)
         DVFDTLP=XFL*ROAV*((TERM+TERMP)*DPDVT+DTERM*DPDV)
         VAFL(I)=A1*VAFLP+A2*VAFL(I)
         VFL(I) =A1*VFLP+A2*VFL(I)
         VFR(I)=A1*VFRP+A2*VFR(I)
         DVFDTL(I)=A1*DVFDTLP+A2*DVFDTL(I)
         DVFDTR(I)=A1*DVFDTRP+A2*DVFDTR(I)
         IF(MMAX.GT.0)THEN
          DO M=1,MMAX
           ROAV=TWO*ROG(I,M)
           IF(I.LT.IMAXP1) ROAV=ROG(I+1,M)+ROG(I,M)
           GFRP=XFL*(TERM+TERMP)*DPDV
           GFLP=-GFRP
           IF(I.EQ.IMAXP1) GFLP=-TWO*GFRP
           GAFLP=ROAV*GFRP
           DGFDTRP=XFL*ROAV*((TERM+TERMP)*DPDVTP+DTERMP*DPDV)
           DGFDTLP=XFL*ROAV*((TERM+TERMP)*DPDVT+DTERM*DPDV)
           GAFL(I,M)=A1*GAFLP+A2*GAFL(I,M)
           GFL(I,M)=A1*GFLP+A2*GFL(I,M)
           GFR(I,M)=A1*GFRP+A2*GFR(I,M)
           DGFDTL(I,M)=A1*DGFDTLP+A2*DGFDTL(I,M)
           DGFDTR(I,M)=A1*DGFDTRP+A2*DGFDTR(I,M)
          END DO
         END IF
        END IF
        ROTOTP=ROTOT
        DVIP=DVI
        TIP=TI
        TERMP=TERM
        DTERMP=DTERM

       END DO

      END IF

!  Calculate dust flux from the surface

      IP=IMAXP1
      IF(ROD0(IP).GT.ZERO)THEN
       TIP=T(IP)
       SUMF=VAPFL+VAPFLA
       RAMP=ABS(SUMF*SQRT(8.D0*RGAS*TIP/(PI*AI)))
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         SUMFG=VGFL(M)+VGFLA(M)
         RAMPG=ABS(SUMFG*SQRT(8.D0*RGAS*TIP/(PI*AG(M))))
         RAMP=MAX(RAMP,RAMPG)
         SUMF=SUMF+SUMFG
        END DO
       END IF
       ACRIT=MIN(TET(IP),CONSTD*RAMP)
       ACOVAM=ACRIT/ADMAX
       IF(ACOVAM.GT.ONE) ACOVAM=ONE
!       IF(ACOVAM.LT.1.D-3) ACOVAM=ZERO 
       ETA(IP)=ROD0(IP)/RHO0(IP)*SQRT(ACOVAM)*BDUSTS
       XMFL(IP)=XMFL(IP)-SUMF*ETA(IP)
      END IF
      ACOVAM=ACOMAX

      RETURN
      END
