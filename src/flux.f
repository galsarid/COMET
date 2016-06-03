!********************
      SUBROUTINE FLUX
!********************

!       Subroutine description:
!       -----------------------
!       Calculate heat fluxes and their derivatives as functions of the temperature.

!       Subroutine mapping:
!       -------------------
!       Calls: THCOND, SURFL, IMPACT.
!       Called by: Main.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION CKA(IMX),CKC(IMX),CKD(IMX)
      DIMENSION DCKADT(IMX),DCKCDT(IMX),DCKDDT(IMX)
      SAVE


      IF(IIMPAC.GT.0)THEN
       CALL IMPACT
      ELSE 
       FIMPACT=ZERO
      END IF
      
      IP=IMAXP1
      CALL THCOND(2,CKA,CKC,CKD,DCKADT,DCKCDT,DCKDDT)
      DO I=2,IP
       SUMX=ROA(I)+ROB(I)+ROC(I)+ROD(I)
       XAI=ROA(I)/SUMX
       XBI=ROB(I)/SUMX
       XCI=ROC(I)/SUMX
       XDI=ROD(I)/SUMX
       IF(ICLAT.EQ.0)THEN
        IF((XAI+XBI).EQ.ZERO)THEN
         COND(I)=XCI*CKC(I)+XDI*CKD(I)
         DCONDT(I)=XCI*DCKCDT(I)+XDI*DCKDDT(I)
        ELSE IF(XCI.EQ.ZERO)THEN
         COND(I)=(XAI+XBI)*CKA(I)+XDI*CKD(I)
         DCONDT(I)=(XAI+XBI)*DCKADT(I)+XDI*DCKDDT(I)
        ELSE 
         F1=(XAI+XBI)*CKA(I)
         DF1=(XAI+XBI)*DCKADT(I)
         F2=XCI*CKC(I)
         DF2=XCI*DCKCDT(I)
         COND(I)=ONE/(ONE/F1+ONE/F2)+XDI*CKD(I)
         DCONDT(I)=DF1/(ONE+F1/F2)**2+DF2/(ONE+F2/F1)**2+XDI*DCKDDT(I)
        END IF
       ELSE
        COND(I)=(XCI+XBI)*CKC(I)+XDI*CKD(I)
        DCONDT(I)=(XCI+XBI)*DCKCDT(I)+XDI*DCKDDT(I)
       END IF
       COND(I)=HERTZF*COND(I)
       DCONDT(I)=HERTZF*DCONDT(I)
      END DO
       
      FL(1)=ZERO

!  Boundary condition

      FSUR=ZERO
      TIP=T(IP)
      RIP=R(IP)
      SAREA=FACR(IP)
      PORI=PO(IP)
      PP=ROV(IP)*TIP/AIB/PORI
      
!  Fluxes through a thin boundary layer (commented out to synchronize with Dina) 
!      IF(IB.GT.0)THEN
!       TNS=ZERO
!       PFL=COND(IP)
!       DRMAX=0.5D0*(RIP-R(IP-1))
!       TETIP=TETAB 
!       CALL SURFL(TIP,PP,JMAX,EMIS,SIG,SFLUX,PI,TNS,XIN,DXIN,
!    1             50  ,HSUB,QTOT,PZERO,EXPT,AIB,TETIP,DRMAX,
!    2             FIN,DFIN,NTIME,1,IB,EPSX1,EPSX2,
!    3             PQV,PFL,PXMFL,PAT,VAFLIP,DVFIPL,FSUR)
!       VAPFL=SAREA*QTOT
!       DO M=1,MMAX
!        VGFL(M)=ZERO
!       END DO
!       FL(IP)=SAREA*FIN
!       DFLL(IP)=SAREA*DFIN
!       VAFLIP=SAREA*XIN*ROC0(IP)/RHO0(IP)
!       DVFIPL=SAREA*DXIN*ROC0(IP)/RHO0(IP)
!      END IF

!  Effective area correction, according to J-F Crifo (1997)
!      QVI=ROC0(IP)/RHOICE
!     1    *PZERO*EXP(EXPT/TIP)*SQRT(AIB/(TWO*PI*TIP))
      QVI=(ROC(IP)+ROB(IP))/RHO(IP)
     1    *PZERO*EXP(EXPT/TIP)*SQRT(AIB/(TWO*PI*TIP))
      HSUBQV=HSUB*QVI
      DHSUBQ=HSUB*QVI*(EXPT/TIP+0.5D0)
      ZGAS=ZERO
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        ZGAS=ZGAS+ZGI(IP,M)   
       END DO
      END IF
      IF((ROC0(IP)*ROB0(IP)).GT.ZERO)THEN
       CORRAC=ROC0(IP)/(ROC0(IP)+ROB0(IP))
      ELSE
       CORRAC=ONE
      END IF
      VAPFL=SAREA*QVI*CORRAC
      VAPFLA=SAREA*QVI*(ONE-CORRAC)*(ONE-ZGAS)
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
!        QVIG=ROSG0(IP,M)/RHOSOL(M)
!     1        *PZEROG(M)*EXP(EXPTG(M)/TIP)*SQRT(AGB(M)/(TWO*PI*TIP))
!         QVIG=ROSG(IP,M)/RHO(IP)*PZEROG(M)*EXP(EXPTG(M)/TIP)
!     1         *SQRT(AGB(M)/(TWO*PI*TIP))
        QVIG=PZEROG(M)*EXP(EXPTG(M)/TIP)*SQRT(AGB(M)/(TWO*PI*TIP))
        QMAX=ROSG(IP,M)*DV(IP)/SAREA/DTIME
        IF(QVIG.GT.QMAX)THEN
         QVIG=QMAX
         DQVIGDT=ZERO
        ELSE
         DQVIGDT=QVIG*(EXPTG(M)/TIP+0.5D0)
        END IF
        HSUBQV=HSUBQV+HSUBG(M)*QVIG
        DHSUBQ=DHSUBQ+HSUBG(M)*DQVIGDT
        VGFL(M)=SAREA*QVIG
        VGFLA(M)=SAREA*QVIG*(ONE-CORRAC)*ZGI(IP,M)
       END DO
      END IF

      FL(IP)=SAREA*(SFLUX-EMIS*SIG*TIP**4-HSUBQV+FSUR+FIMPACT)
      DFLL(IP)=-SAREA*(4.D0*EMIS*SIG*TIP**3+DHSUBQ/TIP)

      DVIP=DV(IP)
      TERMP=DVIP*COND(IP)
      DTERMP=DVIP*DCONDT(IP)

!  End boundary condition

      IMAX=IMAXP1-1
      DO I=IMAX,2,-1
       DVI=DV(I)
       TI=T(I)
       DVPDV=DVI+DVIP
       TERM=DVI*COND(I)
       DTERM=DVI*DCONDT(I)
       XKAPA=FSQ(I)*(TERM+TERMP)/DVPDV
       FL(I)=TWO*XKAPA*(TIP-TI)/DVPDV
       DFLR(I)=TWO*(XKAPA+(TIP-TI)*FSQ(I)*DTERMP/DVPDV)/DVPDV
       DFLL(I)=TWO*(-XKAPA+(TIP-TI)*FSQ(I)*DTERM/DVPDV)/DVPDV
       TIP=TI
       DVIP=DVI
       TERMP=TERM
       DTERMP=DTERM
      END DO
      DFLL(1)=ZERO
      DFLR(1)=ZERO
      DFLR(IMAXP1)=ZERO
      IF(IB.EQ.0) FSUR=FL(IMAXP1)/FACR(IMAXP1)

      RETURN
      END
