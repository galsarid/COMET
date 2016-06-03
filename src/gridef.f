!******************************
      SUBROUTINE GRIDEF(DELTAV)
!******************************

!       Subroutine description:
!       -----------------------
!       Move grid boundaries inward by DELTAV (provided the 5 innermost
!       shells are large enough). First, boundaries are advanced down to
!       the innermost shell; when it becomes too small - down to the second
!       shell, and so on, until the the fifth. Old shells are stored in
!       DVP(I).

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION DVP(IMX),V(IMX),ROGI(LG),ROSGII(LG),ROSGPI(LG),ZGII(LG)
      DIMENSION XMRADI(LR)
      SAVE

!  Define tests:

      utota=zero
      fmasa=zero
      ficea=zero
      frada=zero
      do i=1,imaxp1
       utota=utota+u(i)*dv(i)
       fmasa=fmasa+rho(i)*dv(i)
       ficea=ficea+roice(i)*dv(i)
       if(lmax.gt.0)then
        do l=1,lmax
         frada=frada+xmrad(i,l)
        end do
       end if
      end do

      IMINP=2
      IMIN=1
      IOUT=10
      DO I=IOUT,2,-1
       IF(DELTAV.LT.DV(I))THEN
        IMIN=I
       END IF
      END DO
      IF(IMIN.NE.IMINP.OR.IMIN.EQ.1)THEN
       WRITE(2,*)'WARNING: Too much mass lost during time-step. 
     1            NTIME,IMIN,IMINP,DLETAV,(DV(I),I=2,IOUT) = ',
     2            NTIME,IMIN,IMINP,DELTAV,(DV(I),I=2,IOUT)
       WRITE(17,*)NTIME
       WRITE(17,*)'WARNING: Too much mass lost during time-step. 
     1            IMIN,IMINP,DLETAV,(DV(I),I=2,IOUT) = ',
     2            IMIN,IMINP,DELTAV,(DV(I),I=2,IOUT)
       IF(IMIN.EQ.1)STOP
      END IF

      U(1)=ZERO
      V(1)=ZERO
      DO I=2,IMAXP1
       V(I)=V(I-1)+DV(I)
      END DO
      DO I=IMIN,IMAXP1
       V(I)=V(I)-DELTAV
      END DO
      DO I=2,IMAXP1
       DVP(I)=DV(I)
       DV(I)=V(I)-V(I-1)
       R(I)=(V(I)/PI43)**(ONE/THRE)
       FACR(I)=PI4*R(I)**2
       FSQ(I)=FACR(I)*FACR(I)
      END DO

!  Interpolate to define new (aditive) conserved quantities
!  (energy and component masses in each phase)

      UI=ZERO
      ROAI=ZERO
      ROBI=ZERO
      ROCI=ZERO
      ROVI=ZERO
      RODI=ZERO
      IF(MMAX.GT.0)THEN
       DO M=1,MMAX
        ROGI(M)=ZERO
        ROSGII(M)=ZERO
        ZGII(M)=ZERO
       END DO
      END IF
      IF(LMAX.GT.0)THEN
       DO L=1,LMAX
        XMRADI(L)=ZERO
       END DO
      END IF
      RHOI=ZERO

      DO I=IMAXP1,IMIN+1,-1

       U(I)=(DVP(I)*U(I)+DELTAV*(U(I-1)-UI))/DV(I)
       UI=U(I-1)
       ROB(I)=(DVP(I)*ROB(I)+DELTAV*(ROB(I-1)-ROBI))/DV(I)
       ROBI=ROB(I-1)
       ROC(I)=(DVP(I)*ROC(I)+DELTAV*(ROC(I-1)-ROCI))/DV(I)
       ROCI=ROC(I-1)
       ROV(I)=(DVP(I)*ROV(I)+DELTAV*(ROV(I-1)-ROVI))/DV(I)
       ROVI=ROV(I-1)
       ROD(I)=(DVP(I)*ROD(I)+DELTAV*(ROD(I-1)-RODI))/DV(I)
       RODI=ROD(I-1)
       RHO(I)=ROB(I)+ROC(I)+ROV(I)+ROD(I)
       ROICE(I)=ROB(I)+ROC(I)
       ROAIP=DVP(I)*ROA(I)+DELTAV*(ROA(I-1)-ROAI)

       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         IF(ROAIP.GT.ZERO)THEN
          ZGI(I,M)=(DVP(I)*ROA(I)*ZGI(I,M)+
     1              DELTAV*(ROA(I-1)*ZGI(I-1,M)-ROAI*ZGII(M)))/ROAIP
         ELSE
          ZGI(I,M)=ZERO
         END IF
         ZGII(M)=ZGI(I-1,M)
         ROSG(I,M)=(DVP(I)*ROSG(I,M)+DELTAV*(ROSG(I-1,M)-ROSGII(M)))
     1             /DV(I)
         ROSGII(M)=ROSG(I-1,M)
         ROSGPI(M)=ROSG(I,M)
         ROG(I,M)=(DVP(I)*ROG(I,M)+DELTAV*(ROG(I-1,M)-ROGI(M)))/DV(I)
         ROGI(M)=ROG(I-1,M)
         RHO(I)=RHO(I)+ROSG(I,M)+ROG(I,M)
         ROICE(I)=ROICE(I)+ROSG(I,M)
        END DO
       END IF

       ROA(I)=ROAIP/DV(I)
       ROAI=ROA(I-1)
       RHO(I)=RHO(I)+ROA(I)
       RHOI=RHO(I-1)
       ROICE(I)=ROICE(I)+ROA(I)
       PO(I)=POR(ROA(I),ROB(I),ROC(I),ROD(I),RHOICE,RHODUS,
     1           MMAX,ROSGPI,RHOSOL)

      END DO

!  Calculate temperatures on the new grid

      DO I=IMAXP1,IMIN,-1
       SUMX=ROV(I)/AIB
       TI=T(I)
       UI=U(I)
       UCON=ROICE(I)*BETA+ROV(I)*CV+ROD(I)*CD
       IF(MMAX.GT.0)THEN
        DO M=1,MMAX
         UCON=UCON+ROG(I,M)*CG(M)
         SUMX=SUMX+ROG(I,M)/AGB(M)
        END DO
       END IF
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
       P(I)=T(I)*SUMX/PO(I)
       IF(LMAX.GT.0)THEN
        DO L=1,LMAX
         IF(I.GT.IMIN)THEN
          XMRAD(I,L)=XMRAD(I,L)+
     1               DELTAV*(XMRAD(I-1,L)/DVP(I-1)-XMRADI(L)/DVP(I))
         ELSE
          XMRAD(I,L)=XMRAD(I,L)-DELTAV*XMRADI(L)/DVP(I)
         END IF
         XMRADI(L)=XMRAD(I-1,L)
        END DO
       END IF
      END DO

!  Perform tests:

      utotb=zero
      fmasb=zero
      ficeb=zero
      fradb=zero
      do i=1,imaxp1
       DM(I)=RHO(I)*DV(I)
       utotb=utotb+u(i)*dv(i)
       fmasb=fmasb+rho(i)*dv(i)
       ficeb=ficeb+roice(i)*dv(i)
       if(lmax.gt.0)then
        do l=1,lmax
         fradb=fradb+xmrad(i,l)
        end do
       end if
      end do
      utotb=utotb-utota
      fmasb=fmasb-fmasa
      ficeb=ficeb-ficea
      fradb=fradb-frada
      if(abs(utotb/utota).gt.epsx1.or.abs(fmasb/fmasa).gt.epsx2)then
        write(2,*)'grid',deltav,utota,utotb,fmasa,fmasb,ficea,ficeb,
     1                   frada,fradb
      end if

      RETURN
      END
