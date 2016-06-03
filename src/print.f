!*********************
      SUBROUTINE PRINT
!*********************

!       Subroutine description:
!       -----------------------
!       Print out full data to COMOUT file.

!       Subroutine mapping:
!       -------------------
!       Calls: STRESS, TSTRESS, TABLE, SURFL.
!       Called by: Main, ENDSTP, SURFL, INPUT.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION PRIMAT(IMX,11)
!      DIMENSION SRR(IMX),STT(IMX)
      SAVE
      DATA JPRI/0/

      JPRI=JPRI+1
      FM=ZERO
      IL=2
      IJ=1
      UTOT=ZERO
!      CALL STRESS(IMAXP1,PI43,P,R,SRR,STT,NTIME)
!      CALL TSTRES(IMAXP1,PI43,T,R,SRR,STT)
      RIP=R(IMAXP1)

      DO I=2,IMAXP1
       UTOT=UTOT+U(I)*DV(I)
       FM=FM+DM(I)
       IF(MOD(JPRI,10).NE.0.AND.IMAXP1.GT.300)THEN
        IF(ABS(T(I)-T(2))/T(2).LT.EPSTIM) IL=MAX(I,IMAXP1-ITABLE)
       END IF
       PRIMAT(I,1)=DM(I)
       PRIMAT(I,2)=(RIP-R(I-1))/1.D2
       PRIMAT(I,3)=P(I)
       PRIMAT(I,4)=T(I)
       PRIMAT(I,5)=RHO(I)
       PRIMAT(I,6)=TET(I)
       PRIMAT(I,7)=PO(I)
       PRIMAT(I,8)=GAFL(I,1)/AMG(1)
       PRIMAT(I,9)=GAFL(I,2)/AMG(2)
       PRIMAT(I,10)=VAFL(I)/AMICE
       PRIMAT(I,11)=FL(I)
      END DO

      IF(IMAXP1-IL.GT.200) IJ=2

      WRITE(2,50)TIME,NORB,DIST,UTOT,FM,RIP/1.D5
      WRITE(17,*)NTIME
      WRITE(17,50)TIME,NORB,DIST,UTOT,FM,RIP/1.D5
  50  FORMAT(//7X,8H TIME  =,1PE13.6,16H ORBIT  NUMBER =,I6,
     1         10H DISTANCE=,1PE12.5,6H UTOT=,1PE12.5,6H MTOT=,1PE12.5/,
     2         7H R(KM)=,1PE12.5/)

      WRITE(2,'(6X,11H    DMI    ,11H  DEPTH(M) ,11H  PRESSURE ,
     1             11H  TEMPER   ,11H  DENSITY  ,11H   TETA    ,
     2             11H  POROSITY ,11H  CO FLUX  ,11H  CO2 FLUX ,
     3             11H  VAP FLUX ,11H    FLUX   /)')
      KD=11
      CALL TABLE(KD,IL,IMAXP1,IJ,PRIMAT)
 
      DO I=IL,IMAXP1
       PRIMAT(I,1)=ROA(I)
       PRIMAT(I,2)=ROC(I)
       PRIMAT(I,3)=ROB(I)
       PRIMAT(I,4)=ROV(I)
       PRIMAT(I,5)=XMFL(I)
       PRIMAT(I,6)=COND(I)
       PRIMAT(I,7)=ELHS(I)
       PRIMAT(I,8)=ERHS(I)
       PRIMAT(I,9)=S(I)
       PRIMAT(I,10)=SVG(I)
       PRIMAT(I,11)=DUDT(I)
      END DO

      WRITE(2,'(//6X,11H   ICE-A   ,11H   ICE-C   ,11H   ICE-B   ,
     1               11H   RHO-V   ,11H   DST-F   ,11H   THCON   ,
     2               11H   F-CON   ,11H   F-ADV   ,11H   H-AMC   ,
     3               11H   H-VAP   ,11H   RDUDT   /)')
      KD=11
      CALL TABLE(KD,IL,IMAXP1,IJ,PRIMAT)

      IF(MMAX.GT.0)THEN

       DO I=2,IMAXP1
        DO M=1,MMAX
         PRIMAT(I,M)=ROG(I,M)
        END DO
       END DO

       WRITE(2,'(//6X,11H   CO--V   ,11H   CO2-V   ,11H   HCN-V   ,
     1                11H   NH3-V   ,11H   C2H2V   ,11H   CH4-V   ,
     2                11H   H2S-V   ,11H   CH3OH   ,11H   C2H6V   /)')
       IF(MMAX.LT.11)THEN
        KD=MMAX
       ELSE
        KD=11
       END IF
       CALL TABLE(KD,IL,IMAXP1,IJ,PRIMAT)

       DO I=2,IMAXP1
        DO M=1,MMAX
         PRIMAT(I,M)=ROSG(I,M)
        END DO
       END DO

       WRITE(2,'(//6X,11H   COICE   ,11H   CO2IC   ,11H   HCNIC   ,
     1                11H   NH3IC   ,11H   C2H2I   ,11H   CH4IC   ,
     2                11H   H2SIC   ,11H   CH3OH   ,11H   C2H6I   /)')
       IF(MMAX.LT.11)THEN
        KD=MMAX
       ELSE
        KD=11
       END IF
       CALL TABLE(KD,IL,IMAXP1,IJ,PRIMAT)

      END IF

      DO I=2,IMAXP1
        PRIMAT(I,1)=PSOL(I)
        PRIMAT(I,2)=ALFP(I)
        PRIMAT(I,3)=XM0(I)
        PRIMAT(I,4)=ZKN(I)
        PRIMAT(I,5)=RE(I)
        PRIMAT(I,6)=PCR(I)
      END DO

      WRITE(2,'(//6X,11H   P-SOL   ,11H   ALF-D   ,11H   WATER   ,
     1               11H   KNUDS   ,11H   REYNO   ,11H   PCRIT   ,
     2               11H           ,11H           ,11H           /)')
       KD=6
       CALL TABLE(KD,IL,IMAXP1,IJ,PRIMAT)


      IF(LMAX.GT.0)THEN
       DO I=2,IMAXP1
        DMI=DM(I)
        DO L=1,LMAX
         PRIMAT(I,L)=XMRAD(I,L)/DMI
        END DO
       END DO
       WRITE(2,'(//6X,11H   K40**   ,11H   TH232   ,11H   U238*   ,
     1                11H   U235*   ,11H   AL26*   ,11H   FE60*   ,
     2                11H           ,11H           ,11H           /)')
       IF(LMAX.LT.11)THEN
        KD=LMAX
       ELSE
        KD=11
       END IF
       CALL TABLE(KD,IL,IMAXP1,IJ,PRIMAT)
      END IF

!      IF(JMAX.GT.0.AND.NTIME.GT.0)THEN
!       IF(TET(IMAXP1).GT.ZERO)THEN
!        TETIP=TET(IMAXP1)
!       ELSE
!        TETIP=TETAB
!       END IF
!       CALL SURFL(0.D0,0.D0,JMAX,EMIS,SIG,SFLUX,PI,0.D0,0.D0,0.D0,
!    1             50  ,HSUB,QTOT,PZERO,EXPT,AIB,TETIP,0.D0,
!    2             0.D0,0.D0,NTIME,0     ,IB,EPSX1,EPSX2,
!    3             PQV,PFL,PXMFL,PAT,VAFLIP,DVFIPL,FSUR)
!      END IF

      RETURN
      END

