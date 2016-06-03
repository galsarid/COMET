!***************************************************************
      SUBROUTINE SETEOS(PSOL,RHO,P,DM,DV,IP,G,PI4,PK,RO0,RADIUS,
     1                  EPSX1,EPSHY)
!***************************************************************


!       Subroutine description:
!       -----------------------
!       Calculates the constants of the EOS (PK and RO0) for a given RADIUS, 
!       constrained by the requirement that porosity be larger than 0.02 everywhere


!       Subroutine mapping:
!       -------------------
!       Calls: HYDRO.
!       Called by: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INTEGER NIT
      DIMENSION PSOL(*),RHO(*),P(*),DM(*),DV(*)
      DIMENSION RO(IMX),DPDRO(IMX),DRODV(IMX),V(IMX),FM(IMX)
      DIMENSION A(IMX),B(IMX),C(IMX),D(IMX),XV(IMX),V0(IMX)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0,FOUR=4.D0)
      SAVE
      DATA NIT/100/


      PK1=PK
      DO N=1,NIT
       PK1=PK1/TWO
       CALL HYDRO(PSOL,RHO,P,DM,DV,IP,G,PI4,EPSX1,EPSHY,PK1,RO0)
       V(1)=ZERO
       DO I=2,IP
        V(I)=V(I-1)+DV(I)
       END DO
       R1=(THRE*V(IP)/PI4)**(ONE/THRE)
       WRITE(2,*)'Lower bound', N, PK1, R1, RADIUS
       IF(RADIUS.GE.R1) GO TO 10
      END DO

  10  CONTINUE
      PK2=PK
      DO N=1,NIT
       PK2=PK2*TWO
       CALL HYDRO(PSOL,RHO,P,DM,DV,IP,G,PI4,EPSX1,EPSHY,PK2,RO0)
       V(1)=ZERO
       DO I=2,IP
        V(I)=V(I-1)+DV(I)
       END DO
       R2=(THRE*V(IP)/PI4)**(ONE/THRE)
       WRITE(2,*)'Upper bound', N, PK2, R2, RADIUS
       IF(RADIUS.LE.R2) GO TO 20
      END DO

  20  CONTINUE
      DO N=1,NIT
       PK=SQRT(PK1*PK2)
       CALL HYDRO(PSOL,RHO,P,DM,DV,IP,G,PI4,EPSX1,EPSHY,PK,RO0)
       V(1)=ZERO
       DO I=2,IP
        V(I)=V(I-1)+DV(I)
       END DO
       R3=(THRE*V(IP)/PI4)**(ONE/THRE)
       WRITE(2,*)'Iteration', N, PK, R3, RADIUS
       IF(ABS(RADIUS-R3).LT.EPSHY*RADIUS) GO TO 30
       IF(RADIUS.LT.R3)THEN
        R2=R3
        PK2=PK
       ELSE
        R1=R3
        PK1=PK
       END IF
      END DO
      
      WRITE(2,*)'WARNING: No solution for PK (SETEOS). 
     1           PK,R3,RADIUS = ',PK,R3,RADIUS
      WRITE(17,*)'WARNING: No solution for PK (SETEOS).  
     1           PK,R3,RADIUS = ',PK,R3,RADIUS
      IP=0

  30  RETURN
      END
