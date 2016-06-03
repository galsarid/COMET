!******************************************************************
      SUBROUTINE HYDRO(PSOL,RHO,P,DM,DV,IP,G,PI4,EPSX1,EPSHY,NTIME)
!******************************************************************

!       Subroutine description:
!       -----------------------
!       Calculates hydrostatic structure by adjusting volumes (grid points).

!       Subroutine mapping:
!       -------------------
!       Calls: EOS, SGTSL.
!       Called by: ENDSTP.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INTEGER NIT
      DIMENSION PSOL(*),RHO(*),P(*),DM(*),DV(*)
      DIMENSION RO(IMX),DPDRO(IMX),DRODV(IMX),V(IMX),FM(IMX)
      DIMENSION A(IMX),B(IMX),C(IMX),D(IMX),XV(IMX),V0(IMX)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0,FOUR=4.D0)
      SAVE
      DATA NIT/100/, EPSTIM/1.D2/


      F=G*(PI4/8.1D1)**(ONE/THRE)/TWO
      F43=FOUR/THRE
      DO I=1,IP
       XV(I)=ZERO
       A(I)=ZERO
       B(I)=ZERO
       C(I)=ZERO
       D(I)=ZERO
      END DO
      PSOL(1)=ZERO
      RO(1)=ZERO
      V0(1)=ZERO
      V(1)=ZERO
      FM(1)=ZERO
      DO I=2,IP
       V0(I)=V0(I-1)+DV(I)
       V(I)=V0(I)
       FM(I)=FM(I-1)+DM(I)
      END DO
      V(IP+1)=V(IP)
      FM(IP+1)=FM(IP)
      PSOL(IP+1)=ZERO
      DPDRO(IP+1)=ZERO
      DRODV(IP+1)=ZERO

      DO IT=1,NIT

       KCONV=0
       DO I=2,IP
        RO(I)=DM(I)/(V(I)-V(I-1))
        DRODV(I)=(FM(I)-FM(I-1))/(V(I)-V(I-1))**2
       END DO
       CALL EOS(RO,P,PSOL,DPDRO)
       CONV=ZERO
       DO I=2,IP
        IM=I-1
        A(IM)=DPDRO(I+1)*DRODV(I+1)*V(I+1)
        C(IM)=DPDRO(I)*DRODV(I)*V(I-1)
        B(IM)=-(DPDRO(I+1)*DRODV(I+1)+DPDRO(I)*DRODV(I))*V(I)
     1        -(PSOL(I+1)-PSOL(I))*F43
        D(IM)=PSOL(I+1)-PSOL(I)+F*FM(I)*(FM(I+1)-FM(I-1))/V(I)**F43
        CONV=CONV+D(IM)**2
       END DO
       CONV=SQRT(CONV)/FLOAT(IP)
       CALL SGTSL(IP-1,C,B,A,D,INFO)
       IF(INFO.NE.0)THEN
        WRITE(2,*)'ERROR: Matrix inversion (HYDRO). ',
     1            'INFO element of the diagonal=0 (INFO=',INFO,'). ',
     2            'Exiting.'
        WRITE(17,*)NTIME
        WRITE(17,*)'ERROR: Matrix inversion (HYDRO). ', 
     1             'INFO element of the diagonal=0 (INFO=',INFO,'). ',
     2             'Exiting.'
        STOP
       END IF
       DELMAX=ZERO
       IDELM=1
       XVMAX=ZERO
       DO I=2,IP
        XV(I)=D(I-1)
        ABXV=ABS(XV(I))
        XVMAX=MAX(XVMAX,ABXV)
        IF(ABXV.GT.EPSX1) KCONV=KCONV+1
       END DO
       IF(XVMAX.GT.EPSTIM)THEN
        FACT=EPSTIM/XVMAX
       ELSE
        FACT=ONE
       END IF
       DO I=2,IP
        XV(I)=XV(I)*FACT
        V(I)=V(I)*EXP(XV(I))
        DEL=MAX(DELMAX,ABS(V(I)-V0(I))/V0(I))
        IF(DEL.GT.DELMAX)THEN
         IDELM=I
         DELMAX=DEL
        END IF
       END DO
       IF(KCONV.EQ.0.AND.CONV/PSOL(2).LT.EPSHY)THEN
        DO I=2,IP
         DV(I)=V(I)-V(I-1)
         RHO(I)=DM(I)/DV(I)
        END DO
        WRITE(2,*)'HYDRO Successful. ',
     1            'IT,KCONV,IDELM,DELMAX,CONV/PSOL(2) = ',
     2             IT,KCONV,IDELM,DELMAX,CONV/PSOL(2)
        WRITE(17,*)NTIME
        WRITE(17,*)'HYDRO Successful. ',
     1             'IT,KCONV,IDELM,DELMAX,CONV/PSOL(2) = ',
     2             IT,KCONV,IDELM,DELMAX,CONV/PSOL(2)
        RETURN
       END IF

      END DO

      WRITE(2,*)'Failure in HYDRO. ',
     1          'NIT,KCONV,IDELM,DELMAX,CONV/PSOL(2) = ',
     2           NIT,KCONV,IDELM,DELMAX,CONV/PSOL(2)
      WRITE(17,*)NTIME
      WRITE(17,*)'Failure in HYDRO. ',
     1           'NIT,KCONV,IDELM,DELMAX,CONV/PSOL(2) = ',
     2            NIT,KCONV,IDELM,DELMAX,CONV/PSOL(2)

      RETURN
      END

