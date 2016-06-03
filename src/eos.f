!************************************
      SUBROUTINE EOS(RO,P,PSOL,DPDRO)
!************************************

!       Subroutine description:
!       -----------------------
!       Provides equation of state for the solid matrix

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: HYDRO.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0,FOUR=4.D0)
      DIMENSION RO(*),P(*),PSOL(*),DPDRO(*)
      SAVE

      DATA BM1/1.667D0/, BM2/2.333D0/
      DATA V1/0.333D0/, V2/0.667D0/
      DATA S1/-1.333D0/, S2/3.385D0/, S3/0.667D0/

!  Calculation of Birch-Murnaghan approximation (to leading order):
      IF(EOSF.EQ.O)THEN
       WRITE(2,*)'B-M EOS in use.'
       WRITE(17,*)NTIME
       WRITE(17,*)'B-M EOS in use.'
       DO I=2,IMAXP1
        A=1.5D0*PK
        PSOL(I)=A*((RO(I)/RO0)**BM2-(RO(I)/RO0)**BM1)+P(I)
        DPDRO(I)=A/RO(I)*(BM2*(RO(I)/RO0)**BM2-BM1*(RO(I)/RO0)**BM1)
       END DO
       RETURN
      END IF

!  Calculation of Vinet approximation (to leading order):
      IF(EOSF.EQ.1)THEN
       WRITE(2,*)'Vinet EOS in use.'
       WRITE(17,*)NTIME
       WRITE(17,*)'Vinet EOS in use.'
       DO I=2,IMAXP1
        A=3.0D0*PK
        PSOL(I)=A*((RO(I)/RO0)**V2-(RO(I)/RO0)**V1)+P(I)
        DPDRO(I)=A/RO(I)*(V2*(RO(I)/RO0)**V2-V1*(RO(I)/RO0)**V1)
       END DO
       RETURN
      END IF

!  Calculation of Shanker approximation (to leading order):
      IF(EOSF.EQ.2)THEN
       WRITE(2,*)'Shanker EOS in use.'
       WRITE(17,*)NTIME
       WRITE(17,*)'Shanker EOS in use.'
       DO I=2,IMAXP1
        A=-0.69231D0*PK
        PSOL(I)=A*(RO0/RO(I))**S1*((1.D0-RO0/RO(I))*(S2-RO0/RO(I)))+P(I)
        DPDRO(I)=PSOL(I)/(RO0/RO(I))-A/RO0*(RO0/RO(I))**S3
     1           *(S3+1.0D0-2.0D0*(ROO/RO(I)))
       END DO
       RETURN
      END IF
             
      END

