!*************************************
      SUBROUTINE SGTSL(N,C,D,E,B,INFO)
!*************************************

!       Subroutine description:
!       -----------------------
!       Matrix inverse.

!       Subroutine mapping:
!       -------------------
!       Calls: None. 
!       Called by: Main, COMPOS, HYDRO.

      INTEGER N,INFO
      REAL*8 C(*),D(*),E(*),B(*)

!     SGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND
!     SIDE WILL FIND THE SOLUTION.
!
!     ON ENTRY
!
!        N       INTEGER
!                IS THE ORDER OF THE TRIDIAGONAL MATRIX.
!
!        C       REAL(N)
!                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
!                ON OUTPUT C IS DESTROYED.
!
!        D       REAL(N)
!                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
!                ON OUTPUT D IS DESTROYED.
!
!        E       REAL(N)
!                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
!                ON OUTPUT E IS DESTROYED.
!
!        B       REAL(N)
!                IS THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       IS THE SOLUTION VECTOR.
!
!        INFO    INTEGER
!                = 0 NORMAL VALUE.
!                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
!                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
!                    THIS IS DETECTED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!
!     NO EXTERNALS
!     FORTRAN ABS
!
!     INTERNAL VARIABLES
!
      INTEGER K,KB,KP1,NM1,NM2
      REAL*8 T
!     BEGIN BLOCK PERMITTING ...EXITS TO 100
!
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = 0.0D0
            E(N) = 0.0D0
!
            DO 30 K = 1, NM1
               KP1 = K + 1
!
!              FIND THE LARGEST OF THE TWO ROWS
!
               IF (ABS(C(KP1)) .LT. ABS(C(K))) GO TO 10
!
!                 INTERCHANGE ROW
!
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
!
!              ZERO ELEMENTS
!
               IF (C(K) .NE. 0.0D0) GO TO 20
                  INFO = K
!     ............EXIT
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = 0.0D0
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (C(N) .NE. 0.0D0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
!
!           BACK SOLVE
!
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
!
      RETURN
      END
