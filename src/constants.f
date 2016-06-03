!**************************
      SUBROUTINE  CONSTANTS
!**************************

!       Subroutine description:
!       -----------------------
!       In this subroutine we define the parameters of all different
!       species that the comet is comprised of.

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: Main.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      SAVE

      PI=4.0D0*ATAN(ONE)
      PI4=4.0D0*PI
      PI43=PI4/THRE
      SIG=5.67051D-5
      BOLK=1.380658D-16
      UMAS=1.6605402D-24
      RGAS=BOLK/UMAS
      G=6.67259D-8
      GM=1.328D26
      ONEAU=1.496D13
      YEAR=3.15576D7
!  Maximum time allowed for the simulation to run (age of the SS):
      TIMAX=1.5E17
!  Minimun time allowed (convergence check on T, in MAIN):
      TIMIN=5.0D0
!  Solar luminosity divided by 4*pi:
      SLM=3.06D32

!  Total kinetic energy of the Deep Impact projectile:
      ENTOTDI=1.92D17

!  Dust:

      CD=1.3D7
      RHODUS=3.25D0
!   Silicate dust + organics conductivity
      CONDD=2.0D4
!   Silicate dust conductivity
!      CONDD=1.0D6

!  H2O:

      AI=1.80154D1
      SQV=SQRT(AI)
      AMICE=AI*UMAS
      AIB=AMICE/BOLK

!  Vapor:

      CV=3.0D0/AIB
      PZERO=3.56D13
      EXPT=-6.141667D3
      HSUB=-EXPT*RGAS/AI
      RHOICE=0.917D0
      ALFA=3.75D4
      BETA=9.0D5
      DH2O=2.9D-8
      TSUBV=9.0D1

!  Amorphous-Crystalline ice phase transition:
!  (used in SOURCE)

      HAMC=9.0D8
      ATAU=1.05D13
      EXPTAU=5.370D3
      TMIN=3.D1

!  Water:

      CW=4.186D7
      TMELT=2.7316D2
      HMELT=3.34D9-CW*TMELT+ALFA*TMELT*TMELT+BETA*TMELT
      TML=0.8D0*TMELT
      EXMELT=4.0D1
      CONDW=5.5D4

!  CO (carbon monoxide):

      AG(1)=2.8D1
      PZEROG(1)=1.2631D10
      EXPTG(1)=-7.6416D2
      RHOSOL(1)=1.1D0
      DIAMG(1)=3.4D-8
      TSUB(1)=1.5D1
!  (From Meech & Svoren, Comets II)
!      TSUB(1)=2.5D1

!  CO2 (carbon dioxide):

      AG(2)=4.4D1
      PZEROG(2)=1.079D13
      EXPTG(2)=-3.1484D3
      RHOSOL(2)=1.56D0
      DIAMG(2)=4.25D-8
      TSUB(2)=4.0D1
!  (From Meech & Svoren, Comets II)
!      TSUB(2)=8.0D1

!  HCN (hydrogen cyanide):

      AG(3)=2.7D1
      PZEROG(3)=3.8665D11
      EXPTG(3)=-4.02466D3
      RHOSOL(3)=ONE
      DIAMG(3)=3.6D-8
      TSUB(3)=2.1D1
!  (From Meech & Svoren, Comets II)
!      TSUB(3)=9.5D1

!  NH3 (ammonia):

      AG(4)=1.7D1
      PZEROG(4)=6.1412D12
      EXPTG(4)=-3.6036D3
      RHOSOL(4)=ONE 
      DIAMG(4)=3.6D-8
      TSUB(4)=2.1D1
!  (From Meech & Svoren, Comets II)
!      TSUB(4)=7.8D1

!  C2H2 (acetylene):

      AG(5)=2.6D1
      PZEROG(5)=9.831D11
      EXPTG(5)=-2.6136D3
      RHOSOL(5)=ONE
      DIAMG(5)=2.4D-8
      TSUB(5)=2.1D1
!  (From Meech & Svoren, Comets II)
!      TSUB(5)=5.7D1

!  CH4 (methane):

      AG(6)=1.6D1
      PZEROG(6)=5.97D10
      EXPTG(6)=-1.1902D3
      RHOSOL(6)=ONE
      DIAMG(6)=4.0D-8
      TSUB(6)=1.5D1
!  (From Meech & Svoren, Comets II)
!      TSUB(6)=3.1D1

!  H2S (hydrogen sulfide):

      AG(7)=3.4D1
      PZEROG(7)=1.2631D11
      EXPTG(7)=-2.64842D3
      RHOSOL(7)=ONE
      DIAMG(7)=3.0D-8
      TSUB(7)=2.0D1
!  (From Meech & Svoren, Comets II)
!      TSUB(7)=5.7D1

!  CH3OH (methanol):

      AG(8)=3.2D1
      PZEROG(8)=8.883D11
      EXPTG(8)=-4.632D3
      RHOSOL(8)=ONE
      DIAMG(8)=4.1D-8
      TSUB(8)=4.0D1
!  (From Meech & Svoren, Comets II)
!      TSUB(8)=9.9D1

!  C2H6 (ethane):

      AG(9)=3.0D1
      PZEROG(9)=4.59D9
      EXPTG(9)=-1.938D3
      RHOSOL(9)=ONE
      DIAMG(9)=4.75D-8
      TSUB(9)=2.0D1
!  (From Meech & Svoren, Comets II)
!      TSUB(9)=4.4D1

! Rest of the compounds (if LG>9):

      DO M=10,LG
       AG(M)=ZERO
       PZEROG(M)=ZERO
       EXPTG(M)=ZERO
       RHOSOL(M)=ZERO
       DIAMG(M)=ZERO
       TSUB(M)=ZERO
      END DO

      DO M=1,LG
       HSUBG(M)=-EXPTG(M)*RGAS/AG(M)
       SQG(M)=SQRT(AG(M))
       AMG(M)=AG(M)*UMAS
       AGB(M)=AMG(M)/BOLK
       IF(M.EQ.1)THEN
        CG(M)=2.5D0/AGB(M)
       ELSE
        CG(M)=3.0D0/AGB(M)
       END IF
      END DO

!  Radionuclides Data:

!  40K:

      TAUR(1)=5.7435D16
      HRAD(1)=1.72D16

!  232Th:

      TAUR(2)=6.3115D17
      HRAD(2)=1.65D17

!  238U:

      TAUR(3)=2.0513D17
      HRAD(3)=1.92D17

!  235U:
      
      TAUR(4)=3.2504D16
      HRAD(4)=1.86D17

!  26Al:

      TAUR(5)=3.3451D13
      HRAD(5)=1.48D17

!  60Fe:

      TAUR(6)=6.9427D13
      HRAD(6)=4.92D16


      RETURN
      END
