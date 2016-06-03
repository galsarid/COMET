!***********************************************************
      SUBROUTINE THCOND(IO,CKA,CKC,CKD,DCKADT,DCKCDT,DCKDDT)
!***********************************************************

!       Subroutine description:
!       -----------------------
!       Calculate thermal conductivity coefficients, corrected for porosity.

!       Subroutine mapping:
!       -------------------
!       Calls: None.
!       Called by: FLUX.

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'dimfile.h'
      INCLUDE 'commonfile.h'
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THRE=3.D0)
      DIMENSION CKA(*),CKC(*),CKD(*)
      DIMENSION DCKADT(*),DCKCDT(*),DCKDDT(*)
      SAVE

      DO I=IO,IMAXP1

       TI=T(I)

!  Porosity correction for thermal conductivity in a fractal structure 
!  (from Shoshany et al.2002)
    
!       PFRAC=0.3D0
!       PERC=0.7D0
!       A=4.1D0
!       B=0.22D0    
!       IF(PO(I).GT.PFRAC)THEN
!        RANK=LOG(ONE-PO(I))/LOG(ONE-PFRAC)
!        POROS=(ONE-PFRAC/PERC)**(RANK*(A*PFRAC+B))
!       ELSE
!        POROS=(ONE-PO(I)/PERC)**(A*PO(I)+B)
!       ENDIF

!  Porosity correction of thermal conductivity for a simplified case 
!  (Prialnik et al. Comets II, Smoluchowski 1981)

       POROS=ONE-PO(I)**(TWO/THRE)

       CKA(I)=POROS*3.13D-3*(TWO*ALFA*TI+BETA)*RHOICE
       DCKADT(I)=POROS*3.13D-3*TWO*ALFA*RHOICE
       CKC(I)=POROS*5.67D7/TI*(ONE-XM0(I))+XM0(I)*CONDW
       DCKCDT(I)=-POROS*5.67D7/(TI*TI)*(ONE-XM0(I))
       CKD(I)=POROS*CONDD+EFRAD* 4.D0*SIG*EMIS*TET(I)*TI**3
       DCKDDT(I)=EFRAD*12.D0*SIG*EMIS*TET(I)*TI**2
      
      END DO

      RETURN
      END

