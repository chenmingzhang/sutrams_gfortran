!     SUBROUTINE        OUTVELOCITY            SUTRA-MS VERSION 2004.1

! *** PURPOSE :                                                         
! ***  TO PRINT ELEMENT CENTROID COORDINATES AND VELOCITY COMPONENTS    
! ***  IN A FLEXIBLE, COLUMNWISE FORMAT.  OUTPUT IS TO UNIT fELE.         
!                                                                       
      SUBROUTINE AVGVELOCITY(VMAG, VANG1, VANG2, IN, X, Y, Z, TITLE1, TITLE2) 

      USE ITERAT 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE GRAVEC
      USE PLT1
      USE SutraMSPrecision
      
      USE M_CONTROL

      IMPLICIT NONE

      real (DP) :: &
        VMAG (NE), VANG1 (NE), VANG2 (NEX) 
      real (DP) :: &
        X (NN), Y (NN), Z (NN) 
      real (DP) :: &
        VCOL (NCOLMX), VVAR (7) 
      CHARACTER(1) TITLE1 (80), TITLE2 (80) 
      CHARACTER(15) COLTK6 (7) 
      CHARACTER(1) CPVX, CPVY, CPVZ 
      integer (I4B) :: &
        IN (NIN), IIN (8)
      real (DP) :: &
        TT (999)
      integer (I4B) :: &
        ITT (999) , &
        ISVEL (999) 
      LOGICAL PRINTE 

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        II, III, &
        JT, JTMAX, &
        KT, KTMAX, &
        L, LL, LCP, LCV, &
        M, MM
      REAL (DP) :: &
        CENTRX, CENTRY, CENTRZ, CVA2, &
        DKTM2, DELTK, &
        RN48, &
        TS, &
        VA1, VA2, VECTRX, VECTRY, VECTRZ
      
      
      !LOCAL VARIABLES Chengji 2015-08-10
      REAL (DP):: XET(NE),YET(NE),ZET(NE)
      REAL (DP):: VXT(NE),VYT(NE),VZT(NE)
      INTEGER (I4B):: I                                                                
!.....Calculate headers on time step 1                                  
!.....and create output on time steps greater than or equal to 1        
      IF (IT.EQ.0) RETURN 
!                                                                                                                  
        DO 2013 I=1,NE
         IF(MOD(IT,NSTEP).EQ.1) THEN
           XET(I)=0
           YET(I)=0
           ZET(I)=0
           VXT(I)=0
           VYT(I)=0
           VZT(I)=0
         ENDIF
  2013  CONTINUE 
!                                                                       
!.....  The velocity data for this time step                            
         RN48 = 1D0 / DBLE (N48) 
         DO 50 L = 1, NE 
            CENTRX = 0D0 
            CENTRY = 0D0 
            CENTRZ = 0D0 
            DO 40 II = 1, N48 
               III = II + (L - 1) * N48 
               IIN (II) = IN (III) 
               CENTRX = CENTRX + X (IIN (II) ) 
               CENTRY = CENTRY + Y (IIN (II) ) 
               CENTRZ = CENTRZ + Z (IIN (II) ) 
   40       END DO 
            CENTRX = CENTRX * RN48 
            CENTRY = CENTRY * RN48 
            CENTRZ = CENTRZ * RN48 
            VA1 = 0.017453292D0 * VANG1 (L) 
            LL = MIN (L, NEX) 
            VA2 = 0.017453292D0 * VANG2 (LL) * DKTM2 
            CVA2 = DCOS (VA2) 
            VECTRX = VMAG (L) * DCOS (VA1) * CVA2 
            VECTRY = VMAG (L) * DSIN (VA1) * CVA2 
            VECTRZ = VMAG (L) * DSIN (VA2) 
            VVAR (1) = DBLE (L) 
            VVAR (2) = CENTRX 
            VVAR (3) = CENTRY 
            VVAR (4) = CENTRZ 
            VVAR (5) = VECTRX 
            VVAR (6) = VECTRY 
            VVAR (7) = VECTRZ 
            DO 984 M = 1, NCOLS6 
               VCOL (M) = VVAR (J6COL (M) ) 
  984       END DO
  
 ! ******** SUM UP THE VELOCITY OVER ONE TIDAL CYCLE ********
            XET(L) = XET(L) + CENTRX
            YET(L) = YET(L) + CENTRY
            ZET(L) = ZET(L) + CENTRZ
            VXT(L) = VXT(L) + VECTRX
            VYT(L) = VYT(L) + VECTRY
            VZT(L) = VZT(L) + VECTRZ
  ! *********************************************************
   50    END DO 
      
      DO 11 I=1,NE
       IF(MOD(IT,NSTEP).EQ.0) THEN
         WRITE(7,'(6E15.7)') XET(I)/NSTEP,YET(I)/NSTEP,ZET(I)/NSTEP, &
                             VXT(I)/NSTEP,VYT(I)/NSTEP,VZT(I)/NSTEP
       ENDIF
   11 CONTINUE
!                                                                       
 9999 CONTINUE 
      RETURN 
!                                                                       
      END SUBROUTINE AVGVELOCITY
