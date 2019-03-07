!     SUBROUTINE        O  U  T  N  O  D  E    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO PRINT NODE COORDINATES, PRESSURES OR HEADS, CONCENTRATIONS OR 
! ***  TEMPERATURES, AND SATURATIONS IN A FLEXIBLE, COLUMNWISE FORMAT.  
! ***  OUTPUT IS TO UNIT fNOD.                                            
!                                                                       
      SUBROUTINE OUTNODE2(PVEC, UVEC, SW, IN, X, Y, Z, TITLE1, TITLE2) 
      USE ITERAT 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE GRAVEC
      USE PLT1
      USE SutraMSPrecision
      IMPLICIT NONE
      real (DP) :: &
        PVEC (NN), UVEC (NN, NSPE), SW (NN) 
      real (DP) :: &
        X (NN), Y (NN), Z (NN) 
      real (DP) :: &
        VCOL (NCOLMX), VVAR (7 + NSPE-1) 
      CHARACTER(1) TITLE1 (80), TITLE2 (80) 
      CHARACTER(8) HORP 
      CHARACTER(13) TORC (NSPE) 
      CHARACTER(1) CPHORP, CPTORC, CPSATU 
      integer (I4B) :: &
        IN (NIN), IIN (8)
      real (I4B) :: &
        TT (999)
      integer (I4B) :: &
        ITT (999), ISTORC (999), ISHORP (999), &
        ISSATU (999) 
      LOGICAL PRINTN 
      !LOCAL VARIABLES
      CHARACTER (LEN=15) &
        TCHAR
      INTEGER (I4B) :: &
        LCHAR
      INTEGER (I4B) :: &
        TS, &
        LCHORP, LCTORC, &
        I, &
        JT, JTMAX, &
        K, KK, KT, KTMAX, &
        M
      REAL (DP) :: &
        DKTM2, DELTK
      
      INTEGER (I4B) :: A
  
      DO 980 I = 1, NN    
            KK = 0            
            DO K = 6, (5 + NSPE) 
            KK = KK + 1 
            VVAR (K) = UVEC (I, KK)  
            ENDDO
            A = KK - 1                 ! Chengji 2015-04-14
 !           IF(MOD(I,NN1).EQ.1) THEN   ! Chengji 2015-04-14
 !               WRITE(2009888,'(I8,2E15.7)') I, UVEC(I,A), UVEC(I,KK)
 !           ENDIF
			
  980    END DO                                                        
      RETURN 
!                                                                       
      END SUBROUTINE OUTNODE2
