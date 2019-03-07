!     SUBROUTINE        B  C  T  I  M  E       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  USER-PROGRAMMED SUBROUTINE WHICH ALLOWS THE USER TO SPECIFY:     
! ***   (1) TIME-DEPENDENT SPECIFIED PRESSURES AND TIME-DEPENDENT       
! ***       CONCENTRATIONS OR TEMPERATURES OF INFLOWS AT THESE POINTS   
! ***   (2) TIME-DEPENDENT SPECIFIED CONCENTRATIONS OR TEMPERATURES     
! ***   (3) TIME-DEPENDENT FLUID SOURCES AND CONCENTRATIONS             
! ***       OR TEMPERATURES OF INFLOWS AT THESE POINTS                  
! ***   (4) TIME-DEPENDENT ENERGY OR SOLUTE MASS SOURCES                
!                                                                       
      SUBROUTINE BCTIME (IPBCT, IUBCT, IQSOPT, IQSOUT, GNUP, PM1, GNUPCJ, TDLEVEL, OUTSEEP)  ! CHENGJI 2013-09-02
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE GRAVEC
      USE SutraStorage, ONLY : IPBC, PBC, IUBC, UBC, &
                               QIN, UIN, QUIN, IQSOP, IQSOU, &
                               X, Y, Z, &
                               SpecifiedPBC, &
                               MultiSpeciesBC
      USE M_TIDE
      USE M_SEEPAGE
      USE M_CONTROL
      USE SutraMSPrecision
      USE M_SEASONALT
     
	  
      IMPLICIT NONE

      integer (I4B) :: &
        IPBCT, IUBCT, IQSOPT, IQSOUT
		
!============================================= CHENGJI 2013-09-02	
      REAL(DP):: GNUPCJ, GNUP, PM1
      REAL(DP):: TDLEVEL
      REAL(DP):: SEEPX,SEEPZ
      INTEGER(I4B) :: OUTSEEP ! CHENGJI 2015-08-03
      DIMENSION GNUPCJ(NBCN) 
      DIMENSION PM1(NN)
!================================================================
!================================================= !MT - 04112017	
      REAL(DP):: DRDTSW, SWT
      REAL(DP):: TIMEIND
      REAL(DP):: RHOSWB, RHOFWB
      INTEGER(I4B):: NOMTH, NOYEAR
!================================================================	  
!     LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, IP, IU, IUP, IQP, IQU, &
        K, &
        NSOPI, NSOUI
        
      IF(SEEP.EQ.1) THEN
        OPEN(8,FILE='SEEPAGE_FACE.DAT')
        SEEPX=XSTART
        SEEPZ=ZSTART
      ENDIF

!.....DEFINITION OF REQUIRED VARIABLES                                  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     NN      = EXACT NUMBER OF NODES IN MESH                           
!     NPBC    = EXACT NUMBER OF SPECIFIED PRESSURE NODES                
!     NUBC(K) = EXACT NUMBER OF SPECIFIED CONCENTRATION                 
!               OR TEMPERATURE NODES FOR EACH SPECIES                   
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IT = NUMBER OF CURRENT TIME STEP                                  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     TSEC = TIME AT END OF CURRENT TIME STEP IN SECONDS                
!     TMIN = TIME AT END OF CURRENT TIME STEP IN MINUTES                
!     THOUR = TIME AT END OF CURRENT TIME STEP IN HOURS                 
!     TDAY = TIME AT END OF CURRENT TIME STEP IN DAYS                   
!     TWEEK = TIME AT END OF CURRENT TIME STEP IN WEEKS                 
!     TMONTH = TIME AT END OF CURRENT TIME STEP IN MONTHS               
!     TYEAR = TIME AT END OF CURRENT TIME STEP IN YEARS                 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     PBC(IP)   = SPECIFIED PRESSURE VALUE AT IP(TH) SPECIFIED          
!                 PRESSURE NODE                                         
!     UBC(IP,K) = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY   
!                 INFLOW OCCURRING AT IP(TH) SPECIFIED PRESSURE NODE    
!                 FOR EACH SPECIES                                      
!     IPBC(IP)  = ACTUAL NODE NUMBER OF IP(TH) SPECIFIED PRESSURE NODE  
!                 {WHEN NODE NUMBER I=IPBC(IP) IS NEGATIVE (I<0),       
!                 VALUES MUST BE SPECIFIED FOR PBC AND UBC.}            
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     UBC(IUP,K)  = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE AT     
!                   IU(TH) SPECIFIED CONCENTRATION OR TEMPERATURE NODE  
!                   (WHERE IUP=IU+NPBC)                                 
!                   FOR EACH SPECIES                                    
!     IUBC(IUP,K) = ACTUAL NODE NUMBER OF IU(TH) SPECIFIED              
!                   CONCENTRATION OR TEMPERATURE NODE                   
!                   (WHERE IUP=IU+NPBC)                                 
!                   {WHEN NODE NUMBER I=IUBC(IU) IS NEGATIVE (I<0),     
!                   A VALUE MUST BE SPECIFIED FOR UBC.}                 
!                   FOR EACH SPECIES                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IQSOP(IQP) = NODE NUMBER OF IQP(TH) FLUID SOURCE NODE.            
!                  {WHEN NODE NUMBER I=IQSOP(IQP) IS NEGATIVE (I<0),    
!                  VALUES MUST BE SPECIFIED FOR QIN AND UIN.}           
!     QIN(-I)    = SPECIFIED FLUID SOURCE VALUE AT NODE (-I)            
!     UIN(-I,K)  = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY  
!                  INFLOW OCCURRING AT FLUID SOURCE NODE (-I)           
!                  FOR EACH SPECIES                                     
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IQSOU(IQU,K) = NODE NUMBER OF IQU(TH) ENERGY OR                   
!                    SOLUTE MASS SOURCE NODE                            
!                    {WHEN NODE NUMBER I=IQSOU(IQU) IS NEGATIVE (I<0),  
!                    A VALUE MUST BE SPECIFIED FOR QUIN.}               
!                    FOR EACH SPECIES                                   
!     QUIN(-I,K)  = SPECIFIED ENERGY OR SOLUTE MASS SOURCE VALUE        
!                   AT NODE (-I)                                        
!                   FOR EACH SPECIES                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
!.....ADDITIONAL USEFUL VARIABLES                                       
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     "FUNITS" ARE UNIT NUMBERS FOR INPUT AND OUTPUT FILES              
!         AS ASSIGNED IN THE INPUT FILE, "SUTRA.FIL"                    
!                                                                       
!     X(I), Y(I), AND Z(I) ARE THE X-, Y-, AND Z-COORDINATES OF NODE I  
!     (FOR 2-D PROBLEMS, Z(I) IS THE THICKNESS AT NODE I)               
!                                                                       
!     GRAVX, GRAVY AND GRAVZ ARE THE X-, Y-, AND Z-COMPONENTS OF THE    
!     GRAVITY VECTOR                                                    
!     (FOR 2-D PROBLEMS, GRAVZ = 0)                                     
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
!                                                                       
!.....NSOPI IS ACTUAL NUMBER OF FLUID SOURCE NODES                      
      NSOPI = NSOP - 1 
!.....NSOUI IS ACTUAL NUMBER OF ENERGY OR SOLUTE MASS SOURCE NODES      
      NSOUI = MNSOU - 1 
!                                                                       
!                                                                      
!      READ(5,*) TDLEVEL
!      IF(MOD(IT,354).EQ.0) THEN
!	   REWIND(5)
!      ENDIF
!                                                                       
!                                                                       
!                                                                       
      IF (IPBCT) 50, 240, 240 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (1):  SET TIME-DEPENDENT SPECIFIED PRESSURES OR           
!     CONCENTRATIONS (TEMPERATURES) OF INFLOWS AT SPECIFIED             
!     PRESSURE NODES                                                    
!                                                                       
   50 CONTINUE
      DO 200 IP = 1, NPBC 
         I = SpecifiedPBC(IP)%node
         IF (I) 100, 200, 200 
  100    CONTINUE 
!     NOTE : A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY           
!            TIME STEP IN WHICH PBC( ) CHANGES.                         
!     SpecifiedPBC(IP)%P =  ((          ))                                         
!     DO 150 K=1,NSPE                                                   
! 150 SpecifiedPBC(IP)%U(K) =  ((          ))  

!================================================= CHENGJI 2013-09-02
!	  SpecifiedPBC(IP)%P = 9.8*(1000+SC*714.3)*(TDLEVEL-Y(IABS(I))) 
!==============================!MT - 04112017 - Pressure at seawater boundary ((DRWDT),S)
!=========================================!MT - 04112017 - Seasonal temperature variation
        NOMTH = (TSEC-TSTART)/TIMTH + 1
        NOYEAR = NOMTH / 12.0D0
        TIMEIND = NOMTH - NOYEAR*12.0D0
        IF (BNDSWT.LE.0) THEN
          SWT = SST5*(TIMEIND**5) + SST4*(TIMEIND**4) + SST3*(TIMEIND**3) + SST2*(TIMEIND**2) + SST1*TIMEIND + SST0
!          DRDTSW = 3E-05*(SWT**2) - 0.0064*SWT - 0.0558 !MT - 04112017 - Fitting by excel using NAYAR (2016) data 
!	      SpecifiedPBC(IP)%P = 9.81*(1000 + SC*760 + DRDTSW*(SWT - 4))*(TDLEVEL-Y(IABS(I))) 
        ELSE
          SWT = BNDSWT
!          DRDTSW = 3E-05*(SWT**2) - 0.0064*SWT - 0.0558 !MT - 04112017 - Fitting by excel using NAYAR (2016) data 
!	      SpecifiedPBC(IP)%P = 9.81*(1000 + SC*760 + DRDTSW*(SWT - 4))*(TDLEVEL-Y(IABS(I))) 
        ENDIF
        RHOSWB = (999.9 + 2.034E-2*SWT - 6.162E-3*(SWT**2) + 2.261E-5*(SWT**3) - 4.657E-8*(SWT**4)) + &
                (802*SC - 2.001*SC*SWT + 1.677E-2*SC*(SWT**2) - 3.06E-5*SC*(SWT**3) - 1.613E-5*(SC**2)*(SWT**2))
        RHOFWB = (999.9 + 2.034E-2*FWT - 6.162E-3*(FWT**2) + 2.261E-5*(FWT**3) - 4.657E-8*(FWT**4))
        IF (X(IABS(I)).EQ.X(1)) THEN
            SpecifiedPBC(IP)%P = 9.81*RHOFWB*(FWH-Y(IABS(I))) 
        ELSE
            SpecifiedPBC(IP)%P = 9.81*RHOSWB*(TDLEVEL-Y(IABS(I))) 
        ENDIF
!        WRITE(*,*) IP, I, SWT, FWT, Y(IABS(I)), RHOFWB, SpecifiedPBC(IP)%P
!        PAUSE
!========================================================================================
!	  WRITE(*,*) TDLEVEL, SC
!	  DO 150 K=1,NSPE  
!	  IF(Y(IABS(I)).LE.TDLEVEL) THEN              
!           SpecifiedPBC(IP)%U(K) = SC
!      ELSEIF(Y(IABS(I)).GT.TDLEVEL) THEN
!            SpecifiedPBC(IP)%U(K) = 0
!     ENDIF
!  150 CONTINUE
!============================================= !MT - 06112017 - to consider energy species      
      DO 150 K=1,NSPE
      IF (X(IABS(I)).EQ.X(1)) THEN
        IF(K.EQ.NESP) THEN
            SpecifiedPBC(IP)%U(K) = FWT
        ELSE
            SpecifiedPBC(IP)%U(K) = 0.00
        ENDIF
      ELSE
        IF(K.EQ.NESP) THEN
            IF(Y(IABS(I)).LE.TDLEVEL) THEN              
                SpecifiedPBC(IP)%U(K) = SWT
            ELSEIF(Y(IABS(I)).GT.TDLEVEL) THEN
                SpecifiedPBC(IP)%U(K) = FWT
            ENDIF
        ELSE
            IF(Y(IABS(I)).LE.TDLEVEL) THEN              
                SpecifiedPBC(IP)%U(K) = SC
            ELSEIF(Y(IABS(I)).GT.TDLEVEL) THEN
                SpecifiedPBC(IP)%U(K) = 0
            ENDIF
        ENDIF
     ENDIF
!     WRITE(*,*) IP, I, SWT, K, SpecifiedPBC(IP)%U(K)
!     PAUSE
  150 CONTINUE	
!========================================================================================   	
!=============================================Chengji codes - commented by !MT - 06112017
!      DO 150 K=1,NSPE
!      IF(K.EQ.1) THEN
!        IF(Y(IABS(I)).LE.TDLEVEL) THEN              
!            SpecifiedPBC(IP)%U(K) = SC
!        ELSEIF(Y(IABS(I)).GT.TDLEVEL) THEN
!            SpecifiedPBC(IP)%U(K) = 0
!        ENDIF
!      ELSEIF(K.EQ.2) THEN
!            SpecifiedPBC(IP)%U(K) = 0
!      ENDIF
!  150 CONTINUE	
!========================================================================================   	
      IF(PM1(IABS(I)).GT.0.AND.Y(IABS(I)).GT.TDLEVEL)THEN
	    SpecifiedPBC(IP)%P = 0
	    GNUPCJ(IP)=GNUP  
	      IF(SEEP.EQ.1) THEN
	        IF(Y(IABS(I))+PM1(IABS(I))/10245.GT.SEEPZ) THEN
	         SEEPZ=Y(IABS(I))+PM1(IABS(I))/10245
             SEEPX=X(IABS(I))
            ENDIF
          ENDIF
      ELSEIF(PM1(IABS(I)).GT.0.AND.Y(IABS(I)).LE.TDLEVEL)THEN
	    GNUPCJ(IP) = GNUP
      ELSEIF(PM1(IABS(I)).LT.0.AND.Y(IABS(I)).GT.TDLEVEL)THEN
    	GNUPCJ(IP)= 0
      ELSEIF(PM1(IABS(I)).LT.0.AND.Y(IABS(I)).LE.TDLEVEL)THEN
    	GNUPCJ(IP) = GNUP
      ENDIF
!=====================================================================                 
  200 END DO
  
      IF(SEEP.EQ.1) THEN
        WRITE(8,'(I10,3E16.7)') IT,TDLEVEL,SEEPX,SEEPZ
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  240 IF (IUBCT) 250, 440, 440 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (2):  SET TIME-DEPENDENT SPECIFIED                        
!     CONCENTRATIONS (TEMPERATURES)                                     
!                                                                       
  250 CONTINUE 
      DO 400 K = 1, NSPE 
         DO 350 IU = 1, NUBC (K) 
            IUP = NPBC + IU 
            I = MultiSpeciesBC(K)%SpecifiedU(IU)%node
            IF (I) 300, 400, 400 
  300       CONTINUE 
!       NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY TIME STEP        
!              IN WHICH UBC( ) CHANGES.  IN ADDITION, IF FLUID          
!              PROPERTIES ARE SENSITIVE TO 'U' THEN A FLOW SOLUTION     
!              MUST OCCUR AS WELL
!        IF(K.EQ.1) THEN 
!             MultiSpeciesBC(K)%SpecifiedU(IU)%U=0
!        ELSEIF(K.EQ.2) THEN
!           IF(IT.LE.7200) THEN                                     
!             MultiSpeciesBC(K)%SpecifiedU(IU)%U=0.035
!           ELSE
!             MultiSpeciesBC(K)%SpecifiedU(IU)%U=0
!           ENDIF
!        ENDIF
  350    END DO 
  400 END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  440 IF (IQSOPT) 450, 640, 640 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (3):  SET TIME-DEPENDENT FLUID SOURCES/SINKS,             
!      OR CONCENTRATIONS (TEMPERATURES) OF SOURCE FLUID                 
!                                                                       
  450 CONTINUE 
      DO 600 IQP = 1, NSOPI 
         I = IQSOP (IQP) 
         IF (I) 500, 600, 600 
  500    CONTINUE 
!     NOTE : A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY           
!            TIME STEP IN WHICH QIN( ) CHANGES.                         
!     QIN(-I) =   ((           ))                                       
!     NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY                    
!            TIME STEP IN WHICH UIN( ) CHANGES.                         
!     DO 550 K=1,NSPE                                                   
!       UIN(-I,K) =   ((           ))                                   
! 550 CONTINUE

      IF(TSEC.LE.10*86400) THEN
        QIN(-I)= 0.0124
      ELSE
        QIN(-I)= 0
      ENDIF
!      
      DO 550 K=1,NSPE
      IF(K.EQ.1) THEN                                                  
            UIN(-I,K)=0
      ELSEIF(K.EQ.2) THEN
        IF(TSEC.LE.10*86400) THEN
            UIN(-I,K)=0.1
        ELSE
            UIN(-I,K)=0
        ENDIF
      ENDIF
  550 CONTINUE

                                                          
  600 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  640 IF (IQSOUT) 650, 840, 840 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (4):  SET TIME-DEPENDENT SOURCES/SINKS                    
!     OF SOLUTE MASS OR ENERGY                                          
!                                                                       
  650 CONTINUE 
      DO 800 K = 1, NSPE 
         DO 750 IQU = 1, NSOU (K) 
            I = IQSOU (IQU, K) 
            IF (I) 700, 800, 800 
  700       CONTINUE 
!       NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY                  
!              TIME STEP IN WHICH QUIN( ) CHANGES.                      
!       QUIN(-I,K) =   ((           ))                                  
  750    END DO 
  800 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  840 CONTINUE 
!                                                                       
      RETURN 
      END SUBROUTINE BCTIME 

!                                                                       
!     SUBROUTINE        U  N  S  A  T              SUTRA VERSION 2D3D.1 
!                                                                       
! *** PURPOSE :                                                         
! ***  USER-PROGRAMMED SUBROUTINE GIVING:                               
! ***  (1)  SATURATION AS A FUNCTION OF PRESSURE ( SW(PRES) )           
! ***  (2)  DERIVATIVE OF SATURATION WITH RESPECT TO PRESSURE           
! ***       AS A FUNCTION OF EITHER PRESSURE OR SATURATION              
! ***       ( DSWDP(PRES), OR DSWDP(SW) )                               
! ***  (3)  RELATIVE PERMEABILITY AS A FUNCTION OF EITHER               
! ***       PRESSURE OR SATURATION ( REL(PRES) OR RELK(SW) )            
! ***                                                                   
! ***  CODE BETWEEN DASHED LINES MUST BE REPLACED TO GIVE THE           
! ***  PARTICULAR UNSATURATED RELATIONSHIPS DESIRED.                    
! ***                                                                   
! ***  DIFFERENT FUNCTIONS MAY BE GIVEN FOR EACH REGION OF THE MESH.    
! ***  REGIONS ARE SPECIFIED BY BOTH NODE NUMBER AND ELEMENT NUMBER     
! ***  IN INPUT DATA FILE FOR UNIT fINP.                                  
!                                                                       
      SUBROUTINE UNSAT (SW, DSWDP, RELK, PRES, KREG) 
      
      USE CONTRL
      USE M_SWCC
      USE SutraMSPrecision

      IMPLICIT NONE
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
!     E X A M P L E   C O D I N G   FOR                                 
!     MESH WITH TWO REGIONS OF UNSATURATED PROPERTIES USING             
!     THREE PARAMETER-UNSATURATED FLOW RELATIONSHIPS OF                 
!     VAN GENUCHTEN(1980)                                               
!        RESIDUAL SATURATION, SWRES, GIVEN IN UNITS {L**0}              
!        PARAMETER, AA, GIVEN IN INVERSE PRESSURE UNITS {m*(s**2)/kg}   
!        PARAMETER, VN, GIVEN IN UNITS {L**0}                           
!                                                                       
      INTEGER (I4B) :: &
        KREG
      REAL (DP) :: &
        SW, DSWDP, RELK, PRES

      !LOCAL VARIABLES
      REAL (DP) :: &
        SWRES, AA, VN, SWRM1, AAPVN, VNF, AAPVNN, DNUM, DNOM, SWSTAR 
!      REAL (DP) :: &
!        SWRES1, SWRES2, AA1, AA2, VN1, VN2 
!                                                                       
!     DATA FOR REGION 1:                                                
!      DATA SWRES1 / 0.099E0 /, AA1 / 1.48E-3 /, VN1 / 2.68E0 / 
!      SAVE SWRES1, AA1, VN1 
!     DATA FOR REGION 2:                                                
!      DATA SWRES2 / 0.099E0 /, AA2 / 1.48E-3 /, VN2 / 2.68E0 / 
!      SAVE SWRES2, AA2, VN2 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
!                                                                       
! *** BECAUSE THIS ROUTINE IS CALLED OFTEN FOR UNSATURATED FLOW RUNS,   
! *** EXECUTION TIME MAY BE SAVED BY CAREFUL CODING OF DESIRED          
! *** RELATIONSHIPS USING ONLY INTEGER AND SINGLE PRECISION VARIABLES!  
! *** RESULTS OF THE CALCULATIONS MUST THEN BE PLACED INTO DOUBLE       
! *** PRECISION VARIABLES SW, DSWDP AND RELK BEFORE LEAVING             
! *** THIS SUBROUTINE.                                                  
!                                                                       
!                                                                       
!***********************************************************************
      ! KONG
!      GOTO 1800


!***********************************************************************
!                                                                       
!     SET PARAMETERS FOR CURRENT REGION, KREG                           
      GOTO (10, 20), KREG 
   10 SWRES = REAL(SWRES1)
      AA = REAL(AA1)/9800
      VN = REAL(VN1) 
      GOTO 100 
   20 SWRES = REAL(SWRES2)
      AA = REAL(AA2)/9800
      VN = REAL(VN2) 
  100 CONTINUE 
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!.....SECTION (1):                                                      
!     SW VS. PRES   (VALUE CALCULATED ON EACH CALL TO UNSAT)            
!     CODING MUST GIVE A VALUE TO SATURATION, SW.                       
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     THREE PARAMETER MODEL OF VAN GENUCHTEN(1980)                      
      SWRM1 = 1.E0 - SWRES 
      AAPVN = 1.E0 + (AA * ( - PRES) ) **VN 
      VNF = (VN - 1.E0) / VN 
      AAPVNN = AAPVN**VNF 
      SW = DBLE (SWRES + SWRM1 / AAPVNN) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF (IUNSAT - 2) 600, 1200, 1800 
!***********************************************************************
!***********************************************************************
!.....SECTION (2):                                                      
!     DSWDP VS. PRES, OR DSWDP VS. SW   (CALCULATED ONLY WHEN IUNSAT=1) 
!     CODING MUST GIVE A VALUE TO DERIVATIVE OF SATURATION WITH         
!     RESPECT TO PRESSURE, DSWDP.                                       
!                                                                       
  600 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      DNUM = AA * (VN - 1.E0) * SWRM1 * (AA * ( - PRES) ) ** (VN - 1.E0) 
      DNOM = AAPVN * AAPVNN 
      DSWDP = DBLE (DNUM / DNOM) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      GOTO 1800 
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!.....SECTION (3):                                                      
!     RELK VS. P, OR RELK VS. SW   (CALCULATED ONLY WHEN IUNSAT=2)      
!     CODING MUST GIVE A VALUE TO RELATIVE PERMEABILITY, RELK.          
!                                                                       
 1200 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     GENERAL RELATIVE PERMEABILITY MODEL FROM VAN GENUCHTEN(1980)      
      SWSTAR = (SW - SWRES) / SWRM1 
      RELK = DBLE (SQRT (SWSTAR) * (1.E0 - (1.E0 - SWSTAR** (1.E0 / VNF)) ** (VNF) ) **2)                                                 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
 1800 RETURN
	  
      END
