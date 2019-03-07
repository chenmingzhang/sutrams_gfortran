!     SUBROUTINE        B  U  D  G  E  T       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALCULATE AND OUTPUT FLUID MASS AND SOLUTE MASS OR            
! ***  ENERGY BUDGETS.                                                  
!                                                                       
      SUBROUTINE BUDGET (ML, IBCT, VOL, SW, DSWDP, RHO, QIN, PVEC, &
                         PM1, PBC, QPLITR, IPBC, IQSOP, UVEC, UM1, UM2, UIN, QUIN,    &
                         IQSOU, UBC, IUBC, CS1, CS2, CS3, SL, SR, GNUPCJ) ! CHENGJI 2013-09-02
      USE PARAMS 
      USE CONTRL 
      USE MODSOR 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE SutraStorage, ONLY : SpecifiedPBC, &
                               MultiSpeciesBC
      use SutraMSPrecision
      USE SutraZoneModule
      USE MSErrorHandler

      USE SOLVI ! Chengji 2015-04-14

      implicit none
	  
!============================================= CHENGJI 2013-09-02	
      REAL(DP):: GNUPCJ  
      DIMENSION GNUPCJ(NBCN)
!================================================================  
	  
      CHARACTER(13) UNAME (2) 
      integer (I4B) :: &
        IBCT, ML
      integer (I4B) :: &
        IQSOP (NSOP)
      integer (I4B) :: &
        IQSOU (MNSOU, NSPE)                                               
      integer (I4B) :: &
        IPBC (NBCN), IUBC (NBCN, NSPE)
      real (kind=8) :: &
        QIN (NN), UIN (NN, NSPE)
      real (DP) :: &
        QUIN (NN, NSPE)
      real (DP) :: &
        UBC (NBCN, NSPE),       &
        QPLITR (NBCN), PBC (NBCN)                                         
      real (DP) :: &
        VOL (NN), PVEC (NN), UVEC (NN, NSPE), SW (NN),&
        DSWDP (NN), RHO (NN), PM1 (NN), UM1 (NN, NSPE), UM2 (NN,NSPE), &
        CS1 (NN, NSPE), CS2 (NN), CS3 (NN), SL (NN), SR (NN)
 !     real (DP) :: &    ! COMMENTED BY CHENGJI
 !       STUNT (NSPE), STUPT (NSPE), STUTOT (NSPE) 
      !locals 
      integer (I4B) :: &
        I, K, &
        IP, IQU, IU, IUP, &
        MN
      integer (I4B) :: &
        INEGCT, IQP, NSOPI, NSOUI
      integer (I4B) :: imap
      real (DP) :: &
!        CST, CWT, &    !MT 26102017 - COMMENTED OUT
        CST, CWT, &
        RELK, &
        DNSTOT, FLDTOT, P0FTOT, P1FTOT, P0STOT, P1STOT, &
        QIUTOT, QQUTOT, QPUTOT, QULTOT, SLDTOT, &
        DUDT, EPRSV, ESRV, &
        QU, QPU
	  
      REAL(DP):: STPPOS, STPNEG, STUPOS(NSPE), STUNEG(NSPE), QINPOS, QINNEG, &
				 STPTOT, STUTOT, STFPOS, STFNEG, STFTOT, QINTOT, &
				 QPLPOS, QPLNEG, QPLTOT,      &
				 QFFPOS, QFFNEG, QFFTOT,      &
				 TERM,                        & ! Chengji 2013-09-24
                 TERM11,TERMSPE1,TERMSPE2,    & ! Chengji 2015-04-14
                 ACTFMB,ERFMBA,ERFMBR,        & ! Chengji 2015-08-06
                 QPSOLUTE                       ! Chengji 2016-01-05

      DATA UNAME (1) / 'CONCENTRATION' / , UNAME (2) / ' TEMPERATURE ' / 
      SAVE UNAME 
!                                                                       
!.....SAVE VALUES FOR ENERGY TRANSPORT                                                                                                        

      CST = CS
	  CWT = CW                   
           

!                                                                       
      MN = 2 
      IF (IUNSAT.NE.0) IUNSAT = 1 
      IF (ME.EQ. - 1) MN = 1 
      WRITE (fLST, 10) 
   10 FORMAT(1H1) 
!.....SET UNSATURATED FLOW PARAMETERS, SW(I) AND DSWDP(I)               
      IF (IUNSAT - 1) 40, 20, 40 
   20 DO 30 I = 1, NN 
         IF (PVEC(I) ) 25, 27, 27 
   25    CALL UNSAT (SW(I), DSWDP(I), RELK, PVEC(I), NodeMap(I) ) 
         GOTO 30 
   27    SW(I)    = 1.0D0 
         DSWDP(I) = 0.0D0 
   30 END DO 
!                                                                       
!.....CALCULATE COMPONENTS OF FLUID MASS BUDGET                         
   40 IF (ML - 1) 50, 50, 1000 
   50 CONTINUE 
 
!====================================== Added according to SUTRA, Chengji 2013-09-24
      STPPOS = 0.D0                                            
      STPNEG = 0.D0                                            
      DO K=1,NSPE
	  STUPOS(K) = 0.D0   ! Influx for Kth species                                        
      STUNEG(K) = 0.D0   ! Outflux for Kth species
      ENDDO
      QINPOS = 0.D0                                            
      QINNEG = 0.D0
      
      STPTOT = 0.D0
      STUTOT = 0.D0
      STFPOS = 0.D0
      STFNEG = 0.D0
      STFTOT = 0.D0
      QINTOT = 0.D0
      
	  
!***************************************************************************************** 
!***************************************************************************************** 
!MT - 12112017 - Commented out - Don't need Boundary water flux for SI problem
!	  DO 111 I=1,NN  ! Chengji 2015-04-14
!=========================================================SET DRWDU A FUNCTION OF SALINITY
!         DRWDU(NESP) = DRWDT(I)     !MT 24102017 
!=============================================================================!MT 26102017
!	  IF(MOD(I,NN1).EQ.1) THEN                                       
!      TERM11 = (1-ISSFLO/2)*RHO(I)*VOL(I)*      &
!           (SW(I)*NodeData(NodeMap(i))%sop+NodeData(NodeMap(i))%por*DSWDP(I))*      &
!           (PVEC(I)-PM1(I))/DELTP

!     TERMSPE1 = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SW(I)*DRWDU(1)*VOL(I)*      &
!     (UM1(I,1)-UM2(I,1))/DLTUM1                                
!

!     TERMSPE2 = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SW(I)*DRWDU(2)*VOL(I)*      &
!     (UM1(I,2)-UM2(I,2))/DLTUM1
           
!      WRITE(200988,0307) I,TERM11,TERMSPE1,TERMSPE2
!      ENDIF	  

!  111 CONTINUE
! 0307 FORMAT(I8,4E15.7)
!***************************************************************************************** 
!*****************************************************************************************   
	 
!MT - 12112017 - Commented out to test new implementation of Density treatment  

!      DO 100 I=1,NN 
!=========================================================SET DRWDU A FUNCTION OF SALINITY
!         DRWDU(NESP) = DRWDT(I)     !MT 24102017 
!=============================================================================!MT 26102017                                    
!      TERM = (1-ISSFLO/2)*RHO(I)*VOL(I)*      &
!           (SW(I)*NodeData(NodeMap(i))%sop+NodeData(NodeMap(i))%por*DSWDP(I))*      &
!           (PVEC(I)-PM1(I))/DELTP
!      STPPOS = STPPOS + MAX(0D0, TERM) ! Influx caused by pressure change             
!      STPNEG = STPNEG + MIN(0D0, TERM) ! Outflux caused by pressure change

!     DO K=1,NSPE
!     TERM = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SW(I)*DRWDU(K)*VOL(I)*      &
!     (UM1(I,K)-UM2(I,K))/DLTUM1                                
!     STUPOS(K) = STUPOS(K) + MAX(0D0, TERM)  ! Influx caused by concentration change
!     STUNEG(K) = STUNEG(K) + MIN(0D0, TERM)  ! Outflux caused by concentration change
!     ENDDO
!      
!	  TERM = QIN(I)                                           


!      QINPOS = QINPOS + MAX(0D0, TERM)  ! Influx by flux boundary                       
!      QINNEG = QINNEG + MIN(0D0, TERM)  ! Outflux by flux boundary
!  100 CONTINUE
!  

!      STPTOT = STPPOS + STPNEG  ! Net flux caused by pressure change                  
!      DO K=1,NSPE
!	  STUTOT = STUTOT + STUPOS(K) + STUNEG(K) ! Net flux caused by concentration change         
!      STFPOS = STPPOS + STUPOS(K) ! Influx caused by pressure change and concentration change  
!      STFNEG = STPNEG + STUNEG(K) ! Outflux caused by pressure change and concentration change
!      ENDDO
!      STFTOT = STPTOT + STUTOT  ! Net flux caused by pressure change and concentration change                               
!      QINTOT = QINPOS + QINNEG  ! Net flux caused by flux boundary
!
      QPLPOS = 0.D0                       
      QPLNEG = 0.D0 
      QPLTOT = 0.D0
             
      DO 200 IP=1,NPBC                   
      I=IABS(SpecifiedPBC(IP)%node)                   
      TERM = GNUPCJ(IP)*(SpecifiedPBC(IP)%P-PVEC(I))
      QPLITR(IP) = GNUPCJ(IP)*(SpecifiedPBC(IP)%P-PVEC(I))
      QPLPOS = QPLPOS + MAX(0D0, TERM)  ! Influx caused by pressure-precribed boundary nodes
      QPLNEG = QPLNEG + MIN(0D0, TERM)  ! Outflux caused by pressure-precribed boundary nodes
  200 CONTINUE
                       
      QPLTOT = QPLPOS + QPLNEG   ! Net flux caused by pressure-prescribed boundary nodes        
      QFFPOS = QINPOS + QPLPOS           
      QFFNEG = QINNEG + QPLNEG           
      QFFTOT = QINTOT + QPLTOT 
      
      ACTFMB = 5D-1*(STFPOS - STFNEG + QFFPOS - QFFNEG)
      ERFMBA = STFTOT - QFFTOT
      IF (ACTFMB.NE.0D0) THEN 
        ERFMBR = 1D2*ERFMBA/ACTFMB
      ENDIF
!      WRITE(1,'(I6,6E15.7)') IT,STFPOS,STFNEG,STFTOT,QFFPOS,QFFNEG,QFFTOT
!      WRITE(1,'(I6,15E15.7)') IT,STPPOS,STUPOS(1),STUPOS(2),STPNEG,STUNEG(1),STUNEG(2),   &
!                             STFTOT,QFFPOS,QFFNEG,QFFTOT,ERFMBR
      WRITE(1,'(I6,15E15.7)') IT,QPLPOS, QPLNEG, QPLTOT, QFFPOS,QFFNEG,QFFTOT,ERFMBR
!======================================================================================
	  
!      STPTOT = 0.D0        !!!!!!!!!!!!!!!!!!!!!!!! Original in SUTRA_MS, commented by Chengji 2013-09-24
!      STUTOT = 0.D0 
!      STUNT = 0.D0 
!      STUPT = 0.D0 
!      QINTOT = 0.D0 
!      DO 100 I = 1, NN 
!         STPTOT = STPTOT + (1 - ISSFLO / 2) * RHO(I) * VOL(I) *       &
!           (SW(I) * NodeData(NodeMap(i))%sop + NodeData(NodeMap(i))%por * DSWDP(I) ) * &
!           (PVEC(I) - PM1(I)) / DELTP
!======= STPTOT: total rate of change in stored fluid mass in the region due to pressure changes
		   
!        DO K = 1, NSPE 
!           STUT = (1 - ISSFLO / 2) * NodeData(NodeMap(i))%por * SW(I) * DRWDU (K) * VOL(I) * &
!           ( UM1(I, K) - UM2(I, K) ) / DLTUM1                       
!           STUTOT (K) = STUTOT (K) + STUT 
!======= STUTOT(K): total rate of change in stored fluid mass in the region due to concentration changes for
!======= species K

!           IF (STUT.LT.0.0D0) STUNT (K) = STUNT (K) + STUT  ! Negative rate of change in stored fluid due U change
!           IF (STUT.GT.0.0D0) STUPT (K) = STUPT (K) + STUT  ! Positive rate of change in stored fluid due U c
!         ENDDO 
!         QINTOT = QINTOT + QIN(I) 
!======= QINTOT: total rate of fluid mass change due to all sources and sinks of fluid mass, as well as to
!======= specified fluxes across boundaries
!  100 END DO 
                                                                     
!      QPLTOT = 0.D0     !!!!!!!!!!!!!!!!!!!!!!!! Original in SUTRA_MS, commented by Chengji 2013-09-24
!      DO 200 IP = 1, NPBC 
!         I = IABS (IPBC (IP) ) 
!         QPLITR (IP) = GNUP * (PBC (IP) - PVEC(I) ) 
!         I = IABS ( SpecifiedPBC(IP)%node ) 
!         QPLITR (IP) = GNUPCJ(IP) * ( SpecifiedPBC(IP)%P - PVEC(I) ) 
!         QPLTOT = QPLTOT + QPLITR (IP)
!======= QPLTOT: total rate of fluid mass change in the entire region due to inflows and outflows at all 
!======= specified pressure nodes
!  200 END DO 
!======= CHENGJI 2013-09-03, output total fluid mass flux to '.dat' file
!                                                                       
!.....OUTPUT FLUID MASS BUDGET                                          
!      WRITE (fLST, 300) IT, STPTOT, QINTOT, QPLTOT 
!  300 FORMAT(//11X,'F L U I D   M A S S   B U D G E T      AFTER TIME', &
!     &   ' STEP ',I5,',     IN (MASS/SECOND)'///11X,1PD15.7,5X,         &
!     &   'RATE OF CHANGE IN TOTAL STORED FLUID DUE TO PRESSURE CHANGE', &
!     &   ', INCREASE(+)/DECREASE(-)',                                   &
!     &   //11X,1PD15.7,5X,'TOTAL OF FLUID SOURCES AND SINKS, ',         &
!     &   'NET INFLOW(+)/NET OUTFLOW(-)'/11X,1PD15.7,5X,                 &
!     &   'TOTAL OF FLUID FLOWS AT POINTS OF SPECIFIED PRESSURE, ',      &
!     &   'NET INFLOW(+)/NET OUTFLOW(-)')                                
!                                                                       
!.....OUTPUT FLUID MASS BUDGET CHANGES DUE TO CHANGES                   
!     IN CONCENTRATION                                                  
!      DO 305 K = 1, NSPE
!  305 WRITE (fLST, 310) K, trim(adjustl(SPNAME(K))), STUNT (K), STUPT (K), STUTOT (K)
	  
!  310 FORMAT(/11X,'[SPECIES ',I3,'] - ',A,/11X,1PD15.7,5X,            &
!     &   'NEGATIVE RATE OF CHANGE IN STORED FLUID DUE U CHANGE',        &
!     &   /11X,1PD15.7,5X,'POSITIVE RATE OF CHANGE IN STORED ',          &
!     &   'FLUID DUE U CHANGE',/11X,1PD15.7,5X                           &
!     &   'RATE OF CHANGE IN TOTAL STORED FLUID DUE TO U CHANGE',        &
!     &   ', INCREASE(+)/DECREASE(-)')                                   
!                                                                       
      IF (IBCT.EQ.4) GOTO 600 
      NSOPI = NSOP - 1 
      INEGCT = 0 
      DO 500 IQP = 1, NSOPI 
         I = IQSOP (IQP) 
         IF (I) 325, 500, 500 
  325    INEGCT = INEGCT + 1 
         IF (INEGCT.EQ.1) WRITE (fLST, 350) 
  350 FORMAT(///22X,'TIME-DEPENDENT FLUID SOURCES OR SINKS'//22X,       &
     &   ' NODE',5X,'INFLOW(+)/OUTFLOW(-)'/37X,'  (MASS/SECOND)'//)     
         WRITE (fLST, 450) - I, QIN ( - I) 
  450 FORMAT(18X,I9,10X,1PD15.7) 
  500 END DO 
!                                                                       
  600 IF (NPBC.EQ.0) GOTO 800 
!      WRITE (fLST, 650)   ! CHENGJI 2013-09-03
!  650 FORMAT(///22X,'FLUID SOURCES OR SINKS DUE TO SPECIFIED PRESSURES',&
!     &   //22X,' NODE',5X,'INFLOW(+)/OUTFLOW(-)'/37X,'  (MASS/SECOND)'/)
      DO 700 IP = 1, NPBC 
!         I = IABS ( IPBC (IP) ) 
         I = IABS ( SpecifiedPBC(IP)%node ) 
!         WRITE (fLST, 450) I, QPLITR (IP)
         WRITE (2,'(I8,4E15.7)') I, QPLITR (IP)  ! CHENGJI 2013-09-13, write node-wise fluid mass to the '.dat' file
  700 END DO 
!********************************************************************************************
!********************************************************************************************












!********************************************************************************************
!********************************************************************************************
!.....CALCULATE COMPONENTS OF ENERGY OR SOLUTE MASS BUDGET              
  800 IF (ML - 1) 1000, 5500, 1000 
 1000 CONTINUE 
!                                                                       
      DO 5000 K = 1, NSPE 
!.....GLOBAL COUNTER FOR SPECIES - DEFINED IN MODULE PARAMS             
         KSP = K 
!                                                                       
!.....WRITE SPECIES INFORMATION TO OUTPUT FILE                          
         WRITE (fLST, 1010) K, trim(adjustl(SPNAME(K)))
 1010 FORMAT(//11X,'[SPECIES ',I3,'] - ',A) 
!                                                                       
!.....ZERO VARIABLE FOR MASS BALANCE CALCULATIONS                        
         FLDTOT = 0.D0 
         SLDTOT = 0.D0 
         DNSTOT = 0.D0 
         P1FTOT = 0.D0 
         P1STOT = 0.D0 
         P0FTOT = 0.D0 
         P0STOT = 0.D0 
         QQUTOT = 0.D0 
         QIUTOT = 0.D0 
!.....SET ADSORPTION PARAMETERS                                         
!      IF (ME.EQ. - 1.AND.ADSMOD (K) .NE.'NONE      ') CALL ADSORB (CS1, CS2, CS3, SL, SR, UVEC)                                           
      IF (ME.EQ. - 1.AND.ADSMOD (K) .NE.'NONE      ') CALL ADSORB ()
!                                                                       
!.....SET APPROPRIATE VALUES FOR CS AND CW                              
         IF (K.EQ.NESP) THEN 
            CS = CST 
!            CW = CWT  
         ELSE 
            CS = 0.0D0 
!            CW = 1.0D0 
         ENDIF 
!                                                                       
         DO 1300 I = 1, NN
            DRWDU(NESP) = DRWDT(I)     !MT 24102017  
!			IF (K.EQ.NESP) THEN
!				CW = CWT(I)     			!MT 26102017 	
!			ELSE
!				CW = 1.0
!		    ENDIF
            imap = NodeMap(i)
            ESRV = NodeData(imap)%por * SW(I) * RHO(I) * VOL(I) 
            EPRSV = (1.D0 - NodeData(imap)%por ) * NodeData(imap)%rhos * VOL(I) 
            DUDT = (1 - ISSTRA) * ( UVEC(I, K) - UM1(I, K) ) / DELTU 
            FLDTOT = FLDTOT + ESRV * CW * DUDT 
            SLDTOT = SLDTOT + EPRSV * CS1(I, K) * DUDT            
!======= FLDTOT + SLDTOT = total rate of change in stored solute mass in the region due to change in concentration
!======= FLDTOT + SLDTOT equals to (5.15a) on page 112 of the manual		
!MT - 12112017 - Commented out to test new implementation of Density treatment 

!            DNSTOT = DNSTOT + CW * UVEC(I, K) * (1 - ISSFLO / 2) *      &
!              VOL(I) * ( RHO(I) * (SW(I) * NodeData(imap)%sop + NodeData(NodeMap(i))%por * DSWDP(I) ) * &
!              ( PVEC(I) - PM1(I) ) / DELTP + NodeData(imap)%por * SW(I) *  &
!              DRWDU (K) * ( UM1(I, K) - UM2(I, K) ) / DLTUM1) 
!=============================================================================
!======= DNSTOT = total rate of change in stored solute mass in the region due to change in stored fluid mass
!======= DNSTOT equals to (5.15b) on page 112 of the manual          
!MT - 12112017 - Commented out to test new implementation of Density treatment
!         IF(K.EQ.NSPE) THEN
!          IF(MOD(I,NN1).EQ.1) THEN
!             WRITE(20098,'(2I8,10E15.7)') K, I, DUDT, (PVEC(I)-PM1(I))/DELTP, (UM1(I,K)-UM2(I,K))/DLTUM1,          &         
!                                                ESRV*CW*DUDT, EPRSV*CS1(I,K)*DUDT, CW*UVEC(I,K)*(1-ISSFLO/2)*      &
!                                 VOL(I)*(RHO(I)*(SW(I)*NodeData(imap)%sop+NodeData(NodeMap(i))%por*DSWDP(I))*      &
!                                                             (PVEC(I)-PM1(I))/DELTP+NodeData(imap)%por*SW(I)*      &
!                                                                        DRWDU(K)*(UM1(I,K)-UM2(I,K))/DLTUM1)
!          ENDIF
!         ENDIF
!=============================================================================
!		  
!            P1FTOT = P1FTOT + ESRV * PRODF1 (K) * UVEC(I, K) 
!            P1STOT = P1STOT + EPRSV * PRODS1 (K) * ( SL(I) * UVEC(I, K) + SR(I) )                                                  
!            P0FTOT = P0FTOT + ESRV * PRODF0 (K) 
!            P0STOT = P0STOT + EPRSV * PRODS0 (K) 
            P1FTOT = P1FTOT + ESRV  * ProdSorp(imap)%prodf1(K) * UVEC(I, K) 
!======= P1FTOT is the total rate of first-order solute mass production in the fluid, (5.16a) on page 113			
            P1STOT = P1STOT + EPRSV * ProdSorp(imap)%prods1(K) * ( SL(I) * UVEC(I, K) + SR(I) )  
!======= P1STOT is the total rate of first-order adsorbate production in the fluid, (5.16b) on page 113			
            P0FTOT = P0FTOT + ESRV  * ProdSorp(imap)%prodf0(K) 
            P0STOT = P0STOT + EPRSV * ProdSorp(imap)%prods0(K) 
!======= P0FTOT + P0STOT is the zero-order production of solute and adsorbate mass production in the fluid and solid matrix
!======= (5.17) on page 113
            QQUTOT = QQUTOT + QUIN(I, K) 
!======= QQUTOT is the diffusive-dispersive source of solute mass, (5.20) on page 113
            IF ( QIN(I) ) 1200, 1200, 1250 
 1200       QIUTOT = QIUTOT + QIN(I) * CW * UVEC(I, K) 
            GOTO 1300 
 1250       QIUTOT = QIUTOT + QIN(I) * CW * UIN(I, K)
!======= QIUTOT is the solute mass change at fluid source node, (5.18) on page 113 
 1300    END DO 
!                                                                       
!.....MASS CHANGES DUE TO SPECIFIED PRESSURES                           
         QPUTOT = 0.
         
!         DO 1500 IP = 1, NPBC 
!            IF (QPLITR (IP) ) 1400, 1400, 1450 
! 1400       I = IABS ( SpecifiedPBC(IP)%node ) 
!            QPUTOT = QPUTOT + QPLITR (IP) * CW * UVEC(I, K)
!            QPSOLUTE = QPLITR (IP) * CW * UVEC(I, K)
!            GOTO 1500 
! 1450       QPUTOT = QPUTOT + QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
!            QPSOLUTE = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
 !1500    END DO      
      
         DO 1500 IP = 1, NPBC 
          IF (QPLITR(IP).LE.0) THEN 
            I = IABS ( SpecifiedPBC(IP)%node ) 
            QPUTOT = QPUTOT + QPLITR (IP) * CW * UVEC(I, K)
            QPSOLUTE = QPLITR (IP) * CW * UVEC(I, K)
          ELSEIF (QPLITR(IP).GT.0) THEN 
            QPUTOT = QPUTOT + QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
            QPSOLUTE = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
          ENDIF
            WRITE(4,'(2I8,4E15.7)') K, IABS(SpecifiedPBC(IP)%node), QPSOLUTE
    1500 END DO         
!======= QPUTOT is the solute mass change at specified pressure node, (5.19) on page 113  
 
!                                                                       
!.....MASS CHANGES DUE TO SPECIFIED CONCENTRATIONS                      
         QULTOT = 0.D0 
         IF (NUBC (K) .EQ.0) GOTO 1540 
         DO 1510 IU = 1, NUBC (K) 
            IUP = IU + NPBC 
            I = IABS ( MultiSpeciesBC(K)%SpecifiedU(IU)%node ) 
            QPLITR (IUP) = GNUU (K) * ( MultiSpeciesBC(K)%SpecifiedU(IU)%U - UVEC(I, K) ) 
            QULTOT = QULTOT + QPLITR (IUP) 
!======= QULTOT is the solute mass due to specified concentration conditions, (5.21) on page 114
 1510    END DO 
!======= FLDTOT+SLDTOT+DNSTOT should be close to P1FTOT+P1STOT+P0FTOT+P0STOT+QIUTOT+QPUTOT+QQUTOT+QULTOT                                                                      

 1540    IF (K.EQ.NESP) GOTO 1615 
		
!======= CHENGJI 2013-09-03, output total solute mass to 'dat' file		
		WRITE(3,'(I6,3E15.7)') IT,FLDTOT+SLDTOT+DNSTOT,P1FTOT+P1STOT+P0FTOT+P0STOT+QIUTOT+QPUTOT+QQUTOT+QULTOT
!                                                                       
!.....OUTPUT SOLUTE MASS BUDGET                                         
 1550    WRITE (fLST, 1600) IT, FLDTOT, SLDTOT, DNSTOT, P1FTOT, P1STOT,   &
         P0FTOT, P0STOT, QIUTOT, QPUTOT, QQUTOT, QULTOT                 
 1600 FORMAT(//11X,'S O L U T E   B U D G E T      AFTER TIME STEP ',I5,&
     &   ',   IN (SOLUTE MASS/SECOND)'///11X,1PD15.7,5X,'NET RATE OF ', &
     &   'INCREASE(+)/DECREASE(-) OF SOLUTE DUE TO CONCENTRATION CHANGE'&
     &   /11X,1PD15.7,5X,'NET RATE OF INCREASE(+)/DECREASE(-) OF ',     &
     &   'ADSORBATE'/11X,1PD15.7,5X,'NET RATE OF INCREASE(+)/',         &
     &   'DECREASE(-) OF SOLUTE DUE TO CHANGE IN MASS OF FLUID'//11X,   &
     &   1PD15.7,5X,'NET FIRST-ORDER PRODUCTION(+)/DECAY(-) OF SOLUTE'  &
     &   /11X,1PD15.7,5X,'NET FIRST-ORDER PRODUCTION(+)/DECAY(-) OF ',  &
     &   'ADSORBATE'/11X,1PD15.7,5X,'NET ZERO-ORDER PRODUCTION(+)/',    &
     &   'DECAY(-) OF SOLUTE'/11X,1PD15.7,5X,'NET ZERO-ORDER ',         &
     &   'PRODUCTION(+)/DECAY(-) OF ADSORBATE'/11X,1PD15.7,5X,          &
     &   'NET GAIN(+)/LOSS(-) OF SOLUTE THROUGH FLUID SOURCES AND SINKS'&
     &   /11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF SOLUTE THROUGH ',      &
     &   'INFLOWS OR OUTFLOWS AT POINTS OF SPECIFIED PRESSURE'          &
     &   /11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF SOLUTE THROUGH ',      &
     &   'SOLUTE SOURCES AND SINKS'/11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-)'&
     &  ,' OF SOLUTE AT POINTS OF SPECIFIED CONCENTRATION')             
         GOTO 1645 
!                                                                       
!.....OUTPUT ENERGY BUDGET                                              
 1615    WRITE (fLST, 1635) IT, FLDTOT, SLDTOT, DNSTOT, P0FTOT, P0STOT,   &
         QIUTOT, QPUTOT, QQUTOT, QULTOT                                 
 1635 FORMAT(//11X,'E N E R G Y   B U D G E T      AFTER TIME STEP ',I5,&
     &   ',   IN (ENERGY/SECOND)'///11X,1PD15.7,5X,'NET RATE OF ',      &
     &   'INCREASE(+)/DECREASE(-) OF ENERGY IN FLUID DUE TO TEMPERATURE'&
     & ,' CHANGE'/11X,1PD15.7,5X,'NET RATE OF INCREASE(+)/DECREASE(-) ',&
     &   'OF ENERGY IN SOLID GRAINS'/11X,1PD15.7,5X,'NET RATE OF ',     &
     &   'INCREASE(+)/DECREASE(-) OF ENERGY IN FLUID DUE TO CHANGE IN ',&
     &   'MASS OF FLUID'//11X,1PD15.7,5X,'NET ZERO-ORDER PRODUCTION(+)/'&
     &  ,'LOSS(-) OF ENERGY IN FLUID'/11X,1PD15.7,5X,'NET ZERO-ORDER ', &
     &   'PRODUCTION(+)/LOSS(-) OF ENERGY IN SOLID GRAINS'              &
     &   /11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF ENERGY THROUGH FLUID ',&
     &   'SOURCES AND SINKS'/11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF ',  &
     &   'ENERGY THROUGH INFLOWS OR OUTFLOWS AT POINTS OF SPECIFIED ',  &
     &   'PRESSURE'/11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF ENERGY ',    &
     &   'THROUGH ENERGY SOURCES AND SINKS'/11X,1PD15.7,5X,'NET GAIN(+)'&
     &  ,'/LOSS(-) OF ENERGY THROUGH POINTS OF SPECIFIED TEMPERATURE')  
!                                                                       
 1645    NSOPI = NSOP - 1 
         IF (NSOPI.EQ.0) GOTO 2000 
         IF (ME) 1649, 1649, 1659 
 1649    WRITE (fLST, 1650) 
 1650 FORMAT(///22X,'SOLUTE SOURCES OR SINKS AT FLUID SOURCES AND ',    &
     &   'SINKS'//22X,' NODE',8X,'SOURCE(+)/SINK(-)'/32X,               &
     &   '(SOLUTE MASS/SECOND)'/)                                       
         GOTO 1680 
 1659    WRITE (fLST, 1660) 
 1660 FORMAT(///22X,'ENERGY SOURCES OR SINKS AT FLUID SOURCES AND ',    &
     &   'SINKS'//22X,' NODE',8X,'SOURCE(+)/SINK(-)'/37X,               &
     &   '(ENERGY/SECOND)'/)                                            
 1680    DO 1900 IQP = 1, NSOPI 
            I = IABS (IQSOP (IQP) ) 
            IF ( QIN(I) ) 1700, 1700, 1750 
 1700       QU = QIN(I) * CW * UVEC(I, K) 
            GOTO 1800 
 1750       QU = QIN(I) * CW * UIN(I, K) 
 1800       WRITE (fLST, 450) I, QU 
 1900    END DO 
!                                                                       
 2000    IF (NPBC.EQ.0) GOTO 4500 
         IF (ME) 2090, 2090, 2150 
 2090    WRITE (fLST, 2100) 
 2100 FORMAT(///22X,'SOLUTE SOURCES OR SINKS DUE TO FLUID INFLOWS OR ', &
     &   'OUTFLOWS AT POINTS OF SPECIFIED PRESSURE'//22X,' NODE',8X,    &
     &   'SOURCE(+)/SINK(-)'/32X,'(SOLUTE MASS/SECOND)'/)               
         GOTO 2190 
 2150    WRITE (fLST, 2160) 
 2160 FORMAT(///22X,'ENERGY SOURCES OR SINKS DUE TO FLUID INFLOWS OR ', &
     &   'OUTFLOWS AT POINTS OF SPECIFIED PRESSURE'//22X,' NODE',8X,    &
     &   'SOURCE(+)/SINK(-)'/37X,'(ENERGY/SECOND)'/)                    
 2190    DO 2400 IP = 1, NPBC 
            I = IABS ( SpecifiedPBC(IP)%node ) 
            IF (QPLITR (IP) ) 2200, 2200, 2250 
 2200       QPU = QPLITR (IP) * CW * UVEC(I, K) 
            GOTO 2300 
 2250       QPU = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K) 
 2300       WRITE (fLST, 450) I, QPU 
 2400    END DO 
!                                                                       
         IF (IBCT.EQ.4) GOTO 4500 
         NSOUI = NSOU (K) 
         INEGCT = 0 
         DO 3500 IQU = 1, NSOUI 
            I = IQSOU (IQU, K) 
            IF (I) 3400, 3500, 3500 
 3400       INEGCT = INEGCT + 1 
            IF (ME) 3450, 3450, 3460 
 3450       IF (INEGCT.EQ.1) WRITE (fLST, 3455) 
 3455 FORMAT(///22X,'TIME-DEPENDENT SOLUTE SOURCES AND SINKS'//22X,     &
     &   ' NODE',10X,'GAIN(+)/LOSS(-)'/30X,'  (SOLUTE MASS/SECOND)'//)  
            GOTO 3475 
 3460       IF (INEGCT.EQ.1) WRITE (fLST, 3465) 
 3465 FORMAT(///22X,'TIME-DEPENDENT ENERGY SOURCES AND SINKS'//22X,     &
     &   ' NODE',10X,'GAIN(+)/LOSS(-)'/35X,'  (ENERGY/SECOND)'//)       
 3475       CONTINUE 
            WRITE (fLST, 3490) - I, QUIN ( - I, K) 
 3490 FORMAT(22X,I9,10X,1PD15.7) 
 3500    END DO 
!                                                                       
 4500    IF (NUBC (K) .EQ.0) GOTO 5000 
         IF (K.EQ.NESP) GOTO 4610 
 4600    WRITE (fLST, 4650) 
 4650 FORMAT(///22X,'SOLUTE SOURCES OR SINKS DUE TO SPECIFIED ',        &
     &   'CONCENTRATIONS'//22X,' NODE',10X,'GAIN(+)/LOSS(-)'/30X,       &
     &   '  (SOLUTE MASS/SECOND)'/)                                     
         GOTO 4690 
 4610    WRITE (fLST, 4660) 
 4660 FORMAT(///22X,'ENERGY SOURCES OR SINKS DUE TO SPECIFIED ',        &
     &   'TEMPERATURES'//22X,' NODE',10X,'GAIN(+)/LOSS(-)'/35X,         &
     &   '  (ENERGY/SECOND)'/)                                          
 4690    CONTINUE 
         DO 4700 IU = 1, NUBC (K) 
            IUP = IU + NPBC 
            I = IABS ( MultiSpeciesBC(K)%SpecifiedU(IU)%node ) 
            WRITE (fLST, 450) I, QPLITR (IUP) 
 4700    END DO 
!                                                                       
 5000 END DO 
!.....END OF SPECIES LOOP                                               
!                                                                       
 5500 CONTINUE 
!                                                                       
!.....RESET VALUES FOR ENERGY TRANSPORT FOR SUBSEQUENT CALCULATIONS     
      CS = CST 
!      CW = CWT 
!                                                                             
      RETURN 
      END SUBROUTINE BUDGET                         
