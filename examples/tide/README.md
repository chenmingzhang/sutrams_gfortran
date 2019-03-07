EXAMPLE

1. RUN THE MODEL IN FOLDER “INPUT” 
2. MOVE THESE FILES TO THE SAME FOLDER OF THE MATLAB POST-PROCESSING CODES
    a. 'M5_25.NOD'
    b. 'M5_25.ELE'
    c. 'NODEWISE_FLUID_MSSS.DAT'
    d. 'NODEWISE_SOLUTE_MASS.DAT'
3. PREPARE THE DATA
   Run this in matlab command: 
LABTIDE_average(114,0,421)
   You should have 4 files generated in the same folder, including:
    a. 'NODDATA.mat'
    b. 'ELEDATA.mat'
    c. 'FLUIDMASS.mat'
    d. 'SALTMASS.mat'   
4. PLOT CONTOUR
    figure;
    subplot(3,1,1); LABTIDE_contour(4,115); text(0,0.95,'a) Pressure');
    subplot(3,1,2); LABTIDE_contour(5,115); text(0,0.95,'b) Species #1');
    subplot(3,1,3); LABTIDE_contour(6,115); text(0,0.95,'c) Species #2');   
5. PLOT VELOCITY VECTORS
    figure;
    LABTIDE_contour(6,115); hold on
    LABTIDE_quiver(115,5);  
6. PLOT STREAMLINE
    Figure;
    LABTIDE_contour(6,115); hold on
    LABTIDE_streamline(115,1.85,0.5);   
    LABTIDE_streamline(115,1.955,0.4);
    LABTIDE_streamline(115,1.955,0.2);
    LABTIDE_streamline(115,1.4,0.72);
    LABTIDE_streamline(115,0.2,0.2);    

