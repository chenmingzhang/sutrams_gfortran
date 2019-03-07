% Calculate average size of the USP in tidal conditions
function SLPTIDE_average(ts,nfwbc,nswbc)
    
    nn = 76466; % Number of nodes
    ne = 75900; % Number of elements
	
	% read NODE and ELEMENT data from Sutra MS output
    nod = readNOD_MS('M5_25','outputnumber',ts); % readNOD_MS('filename','outputnumber',number of timestep)
    ele = readELE_MS('M5_25','outputnumber',ts); % readELE_MS('filename','outputnumber',number of timestep)
    
    % Parse salinity/temperature/pressure data of each time step
    for i = 1:ts
        nod_data = nod(i);
        p_nod(:,i) = nod_data.terms{4}; % pressure data
        t_nod(:,i) = nod_data.terms{5}; % data of species #1
        s_nod(:,i) = nod_data.terms{6}; % data of species #2
    end
        x_nod = nod_data.terms{2};		% x coordinate
        z_nod = nod_data.terms{3};		% z coordinate
	
	% Calculate average over one tidal cycle
    for i = 1:nn
        p_nod(i,ts+1) = mean(p_nod(i,1:ts)); 
        t_nod(i,ts+1) = mean(t_nod(i,1:ts)); 
        s_nod(i,ts+1) = mean(s_nod(i,1:ts)); 
    end    
    % save
    save ('NODDATA.mat','x_nod','z_nod','p_nod','t_nod','s_nod','nn');

    % Parse velocity data of each time step
    for i = 1:ts
        ele_data = ele(i);
        vx_ele(:,i) = ele_data.terms{4};
        vz_ele(:,i) = ele_data.terms{5};
    end    
        xcoor(:,1) = ele_data.terms{2};
        zcoor(:,1) = ele_data.terms{3};
		
	% Calculate average over one tidal cycle
    for i = 1:ne
        vx_ele(i,ts+1) = mean(vx_ele(i,1:ts)); 
        vz_ele(i,ts+1) = mean(vz_ele(i,1:ts)); 
    end        
    % save
    save ('ELEDATA.mat','xcoor','zcoor','vx_ele','vz_ele','ne');
    
    % FLUIDMASS (UNIT: M3/M/DAY)
    fluidmass = load('NODEWISE_FLUID_MSSS.DAT');
    flmsfreshw = zeros(nfwbc,ts+1);
	flmssaltw = zeros(nswbc,ts+1);
    for i = 1:ts
        flmsfreshw(1:nfwbc,i) = fluidmass((i-1)*(nfwbc+nswbc)+1:(i-1)*(nfwbc+nswbc)+nfwbc,2)*864;
        flmssaltw(1:nswbc,i) = fluidmass((i-1)*(nfwbc+nswbc)+nfwbc+1:(i-1)*(nfwbc+nswbc)+nfwbc+nswbc,2)*864;
    end    
    for i = 1:nfwbc
        flmsfreshw(i,ts+1) = mean(flmsfreshw(i,1:ts)); 
    end    
    for i = 1:nswbc
        flmssaltw(i,ts+1) = mean(flmssaltw(i,1:ts)); 
    end
	% save
    save('FLUIDMASS.mat','flmsfreshw','flmssaltw');
    
    % SALT MASS (UNIT: KG/M/DAY)
    solutemass = load('NODEWISE_SOLUTE_MASS.DAT');
	saltmass = zeros(nswbc,ts+1);
    for i = 1:ts
        saltmass(1:nswbc,i) = solutemass((2*i-1)*(nfwbc+nswbc)+nfwbc+1:2*i*(nfwbc+nswbc),3)*864000; 
    end
    for i = 1:nswbc
        saltmass(i,ts+1) = mean(saltmass(i,1:ts)); 
    end
	% save
    save('SALTMASS.mat','saltmass');
    
    
    
    
    
    
    
    
    
    
    
    