function SLPTIDE_contour(param,ts)
load('NODDATA.mat')
ctr_meshx = 0.004; ctr_meshz = 0.004;
width = 1.96; height = 0.872; x1 = 1.16;
x2 = 1.96; z1 = 0.472;
meshnodx = width / ctr_meshx;
meshnodz = height / ctr_meshz;

xmatrix=linspace(0,width,meshnodx);
zmatrix=linspace(0,height,meshnodz);
[xnods,znods]=meshgrid(xmatrix,zmatrix);
for i=1:meshnodz
	for j=(x1/ctr_meshx):(x2/ctr_meshx)                                                     
		j1 = x1/ctr_meshx;
		j2 = x2/ctr_meshx - x1/ctr_meshx;
		znods(i,j)=(height-(height - z1)*(j-j1)/j2)*(i-1)/(meshnodz - 1);    
	end                                                                  
end
for i=1:meshnodz
	for j=(x2/ctr_meshx + 1):meshnodx                                                     
		znods(i,j)=(i-1)*z1/meshnodz;
	end                                                                  
end
tpgrid=griddata(x_nod,z_nod,t_nod(:,ts),xnods,znods,'linear');   
ppgrid=griddata(x_nod,z_nod,p_nod(:,ts),xnods,znods,'linear');   
spgrid=griddata(x_nod,z_nod,s_nod(:,ts),xnods,znods,'linear');   
%
if param == 4
    paramplot = ppgrid;
elseif param == 5
    paramplot = tpgrid;
elseif param == 6
    paramplot = spgrid; 
end    
    contourf(xnods,znods,paramplot,100,'linestyle','none'); 
    axis([0 width 0 height]);
    xlabel('{\itx} (m)'); ylabel('{\itz} (m)'); grid on
    colorbar('Color','w','FontSize',11,'Location','West'); colormap('jet')
