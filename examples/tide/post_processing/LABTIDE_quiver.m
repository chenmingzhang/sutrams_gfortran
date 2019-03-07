function LABSLP_quiver(ts,scale)
load('ELEDATA.mat');
xa1=linspace(1.16,1.96,1001);xa2=linspace(0,1.159,1000);xav=[xa2 xa1];
zav=linspace(0,0.872,500);
[xav1,zav1]=meshgrid(xav,zav);
for i=1:500
    for j=1:1000
        zav1(i,j)=0.872*(i-1)/199;
    end
    for j=1001:2001 
        zav1(i,j)=(0.872-(0.872-0.472)*(j-1001)/1000)*(i-1)/199;
    end 
end
for i = 1:ne
    vxele(i,1) = vx_ele(i,t); 
    vzele(i,1) = vz_ele(i,t); 
end
    vx=griddata(xcoor,zcoor,vxele,xav1,zav1,'linear');
    vz=griddata(xcoor,zcoor,vzele,xav1,zav1,'linear');

quiver(xav1(1:10:end,1:25:end),zav1(1:10:end,1:25:end),vx(1:10:end,1:25:end),vz(1:10:end,1:25:end),scale,'Color',[0.65 0.65 0.65]);
