function SLPTIDE_streamline(ts,startx,startz)
load('ELEDATA.mat');
xa1=linspace(1.16,1.96,1001);xa2=linspace(0,1.159,500);xav=[xa2 xa1];
zav=linspace(0,0.872,500);
[xav1,zav1]=meshgrid(xav,zav);
for i=1:500
    for j=1:1501
        zav1(i,j)=0.872*(i-1)/499;
    end
end
    vx=griddata(xcoor,zcoor,vx_ele(:,ts),xav1,zav1,'linear');
    vz=griddata(xcoor,zcoor,vz_ele(:,ts),xav1,zav1,'linear');
%
stream_ = stream2(xav1,zav1,vx,vz,startx,startz); 
stream = stream_;
for j = 1:100  % Matlab stops every 10000 records, so the loop is used to continue the calculation 
    TF = isnan(stream{1,1}(length(stream{1,1}),1));
    if TF == 0
        startx = stream{1,1}(length(stream{1,1}),1);startz = stream{1,1}(length(stream{1,1}),2);
        stream_ = stream2(xav1,zav1,vx,vz,startx,startz);            
        stream{1,1} = [stream{1,1};stream_{1,1}];
    else
        break
    end        
end
stream{1,1}(any(isnan(stream{1,1}), 2), :) = [];
plot(stream{1,1}(:,1),stream{1,1}(:,2),'w','LineWidth',1.5);
xlabel('{\itx} (m)'); ylabel('{\itz} (m)'); grid on
save('streamline.mat','stream');