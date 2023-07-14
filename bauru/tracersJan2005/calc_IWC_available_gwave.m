clear meanW meanT maxSatIce maxSatTot

if ~exist('TT')
    TT=potemp(1).p./(1e5./pressure(1).p).^0.286;   
end

[xi xi2]=findheight(x,200e3,500e3);

[zi zi2]=findheight(z,15e3,22e3);

%temp=windVW(1).W(zi:zi2,xi:xi2,:);
temp2=TT(zi:zi2,xi:xi2,:);
temp3=vap(1).v(zi:zi2,xi:xi2,1:83);
for i=1:size(temp2,1)
    for j=1:size(temp2,3)
        %ii=find(temp(i,:,j)>0);
        %meanW(i,j)=mean(temp(i,ii,j));
        %maxW(i,j)=mean(temp(i,ii,j));
        
        %meanVapW(i,j)=mean(temp3(i,ii,j));
        %maxVapW(i,j)=max(temp3(i,ii,j));
        
        P=pressure(1).p(i+zi-1,xi:xi2,j);
        satmr=(SatVapPress(temp2(i,:,j),'teten','ice',P)./P * 18/28.79);
        sat=temp3(i,:,j)./satmr;
        totice=icemr(1).i(i+zi-1,xi:xi2,j) + snowmr(1).i(i+zi-1,xi:xi2,j) + graupelmr(1).i(i+zi-1,xi:xi2,j);
        satTotW= ( temp3(i,:,j)+totice/1e3 )./satmr; %(IWC + vapour mr)/icesatmr
       
% [v k]=max(temp3(i,:,j));
 %       P=pressure(1).p(i+zi-1,k+xi-1,j);
  %      t=TT(i+zi-1,k+xi-1,j);
        %maxSatIce(i,j)=v./(SatVapPress(t,'teten','ice',P)/P * 18/28.79);
        maxSatIce(i,j)=max(sat);
        maxSatTot(i,j)=max(satTotW);
        
        meanT(i,j)=mean(temp2(i,:,j))-273.15;
        minT(i,j)=min(temp2(i,:,j))-273.15;
    end
end

