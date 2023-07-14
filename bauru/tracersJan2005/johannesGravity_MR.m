minice=1e-3; %threshold ice for means in g/kg

[xi xi2]=findheight(x,200e3,500e3);
[xim xim2]=findheight(xmpc,200e3,500e3);

[zi zi2]=findheight(z,15e3,22e3);
[zim zim2]=findheight(zmpc,15e3,22e3);

%for mpc
%maxNCmpc=max(icemrmpc(xim:xim2,1,zim:zim2,:).*rhompc(xim:xim2,1,zim:zim2,:),[],1);
maxMRImpc=max(icemrmpc(xim:xim2,1,zim:zim2,:),[],1)*1e3;
maxMRImpc=squeeze(max(maxNCmpc,[],3)); %timeseries of max number conc

%meanNCmpc=mean(icemrmpc(xim:xim2,1,zim:zim2,:).*rhompc(xim:xim2,1,zim:zim2,:),1);
meanMRImpc=mean(icemrmpc(xim:xim2,1,zim:zim2,:),1)*1e3;
meanMRImpc=squeeze(mean(meanNCmpc,3)); %timeseries of mean number conc

%mean only for points where no. conc>0
for i=1:size(icemrmpc,4)
%[i1 i2 i3 i4]=ind2sub(size(c2),find(c2>0)); %ind2sub converts linear indices into multidimensional ones
%ii=find(icemrmpc(xim:xim2,1,zim:zim2,i).*rhompc(xim:xim2,1,zim:zim2,i)>0);
ii=find(icemrmpc(xim:xim2,1,zim:zim2,i)>0);
%temp=icemrmpc(xim:xim2,1,zim:zim2,i).*rhompc(xim:xim2,1,zim:zim2,i);
temp=icemrmpc(xim:xim2,1,zim:zim2,i)*1e3;
meanMRImpc2(i)=mean(temp(ii));

temp=squeeze(temp);
for ih=1:zim2-zim+1
i2=find(temp(:,ih)>minice);
timHmpcMR(ih,i)=mean(temp(i2,ih),1);
end
end

%load sdla_xxx.mat where xxx=e.g. icenc for various HM fields
%load sdla_Rho_Johannes.mat 
%maxMRIlem=max(icenc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),[],1);
maxMRIlem=max(icemr(1).i(zi:zi2,xi:xi2,:),[],1);
maxMRIlem=squeeze(max(maxMRIlem,[],2));

%maxMRSlem=max(snownc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),[],1);
maxMRSlem=max(snowmr(1).i(zi:zi2,xi:xi2,:),[],1);
maxMRSlem=squeeze(max(maxMRSlem,[],2));

%maxNGlem=max(graupelmr(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),[],1);
maxMRGlem=max(graupelmr(1).i(zi:zi2,xi:xi2,:),[],1);
maxMRGlem=squeeze(max(maxMRGlem,[],2));


%meanMRlem=mean(icemr(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),1);
meanMRIlem=mean(icemr(1).i(zi:zi2,xi:xi2,:),1);
meanMRIlem=squeeze(mean(meanMRIlem,2));

%meanMRSlem=mean(snowmr(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),1);
meanMRSlem=mean(snowmr(1).i(zi:zi2,xi:xi2,:),1);
meanMRSlem=squeeze(mean(meanMRSlem,2));

%meanNGlem=mean(graupelmr(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),1);
meanMRGlem=mean(graupelmr(1).i(zi:zi2,xi:xi2,:),1);
meanMRGlem=squeeze(mean(maxMRGlem,2));

%calculates mean of ice, snow and graupel only for points where mass>minice g/kg

for i=1:size(icemr(1).i,3)
    %[i1 i2 i3]=ind2sub(size(icemr(1).i),find(icemr(1).i>1)); %ind2sub converts linear indices into multidimensional ones
    temp=icemr(1).i(zi:zi2,xi:xi2,i);
    ii=find(temp>minice);
    temp(temp==1)=0;
    %temp=temp.*rhoLEM(zi:zi2,xi:xi2,i);

    
    temp2=snowmr(1).i(zi:zi2,xi:xi2,i);
    ii2=find(temp2>minice);
    temp2(temp2==1)=0;
    %temp2=temp2.*rhoLEM(zi:zi2,xi:xi2,i);


        
    temp3=graupelmr(1).i(zi:zi2,xi:xi2,i);
    ii3=find(temp3>minice);
    temp3(temp3==1)=0;
    %temp3=temp3.*rhoLEM(zi:zi2,xi:xi2,i);


        
    iii=unique([ii;ii2;ii3]);
    meanMRIlem2(i)=squeeze(mean(temp(iii)+temp2(iii)+temp3(iii)));
    
    temp=icemr(1).i(zi:zi2,xi:xi2,i)+snowmr(1).i(zi:zi2,xi:xi2,i)+graupelmr(1).i(zi:zi2,xi:xi2,i);
    for ii=1:zi2-zi+1
        i2=find(temp(ii,:)>minice);
        timHlemMR(ii,i)=mean(temp(ii,i2),2);
    end
    
end
