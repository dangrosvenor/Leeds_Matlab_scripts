[xi xi2]=findheight(x,200e3,500e3);
[xim xim2]=findheight(xmpc,200e3,500e3);

[zi zi2]=findheight(z,16e3,22e3);
[zim zim2]=findheight(zmpc,16e3,22e3);

%for mpc
%maxNCmpc=max(icempc(xim:xim2,1,zim:zim2,:).*rhompc(xim:xim2,1,zim:zim2,:),[],1);
maxNImpc=max(icempc(xim:xim2,1,zim:zim2,:),[],1);
maxNImpc=squeeze(max(maxNCmpc,[],3)); %timeseries of max number conc

%meanNCmpc=mean(icempc(xim:xim2,1,zim:zim2,:).*rhompc(xim:xim2,1,zim:zim2,:),1);
meanNImpc=mean(icempc(xim:xim2,1,zim:zim2,:),1);
meanNImpc=squeeze(mean(meanNCmpc,3)); %timeseries of mean number conc

%mean only for points where no. conc>0
for i=1:size(icempc,4)
%[i1 i2 i3 i4]=ind2sub(size(c2),find(c2>0)); %ind2sub converts linear indices into multidimensional ones
%ii=find(icempc(xim:xim2,1,zim:zim2,i).*rhompc(xim:xim2,1,zim:zim2,i)>0);
ii=find(icempc(xim:xim2,1,zim:zim2,i)>0);
%temp=icempc(xim:xim2,1,zim:zim2,i).*rhompc(xim:xim2,1,zim:zim2,i);
temp=icempc(xim:xim2,1,zim:zim2,i);
meanNImpc2(i)=mean(temp(ii));

end

%load sdla_xxx.mat where xxx=e.g. icenc for various HM fields
%load sdla_Rho_Johannes.mat 
%maxNIlem=max(icenc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),[],1);
maxNIlem=max(icenc(1).i(zi:zi2,xi:xi2,:),[],1);
maxNIlem=squeeze(max(maxNIlem,[],2))*1e-6;

%maxNSlem=max(snownc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),[],1);
maxNSlem=max(snownc(1).i(zi:zi2,xi:xi2,:),[],1);
maxNSlem=squeeze(max(maxNSlem,[],2))*1e-6;

%maxNGlem=max(graupelnc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),[],1);
maxNGlem=max(graupelnc(1).i(zi:zi2,xi:xi2,:),[],1);
maxNGlem=squeeze(max(maxNGlem,[],2))*1e-6;


%meanNIlem=mean(icenc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),1);
meanNIlem=mean(icenc(1).i(zi:zi2,xi:xi2,:),1);
meanNIlem=squeeze(mean(meanNIlem,2))*1e-6;

%meanNSlem=mean(snownc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),1);
meanNSlem=mean(snownc(1).i(zi:zi2,xi:xi2,:),1);
meanNSlem=squeeze(mean(meanNSlem,2))*1e-6;

%meanNGlem=mean(graupelnc(1).i(zi:zi2,xi:xi2,:).*rhoLEM(zi:zi2,xi:xi2,:),1);
meanNGlem=mean(graupelnc(1).i(zi:zi2,xi:xi2,:),1);
meanNGlem=squeeze(mean(maxNGlem,2))*1e-6;

%calculates mean of ice, snow and graupel only for points where no conc>1
for i=1:size(icenc(1).i,3)
    %[i1 i2 i3]=ind2sub(size(icenc(1).i),find(icenc(1).i>1)); %ind2sub converts linear indices into multidimensional ones
    temp=icenc(1).i(zi:zi2,xi:xi2,i);
    ii=find(temp>1);
    temp(temp==1)=0;
    %temp=temp.*rhoLEM(zi:zi2,xi:xi2,i);

    
    temp2=snownc(1).i(zi:zi2,xi:xi2,i);
    ii2=find(temp2>1);
    temp2(temp2==1)=0;
    %temp2=temp2.*rhoLEM(zi:zi2,xi:xi2,i);


        
    temp3=graupelnc(1).i(zi:zi2,xi:xi2,i);
    ii3=find(temp3>1);
    temp3(temp3==1)=0;
    %temp3=temp3.*rhoLEM(zi:zi2,xi:xi2,i);


        
    iii=unique([ii;ii2;ii3]);
    meanNIlem2(i)=squeeze(mean(temp(iii)+temp2(iii)+temp3(iii)))*1e-6;
    
    temp=icenc(1).i(zi:zi2,xi:xi2,i)+snownc(1).i(zi:zi2,xi:xi2,i)+graupelnc(1).i(zi:zi2,xi:xi2,i);
    for ii=1:zi2-zi+1
        i2=find(temp(ii,:)>1);
        timHlemNC(ii,i)=mean(temp(ii,i2),2);
    end
    
end

