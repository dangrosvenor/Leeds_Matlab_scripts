%routines for manipulating Johnanne's data

%run JohannesLoad first

%c contains ice no.conc in #/kg * 1e-6 (per cc)
%dimensions 703 1 123 83
%rho contains air density kg/m^3

c2=c.*rho; %number conc per cm^3


maxNCmpc=max(c2,1);
maxNCmpc=squeeze(max(maxNCmpc,[],3)); %timeseries of max number conc

meanNCmpc=mean(c2,1);
meanNCmpc=squeeze(mean(meanNCmpc,3)); %timeseries of mean number conc

for i=1:size(c2,4)
%[i1 i2 i3 i4]=ind2sub(size(c2),find(c2>0)); %ind2sub converts linear indices into multidimensional ones
ii=find(c2(:,:,:,i)>0);
temp=c2(:,:,:,i);
meanNImpc2(i)=mean(temp(ii));

end

%load sdla_xxx.mat where xxx=e.g. icenc for various HM fields
%load sdla_Potemp_Johannes.mat & sdla_Pressure_Johannes.mat
rhoLEM=pressure(1).p.^(1-0.286).*28.97e-3.*1e5^0.286 ./8.31 ./potemp(1).p;
maxNIlem=max(icenc(1).i.*rhoLEM,[],1);
maxNIlem=squeeze(max(maxNIlem,[],2))*1e-6;

maxNSlem=max(snownc(1).i.*rhoLEM,[],1);
maxNSlem=squeeze(max(maxNSlem,[],2))*1e-6;

maxNGlem=max(graupelnc(1).i.*rhoLEM,[],1);
maxNGlem=squeeze(max(maxNGlem,[],2))*1e-6;


meanNIlem=mean(icenc(1).i.*rhoLEM,1);
meanNIlem=squeeze(mean(meanNIlem,2))*1e-6;

meanNSlem=mean(snownc(1).i.*rhoLEM,1);
meanNSlem=squeeze(mean(meanNSlem,2))*1e-6;

meanNGlem=mean(graupelnc(1).i.*rhoLEM,1);
meanNGlem=squeeze(mean(maxNGlem,2))*1e-6;

for i=1:size(icenc(1).i,3)
    %[i1 i2 i3]=ind2sub(size(icenc(1).i),find(icenc(1).i>1)); %ind2sub converts linear indices into multidimensional ones
    temp=icenc(1).i(:,:,i);
    ii=find(temp>1);
    temp(temp==1)=0;
    temp=temp.*rhoLEM(:,:,i);
    
    temp2=snownc(1).i(:,:,i);
    ii2=find(temp2>1);
    temp2(temp2==1)=0;
    temp2=temp2.*rhoLEM(:,:,i);
        
    temp3=graupelnc(1).i(:,:,i);
    ii3=find(temp3>1);
    temp3(temp3==1)=0;
    temp3=temp3.*rhoLEM(:,:,i);

    iii=unique([ii;ii2;ii3]);
    meanNIlem2(i)=squeeze(mean(temp(iii)+temp2(iii)+temp3(iii)))*1e-6;
    
    
end

for i=1:size(icenc(1).i,3)
    %[i1 i2 i3]=ind2sub(size(icenc(1).i),find(icenc(1).i>1)); %ind2sub converts linear indices into multidimensional ones
    temp=icenc(1).i(:,:,i);
    ii=find(temp>1);
    temp(temp==1)=0;
    temp=temp.*rhoLEM(:,:,i);
	meanNIlem2(i)=squeeze(mean(temp(ii)))*1e-6;
end