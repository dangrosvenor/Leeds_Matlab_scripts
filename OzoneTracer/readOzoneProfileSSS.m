comp='pc';
%comp='laptop';

switch comp
case 'pc'
	pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/O3SSS/flight_SF4/SS040224.b20';
case 'laptop'
	pat='c:/documents and settings/g/my documents/baloonDATA/O3SSS/flight_SF4/SS040224.b20';
end

%pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/MICROLIDAR/flight_SF4/040224_MULID_AMES.txt';
 fid=fopen(pat,'rt');

 clear diracd
 firstvals=[-1 2.601]; %first values of data table SF4
 %firstvals=[-1 2.469]; %SF1
 %firstvals=[-1 1.676]; %SF2
 %firstvals=[72537 2651]; %SF3
 
 nocols=8; %number of columns in data table =8 for ozone
 
 Mx=48; %molecular weight of ozone
 Ma=28.97; %of air
 %temp=fscanf(fid,'%g',[1]);
 cont=1;
 while (cont==1)
     fscanf(fid,'%s',[1]);
     for i=1:length(firstvals)
        temp=fscanf(fid,'%g',[1]);
        if temp==firstvals(i);
            diracd(i,1)=temp;
            cont=0;
        else
            cont=1;
        end
    end
 end


 if nocols>1
     diracd(i+1:nocols,1)=fscanf(fid,'%g',[nocols-i]);
 end

 dirac2=fscanf(fid,'%g',[nocols Inf]);
 diracd(:,2:size(dirac2,2)+1)=dirac2;    
 fclose(fid);
 
 
 dirac2(9,:)=dirac2(2,:)*1e-3*Mx./(dirac2(3,:)*100*Ma); %conversion to mixing ratio from partial pressure
 dirac2(10,:)=dirac2(9,:)*1e6*Ma./Mx; %conversion to ppmv
%  nans=find(diracd==999);
%  diracd(nans)=NaN;
%  nans=find(diracd==99999);
%  diracd(nans)=NaN;

'done read SSS'

