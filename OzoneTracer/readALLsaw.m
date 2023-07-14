%read in all dmi sonde data and put into data(i).dmi



Mx=48; %molecular weight of ozone
Ma=28.97; %of air

pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/H2OSAW';
%pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/MICROLIDAR/flight_SF4/040224_MULID_AMES.txt';

firstvals={'UT.', '---'};
nocols=12;

lis=dir(pat);

for i=3:length(lis)
    pat2=strcat(pat,'/',lis(i).name);
    lis2=dir(pat2);
    pat2=strcat(pat2,'/',lis2(3).name);
    
    lis3(i-2).name=lis2(3).name;
    
    if strcmp(lis2(3).name(3:4),'OS')==1;
        nocols=15;
    else
        nocols=12;
    end
    
    data(i-2).sss=readtable(pat2,firstvals,nocols);
    
    %data(i-2).sss(9,:)=data(i-2).sss(2,:)*1e-3*Mx./(data(i-2).sss(3,:)*100*Ma); %conversion to mixing ratio from partial pressure
    %data(i-2).sss(10,:)=data(i-2).sss(9,:)*1e6*Ma./Mx; %conversion to ppmv

    nans=find(data(i-2).sss==-99);
    data(i-2).sss(nans)=NaN;
end
 
'finished so I have'
 

 
 
 
%  nans=find(diracd==999);
%  diracd(nans)=NaN;
%  nans=find(diracd==99999);
%  diracd(nans)=NaN;

