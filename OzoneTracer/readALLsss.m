%read in all dmi sonde data and put into data(i).dmi



Mx=48; %molecular weight of ozone
Ma=28.97; %of air

pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/O3SSS';
%pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/MICROLIDAR/flight_SF4/040224_MULID_AMES.txt';

firstvals={'Pot.temp./K'};
nocols=8;

lis=dir(pat);

for i=3:length(lis)
    pat2=strcat(pat,'/',lis(i).name);
    lis2=dir(pat2);
    pat2=strcat(pat2,'/',lis2(3).name);
       
    data(i-2).sss=readtable(pat2,firstvals,nocols);
    
    data(i-2).sss(9,:)=data(i-2).sss(2,:)*1e-3*Mx./(data(i-2).sss(3,:)*100*Ma); %conversion to mixing ratio from partial pressure
    data(i-2).sss(10,:)=data(i-2).sss(9,:)*1e6*Ma./Mx; %conversion to ppmv

    
    
end
 
 
 

 
 
 
%  nans=find(diracd==999);
%  diracd(nans)=NaN;
%  nans=find(diracd==99999);
%  diracd(nans)=NaN;

