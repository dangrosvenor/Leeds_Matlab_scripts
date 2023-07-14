%read in all dmi sonde data and put into data(i).dmi
comp='pc';
%comp='laptop';

switch comp
case 'pc'
	pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/dmibackscatter';
case 'laptop'
	pat='c:/documents and settings/g/my documents/baloonDATA/dmibackscatter';
end


Mx=48; %molecular weight of ozone
Ma=28.97; %of air

%	pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/dmibackscatter';
%pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/MICROLIDAR/flight_SF4/040224_MULID_AMES.txt';

firstvals={'UTC', '---'};
nocols=8;

lis=dir(pat);

for i=3:length(lis)
    pat2=strcat(pat,'/',lis(i).name);
    lis2=dir(pat2);
    pat2=strcat(pat2,'/',lis2(3).name);
       
    data(i-2).dmi=readtable(pat2,firstvals,nocols);
    
    data(i-2).dmi(9,:)=data(i-2).dmi(8,:)*1e-6*48/28.97; %conversion from ppmv to mixing ratio
    
    data(i-2).dmi(1,:)=data(i-2).dmi(1,:)*1e3; %convert from km to m
end
 
 
 'done read DMI'

 
 
 
%  nans=find(diracd==999);
%  diracd(nans)=NaN;
%  nans=find(diracd==99999);
%  diracd(nans)=NaN;

