pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/DMIbackscatter/flight_10_24Feb/ba040224.dat';
 fid=fopen(pat,'rt');
 
 %temp=fscanf(fid,'%g',[1]);
 cont=1;
 while (cont==1)
     fscanf(fid,'%s',[1]);
     temp=fscanf(fid,'%g',[1]);
     if temp==0.096;cont=0;end
 end

 diracd(1,1)=temp;
 diracd(2:8,1)=fscanf(fid,'%g',[7]);
 

 dirac2=fscanf(fid,'%g',[8 Inf]);
 diracd(:,2:size(dirac2,2)+1)=dirac2;    
 fclose(fid);
 
%  nans=find(diracd==999);
%  diracd(nans)=NaN;
%  nans=find(diracd==99999);
%  diracd(nans)=NaN;

