pat='c:/documents and settings/login/my documents/hibiscus/baloonDATA/DIRAC/flight_SF4/20040224.b20';
 fid=fopen(pat,'rt');
 
 %temp=fscanf(fid,'%g',[1]);
 cont=1;
 while (cont==1)
     fscanf(fid,'%s',[1]);
     temp=fscanf(fid,'%g',[1]);
     if temp==67680;cont=0;end
 end

 diracd(1,1)=temp;
 diracd(2:9,1)=fscanf(fid,'%g',[8]);
 

 dirac2=fscanf(fid,'%g',[9 Inf]);
 diracd(:,2:size(dirac2,2)+1)=dirac2;    
 fclose(fid);
 
 nans=find(diracd==999);
 diracd(nans)=NaN;
 nans=find(diracd==99999);
 diracd(nans)=NaN;