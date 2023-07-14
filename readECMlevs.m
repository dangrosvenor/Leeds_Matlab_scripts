 pat='c:/pcgrads/win32e/ecmwf/60 model levels3.txt';
 fid=fopen(pat,'rt');
 
 temp=fscanf(fid,'%g',[1]);
 while length(temp)==0
     fscanf(fid,'%s',[1]);
     temp=fscanf(fid,'%g',[1]);
 end
 
 lev(1,1)=temp;
 lev(2:4,1)=fscanf(fid,'%g',[3]);
 

 lev(1:5,2:61)=fscanf(fid,'%g',[5 60]);
     
 fclose(fid);
 
 
 a=lev(2,:);
 b=lev(3,:);
 psurf=900;
 
 p=a+b.*psurf;