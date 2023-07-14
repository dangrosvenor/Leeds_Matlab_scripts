direc(1).d='c:\program files\pcgrads\win32e\07.02-00\vslicesen1';
oridefile=1;
readGRADS;
pr(2).p=pr(1).p;

direc(1).d='c:\program files\pcgrads\win32e\07.02-00\vslicelat1';
oridefile=1;
readGRADS;


direc(1).d='c:\cygwin\home\user\lesexe\les2.3\fluxdata108.f';

clear gradsen gradlat

fid=fopen(direc(1).d,'w');
np=108;
n=26;

m=(np-1)/(n-1);
mm=round(np/5);

meansen=-mean(pr(2).p);
meanlat=-mean(pr(1).p);

for i=1:n-1
    ns(i)=1+round((i-1)*m);
end

ns(n)=np;


fprintf(fid,'      REAL SFLUX(%d)\n',np);
fprintf(fid,'      REAL LFLUX(%d)\n',np);

for i=1:np
    gradsen(i)=interp1(ns,pr(2).p,i,'linear');
    gradlat(i)=interp1(ns,pr(1).p,i,'linear');
end

for i=np-mm:np
    gradsen(i)=interp1([np-mm np],[gradsen(np-mm) gradsen(1)],i,'linear');
    gradlat(i)=interp1([np-mm np],[gradlat(np-mm) gradlat(1)],i,'linear');
end

fprintf(fid,'\n');
for i=1:np
    fprintf(fid,'      SFLUX(%d)=%.2f\n',i,(-gradsen(i)/meansen));
end

for i=1:np
    fprintf(fid,'      LFLUX(%d)=%.2f\n',i,(-gradlat(i)/meanlat));
end


    
fclose(fid);