%creates flux and rougness lenght namelists
twod=1; %flag for 2d
ismooth=1;

%direc(1).d='c:\program files\pcgrads\win32e\07.02-00\vslicesen1';
%direc(2).d='c:\program files\pcgrads\win32e\07.02-00\vslicelat1';
%direc(3).d='c:\program files\pcgrads\win32e\07.02-00\vsliceZ0';

oridefile=0;
readGRADS;


if twod==1
    nps=200;
%     lata=-19.015;
%     lona=-52.582;
%     latb=-25.4184;
%     lonb=-45.7704;

    lata=-20.4523;
    lona=-51.1057;
    latb=-23.6318;
    lonb=-47.6462;


%     latb=-25.367;
%     lonb=-45.8218;
    
    slat=lata+((latb-lata)/nps)*(0:nps);
    slon=lona+((lonb-lona)/nps)*(0:nps);
    %slon=lon2(1)+(0:nps)*(lon2(end)-lon2(1))/nps;
    %slat=lat(1)+(nps:-1:0)*(lat(end)-lat(1))/nps;
    ziCBNT=interp2(lon,lat,pr(4).p,slon,slat);
    ziSEN=interp2(lon,lat,pr(1).p,slon,slat);
    ziLAT=interp2(lon,lat,pr(2).p,slon,slat);
    ziZ0=interp2(lon,lat,pr(3).p,slon,slat);
    
    
end



outfile='c:\cygwin\home\user\lesexe\les2.3\fluxvardata1000km-500pts.f';
outfile2='c:\cygwin\home\user\lesexe\les2.3\roughvardata1000km-500pts.f';

clear gradsen gradlat z0 z0th

fid=fopen(outfile,'w');
fid2=fopen(outfile2,'w');


np=500;  %2187;
n=201; %size input array

m=(np-1)/(n-1);
mm=round(np/10);

mm1=round(1*np/n);
mm1=1;
mm2=round(185*np/n);
mm2=485;



for i=1:n-1
    ns(i)=1+round((i-1)*m);
end

ns(n)=np;

fprintf(fid,'&FLUXVAR\n',np);
fprintf(fid2,'&ROUGHVAR\n',np);

%fprintf(fid,'      REAL SFLUX(%d)\n',np);
%fprintf(fid,'      REAL LFLUX(%d)\n',np);



for i=1:np
    gradsen(i)=interp1(ns,ziSEN',i,'linear');
    gradlat(i)=interp1(ns,ziLAT',i,'linear');
    z0(i)=interp1(ns,ziZ0',i,'linear');
    z0th(i)=z0(i)/10;
end

    

if ismooth==1
	for i=mm2:np    %for smoothing end values toward start ones - no need as discontinuity exists in reality
        gradsen(i)=interp1([np-mm2 np+mm1],[gradsen(np-mm2) gradsen(mm1)],i,'linear');
        gradlat(i)=interp1([np-mm2 np+mm1],[gradlat(np-mm2) gradlat(mm1)],i,'linear');
        z0(i)=interp1([np-mm2 np+mm1],[z0(np-mm2) z0(mm1)],i,'linear');
        z0th(i)=z0(i)/10;
	end
	
	for i=1:mm1-1    %for smoothing end values toward start ones - no need as discontinuity exists in reality
        gradsen(i)=interp1([0 mm1-1],[gradsen(np) gradsen(mm1)],i,'linear');
        gradlat(i)=interp1([0 mm1-1],[gradlat(np) gradlat(mm1)],i,'linear');
        z0(i)=interp1([0 mm1-1],[z0(np) z0(mm1)],i,'linear');
        z0th(i)=z0(i)/10;
	end

end

meansen=-mean(gradsen);
meanlat=-mean(gradlat);



fprintf(fid,'\n');
for i=1:np
    fprintf(fid,'SFLUX(%d)=%.2f,\n',i,(-gradsen(i)/meansen));
end

for i=1:np
    fprintf(fid,'LFLUX(%d)=%.2f,\n',i,(-gradlat(i)/meanlat));
end

for i=1:np
    fprintf(fid2,'Z0var(%d)=%.2f,\n',i,z0(i));
end

for i=1:np
    fprintf(fid2,'Z0THvar(%d)=%.3f,\n',i,z0th(i));
end

fprintf(fid,'&END',np);
fprintf(fid2,'&END',np);
    
fclose(fid);
fclose(fid2);