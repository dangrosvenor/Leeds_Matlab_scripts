time=1; %chosen 07 UTC
[ilat,ilon]=getind_latlon_quick(lat2d.var,lon2d.var,24.08,38.06,0.1);
qv = nc{'QVAPOR'}(time,:,ilat,ilon);
tc = WRFUserARW(nc,'tc',time,ilat,ilon);
p = nc{'P'}(time,:,ilat,ilon) + nc{'PB'}(time,:,ilat,ilon);
z = WRFUserARW(nc,'Z',time,ilat,ilon);
u2d = 0.5* (nc{'U'}(time,:,:,1:end-1) + nc{'U'}(time,:,:,2:end) );
v2d = 0.5* (nc{'V'}(time,:,1:end-1,:) + nc{'V'}(time,:,2:end,:) ); 
elevation = nc{'HGT'}(time,ilat,ilon);

u = u2d(:,ilat,ilon);
v = v2d(:,ilat,ilon);

tc=tc+273.15;
tc = tc * 1.02; %make temperature 5 % higher


f=1e6*28.97/18;
qs = SatVapPress(tc,'goff','liq',p,1) / f; % divide by f to give qs in kg/kg

rh = 100 * qv / qs;

filename = '/home/mbexddg5/OBS_DOMAIN201';

fid = fopen(filename,'wt');
fprintf(fid,' %s%s%s%s0000\n',Times(time,1:4),Times(time,6:7),Times(time,9:10),Times(time,12:13));
fprintf(fid,'  %7.2f   %.2f   \n',lat2d.var(ilat,ilon),lon2d.var(ilat,ilon));


numz = length(z);
id='ID'; id(end+1)='_'; id(end+1:40)='X';
namef='NAMEF'; namef(end+1)='_'; namef(end+1:40)='X';
fprintf(fid,'  %s   %s   \n',id,namef);

platform='PLATFORM'; platform(end+1)='_'; platform(end+1:16)='X';
source='SOURCE'; source(end+1)='_'; source(end+1:16)='X';
fprintf(fid,'  %s  %s      %.0f.     T     F     %d\n',platform,source,elevation,numz);

rh = rh*1.2;
rh(rh>100)=100;

qc=10;

for i = 1:length(z)

	fprintf(fid,' %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f \n',p(i),qc,z(i),qc,tc(i),qc,u(i),qc,v(i),qc,rh(i),qc);

end

fclose(fid);



