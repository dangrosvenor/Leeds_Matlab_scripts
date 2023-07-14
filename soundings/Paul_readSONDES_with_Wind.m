function SONDE=Paul_readSONDES_with_Wind(fileNames)

fid=fopen(fileNames{1},'r');

% Read in data from RS vaisala sondes
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);

str=fgetl(fid); % Start time
ind1=findstr(str,'Started at');
ind2=findstr(str,'UTC');

startTime = datenum(str(ind1+11:ind2-1));
for i=1:34
    fgetl(fid);
end

fclose(fid);

% Now read in data
[min,sec,hPa,gpm,T,RH,Dewp]=textread(fileNames{1},...
    '%d %d %f %f %f %f %f %*[^\n]','headerlines',42);

SONDE.TIME=startTime + datenum(0,0,0,0,min,sec);
SONDE.P=hPa.*100;
SONDE.Z=gpm;
SONDE.T=T+273;
SONDE.RH=RH;
SONDE.TD=Dewp+273;

SONDE.TH = SONDE.T.*(1000./hPa).^0.286;

SONDE.QV = 0.622.*SatVapPress(SONDE.TD,'buck2','liq')./(SONDE.P-SatVapPress(SONDE.TD,'buck2','liq'));
SONDE.QVsat = 0.622.*SatVapPress(SONDE.T,'buck2','liq')./(SONDE.P-SatVapPress(SONDE.T,'buck2','liq'));



% Now read in data
[min,sec,hPa,gpm,m,vel,direc]=textread(fileNames{2},...
    '%d %d %f %f %f %f %f %*[^\n]','headerlines',42);

SONDE.TIME2=startTime + datenum(0,0,0,0,min,sec);
SONDE.P2=hPa.*100;
SONDE.Z21=gpm;
SONDE.Z22=m;
SONDE.vel=vel;
SONDE.direc=direc;
SONDE.U = vel.*sin(direc./180.*pi);
SONDE.V = vel.*cos(direc./180.*pi);


