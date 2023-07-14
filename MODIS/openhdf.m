pat='c:/documents and settings/login/my documents/satellite_dat/MODIS';
lis=dir(pat);

nfile=3;


pat=[pat '/' lis(nfile).name];
sd_id=hdfsd('start',pat,'read');


%stat = hdfsd('end',sd_id); %command to close