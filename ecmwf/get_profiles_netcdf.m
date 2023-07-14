comp='pc';
%comp='laptop';

ecmwf_case='24thFeb_profiles';
ecmwf_case='24thFeb_surfaces2';
ecmwf_case='MilesCity_profiles';
%   inqnc(file) gives you the model fields available (type -1 to exit from this)

switch comp
case 'laptop'
    direcRoot='c:/documents and settings/G/my documents/';
case 'pc'
    direcRoot='c:/documents and settings/login/my documents/';
end

direc=[direcRoot 'ecmwf/'];

exdirA=[direc 'model_levels.txt'];
dat=dlmread(exdirA,' ');
ecmwf(1).p=dat(2:61,4);

switch ecmwf_case
case 'MilesCity_profiles'
    direc='c:/documents and settings/login/my documents/LEEDS_MMOCCA/ecmwf/';
    
	exdirA=[direc 'ecm-e40_1deg-model-levs_1981071918-1981072000_q.nc'];
	
	ecmwf(1).q=getnc(exdirA,'q');
	ecmwf(1).lon=getnc(exdirA,'longitude');
	ecmwf(1).lat=getnc(exdirA,'latitude');
	ecmwf(1).zi=getnc(exdirA,'level');
    ecmwf(1).time=getnc(exdirA,'time');
	
	
	exdirA=[direc 'ecm-e40_1deg-model-levs_1981071918-1981072000_t.nc'];
	ecmwf(1).t=getnc(exdirA,'t');
	
	exdirA=[direc 'ecm-e40_1deg-model-levs_1981071918-1981072000_u.nc'];
	ecmwf(1).u=getnc(exdirA,'u');
	
	exdirA=[direc 'ecm-e40_1deg-model-levs_1981071918-1981072000_v.nc'];
	ecmwf(1).v=getnc(exdirA,'v');
    
case '24thFeb_profiles'
	exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_q.nc'];
	
	ecmwf(1).q=getnc(exdirA,'q');
	ecmwf(1).lon=getnc(exdirA,'longitude');
	ecmwf(1).lat=getnc(exdirA,'latitude');
	ecmwf(1).zi=getnc(exdirA,'level');
    ecmwf(1).time=getnc(exdirA,'time');
	
	
	exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_t.nc'];
	ecmwf(1).t=getnc(exdirA,'t');
	
	exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_u.nc'];
	ecmwf(1).u=getnc(exdirA,'u');
	
	exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_v.nc'];
	ecmwf(1).v=getnc(exdirA,'v');
	
	
    
case '24thFeb_surfaces'
    exdirA=[direc 'near_surf/ecm-op_1deg-model-levs_2004022500_q.nc'];	
	ecmwf(1).q=getnc(exdirA,'q');
	ecmwf(1).lon=getnc(exdirA,'longitude');
	ecmwf(1).lat=getnc(exdirA,'latitude');
	ecmwf(1).zi=getnc(exdirA,'level');
	
	
	exdirA=[direc 'near_surf/ecm-op_1deg-model-levs_2004022500_t.nc'];
	ecmwf(1).t=getnc(exdirA,'t');
	
    

case '24thFeb_surfaces2'
    exdirA=[direc 'near_surf/ecm-op_1deg-model-levs_2004022418-2004022500_q.nc'];	
	ecmwf(1).q=getnc(exdirA,'q');
	ecmwf(1).lon=getnc(exdirA,'longitude');
	ecmwf(1).lat=getnc(exdirA,'latitude');
	ecmwf(1).zi=getnc(exdirA,'level');
	
	
	exdirA=[direc 'near_surf/ecm-op_1deg-model-levs_2004022418-2004022500_t.nc'];
	ecmwf(1).t=getnc(exdirA,'t');
    
    exdirA=[direc 'near_surf/ecm-op_1deg-model-levs_2004022418-2004022500_u.nc'];
	ecmwf(1).u=getnc(exdirA,'u');
    
    exdirA=[direc 'near_surf/ecm-op_1deg-model-levs_2004022418-2004022500_v.nc'];
	ecmwf(1).v=getnc(exdirA,'v');
	
    
end



parr=flipud(ecmwf(1).p);
if length(size(ecmwf(1).t))==4
    lenz=size(ecmwf(1).t,2);
    z = hyd_pressure_height( parr(1:lenz ) , ecmwf(1).t(1,:,1,1),0);
else
    lenz=size(ecmwf(1).t,1);
    z = hyd_pressure_height( parr(1:lenz ) , ecmwf(1).t(:,1,1),0);

end

ecmwf(1).z=z';

ecmwf(1).p2=flipud(ecmwf(1).p)*100; %flip so that is in correct orientation (pressure decreasing) and is in Pa

% run mean_vap_BL_ecmwf for mean boundary layer calcs
lon=ecmwf(1).lon-360;
lat=ecmwf(1).lat;
   
'done read ecmwf'