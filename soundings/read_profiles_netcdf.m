comp='pc';
%comp='laptop';


switch comp
case 'laptop'
    direc='c:/documents and settings/G/my documents/ecmwf/';
case 'pc'
    direc='c:/documents and settings/login/my documents/ecmwf/';
end

exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_q.nc'];

ecmwf(1).q=getnc(exdirA,'q');
ecmwf(1).lon=getnc(exdirA,'longitude');
ecmwf(1).lat=getnc(exdirA,'latitude');
ecmwf(1).zi=getnc(exdirA,'level');


exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_t.nc'];
ecmwf(1).t=getnc(exdirA,'t');

exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_u.nc'];
ecmwf(1).u=getnc(exdirA,'u');

exdirA=[direc 'ecm-op_1deg-model-levs_2004022418_v.nc'];
ecmwf(1).v=getnc(exdirA,'v');

exdirA=[direc 'model_levels.txt'];
dat=dlmread(exdirA,' ');
ecmwf(1).p=dat(2:61,4);

'done read ecmwf'