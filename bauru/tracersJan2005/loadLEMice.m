%exdir2='C:\matlabR12\work\bauru\tracersJan2005\force+3_3th3qv\diags\';
exdir2='G:\runs\dmi1715_5ppmv_25km_3\';

loadcase='pt';

switch loadcase
case 'nc'
	exname=[exdir2 'sdla_icenc_2.mat'];
	load(exname);
	exname=[exdir2 'sdla_snownc.mat'];
	load(exname);
	exname=[exdir2 'sdla_graupelnc.mat'];
	load(exname);
case 'mr'
	exname=[exdir2 'sdla_icemr.mat'];
	load(exname);
	exname=[exdir2 'sdla_snowmr.mat'];
	load(exname);
	exname=[exdir2 'sdla_graupelmr.mat'];
	load(exname);
case 'pt'
	exname=[exdir2 'sdla_Pressure_Johannes_2.mat'];
	load(exname);
	exname=[exdir2 'sdla_PoTemp_Johannes_2.mat'];
	load(exname);
case 'rho'
    exname=[exdir2 'sdla_Rho_Johannes.mat'];
end