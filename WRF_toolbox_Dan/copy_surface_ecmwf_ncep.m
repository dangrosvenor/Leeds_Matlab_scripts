%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_sfRUC_v3/';
dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_sfRUC_v3/';
dir_ecmwf = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08/';
%dir_ncep = '/home/mbexddg5/mobile4/Database4/Dans_Model/WRF/ant_jan06_ecmwf_Nov08_met_em/';

%file_ncep = 'met_em_d03-06_18-00-00';
%file_ecmwf = 'met_em_d03-01-06_18_ecmwf';

files = dir([dir_ecmwf '/met_em*']);

for i=1:length(files)
	file_ncep = [dir_ncep files(i).name]
	file_ecmwf = [dir_ecmwf files(i).name]
	nc_ncep = netcdf(file_ncep); %open ncep file for reading
	nc_ecmwf = netcdf(file_ecmwf,'write'); %open ecmwf file for writing

	copy(nc_ncep{'ST010200'},nc_ecmwf,1,1);
	copy(nc_ncep{'ST000010'},nc_ecmwf,1,1);
	copy(nc_ncep{'SM010200'},nc_ecmwf,1,1);
	copy(nc_ncep{'SM000010'},nc_ecmwf,1,1);

	name(nc_ecmwf{'SM'},'DELETED_SM');
	name(nc_ecmwf{'ST'},'DELETED_ST');
	copy(nc_ncep{'SM'},nc_ecmwf,1,1);
	copy(nc_ncep{'ST'},nc_ecmwf,1,1);

%	copy(nc_ncep{'SKINTEMP'},nc_ecmwf,1,1);
%	copy(nc_ncep{'SOILCBOT'},nc_ecmwf,1,1); %same
%	copy(nc_ncep{'SOILCTOP'},nc_ecmwf,1,1); %same
%	copy(nc_ncep{'SOILTEMP'},nc_ecmwf,1,1); %same
%	copy(nc_ncep{'LU_INDEX'},nc_ecmwf,1,1); %same
%	copy(nc_ncep{'LANDUSEF'},nc_ecmwf,1,1); %same
%	copy(nc_ncep{'LANDMASK'},nc_ecmwf,1,1); %same
	copy(nc_ncep{'SOILHGT'},nc_ecmwf,1,1);  %large diffs
	copy(nc_ncep{'LANDSEA'},nc_ecmwf,1,1);  %large diffs - ecmwf looks better but includes the ice shelf - ncep doesn't
	copy(nc_ncep{'SEAICE'},nc_ecmwf,1,1); %bit different
	copy(nc_ncep{'SLOPECAT'},nc_ecmwf,1,1); %large diffs
	copy(nc_ncep{'SNOALB'},nc_ecmwf,1,1); %large diffs
	copy(nc_ncep{'GREENFRAC'},nc_ecmwf,1,1); %large diffs
	copy(nc_ncep{'ALBEDO12M'},nc_ecmwf,1,1); %large diffs
	copy(nc_ncep{'SPECHUMD'},nc_ecmwf,1,1); %not present for ECMWF
%	copy(nc_ncep{'SLPY'},nc_ecmwf,1,1); %
%	copy(nc_ncep{'SLPX'},nc_ecmwf,1,1); %both X and Y very similar
%	copy(nc_ncep{'HGT_M'},nc_ecmwf,1,1); %same
	copy(nc_ncep{'SNOW'},nc_ecmwf,1,1); %very different - max 0f 10 kg/m2 for ecmwf and 48000 for ncep!! - ncep sounds more sensible -  48,000 = 48 m depth, 10 = 1 cm!

	name(nc_ecmwf{'SM100255'},'DELETED_1');
	name(nc_ecmwf{'SM028100'},'DELETED_2');
	name(nc_ecmwf{'SM007028'},'DELETED_3');
	name(nc_ecmwf{'SM000007'},'DELETED_4');
	name(nc_ecmwf{'ST100255'},'DELETED_5');
	name(nc_ecmwf{'ST028100'},'DELETED_6');
	name(nc_ecmwf{'ST007028'},'DELETED_7');
	name(nc_ecmwf{'ST000007'},'DELETED_8');
	name(nc_ecmwf{'SST'},'DELETED_9');

	close(nc_ncep);
	close(nc_ecmwf);

end




