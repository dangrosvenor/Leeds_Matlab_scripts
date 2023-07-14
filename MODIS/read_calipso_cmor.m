function out2=read_calipso_cmor(nc_dir,nc_inst_file,var_str,ilat,ilon,indices_str)

nc_inst=netcdf([nc_dir nc_inst_file],'nowrite');

out = eval(['nc_inst{''' var_str '''}' indices_str ';']); %
%double clhcalipso(time, lat, lon) ; _FillValue = 1.00000002004088e+20
out(out>1e19)=NaN;

%Move the data so that the longitudes are consistent with GCMs, etc.
out2(:,1:90)=out(:,91:180);
out2(:,91:180)=out(:,1:90);
