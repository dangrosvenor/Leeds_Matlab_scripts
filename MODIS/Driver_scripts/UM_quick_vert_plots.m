% um_str='u-aj091'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No-mix-tr-bl-210cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str='u-aj097'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
% 
% um_str='u-aj470'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-600cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str='u-aj471'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3-CS-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='o'; idat=idat+1;
% 
% um_str='u-aj472'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No-mix-tr-bl-210cm3-CS-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str='u-aj475'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3-rain-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.8 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
% 
% um_str='u-aj476'; fileUM{idat} = [um_str '/' um_str '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No-mix-tr-bl-210cm3-rain-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% 



run = 'u-aj471'; %CS off, procesing, rain on, 210 per cc
run = 'u-aj475'; %CS on, procesing, rain off, 210 per cc

it = 2;
iz = [1 21];


file_Nd=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_Nd_VOCALS_1p0_L70_ukv_.pp.nc'];
file_NR=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_NR_VOCALS_1p0_L70_ukv_.pp.nc'];
file_Naccum=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_accum_num_VOCALS_1p0_L70_ukv_.pp.nc'];
file_Maccum=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_accum_mass_VOCALS_1p0_L70_ukv_.pp.nc'];
file_Ncoarse=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_coarse_num_VOCALS_1p0_L70_ukv_.pp.nc'];
file_Mcoarse=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_coarse_mass_VOCALS_1p0_L70_ukv_.pp.nc'];
file_NR=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_NR_VOCALS_1p0_L70_ukv_.pp.nc'];
%file_qR=['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' run '/' run '_qR_VOCALS_1p0_L70_ukv_.pp.nc'];

z_um = netcdf_Dan2(file_Nd,'hybrid_ht');
Nd_02 = netcdf_Dan2(file_Nd,'Nd',[1 1 iz(1) it],[Inf Inf iz(end)-iz(1)+1 1]);
NR_02 = netcdf_Dan2(file_NR,'NR',[1 1 iz(1) it],[Inf Inf iz(end)-iz(1)+1 1]);
Naccum_02 = netcdf_Dan2(file_Naccum,'accum_num',[1 1 iz(1) it],[Inf Inf iz(end)-iz(1)+1 1]);
Maccum_02 = netcdf_Dan2(file_Maccum,'accum_mass',[1 1 iz(1) it],[Inf Inf iz(end)-iz(1)+1 1]);
Ncoarse_02 = netcdf_Dan2(file_Ncoarse,'coarse_num',[1 1 iz(1) it],[Inf Inf iz(end)-iz(1)+1 1]);
Mcoarse_02 = netcdf_Dan2(file_Mcoarse,'coarse_mass',[1 1 iz(1) it],[Inf Inf iz(end)-iz(1)+1 1]);

figure; pcolor([1:600],z_um(iz(1):iz(end)),squeeze(Nd_02(:,300,:)+NR_02(:,300,:)+Naccum_02(:,300,:)+Ncoarse_02(:,300,:))); colorbar; shading flat
figure; pcolor([1:600],z_um(iz(1):iz(end)),squeeze(Maccum_02(:,300,:)+Mcoarse_02(:,300,:))); colorbar; shading flat
