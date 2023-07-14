%Pick a sub-region if desired
LAT_val = [45 80]; LON_val = [-60 30]; %Iceland region in Florent's plots
    %to avoid the rims
    
gcm_str = 'UM';

%--- Load and process the data
%dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/longer_timeseries/';
dirUM = '/home/disk/eos8/d.grosvenor/UM/Florent_Iceland/';

for idat=1:99
   flag{idat} = 'nc';
   fileUM_rho{idat} = ''; 
end

idat=1;

%fileUM{idat} = 'xkqkv_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'aoxaaa.pa20140901.pp.nc'; labs_UM(idat).l = 'Florent Iceland'; fileUM_W{idat} = 'xkqkv_W_3D_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;

idat_UM=1;
%Just working with output from one run for this
filename_W = [dirUM fileUM_W{idat_UM}]; 
filename = [dirUM fileUM{idat_UM}]; 

%read in UM data - for reading in one variable, set it to file_lwp
clear vars_in
vars_in.var = 'reff';
vars_in.flag = flag{idat_UM};
%vars_in.file_lwp = filename;
%vars_in.file_rho = filename_rho;
vars_in.file_lwp = filename;
vars_in.pole_lat = pole_lat;
vars_in.pole_lon = pole_lon;
vars_in.time_in = [];
vars_in.iz = -1; %Set to -1 for all levels


[re_temp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
vars_in.var = 'reff_weight';
[re_weight,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

%Divide reff by the weights
re_UM = re_temp ./ re_weight;


%% Tau
vars_in.var = 'COD';
[tau_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

%% Cloud top Nd for comparison with the calc
vars_in.var = 'Nd_ctop';
[Nd_ctop_UM_temp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
vars_in.var = 'Nd_ctop_weight';
[Nd_ctop_weight,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
%Divide by the weights
Nd_ctop_UM = Nd_ctop_UM_temp ./ Nd_ctop_weight;

%% Temperature and qcl to estimate cloud top temp
vars_in.var = 'T';
[T_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
vars_in.var = 'qcl';
[qcl_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

[max_qL,imax_qL]= max(qcl_UM,[],1);
for ilat=1:size(tau_UM,1)
    for ilon=1:size(tau_UM,2)
        ctt_UM(ilat,ilon) = T_UM(imax_qL(i),ilat,ilon);
    end
end
Tsurf = squeeze(T_UM(1,:,:));

%Restrict to clouds that are below say 3km e.g. Tsurf minus 27 degrees or
%so (9K per km, perhaps a bit less since may not be dry adiabatic)
%Probably not a good measure of CTH - the difference is negative in some
%places... perhaps not so important anyway
ihigh = find(Tsurf-ctt_UM>20);





N_UM = MODIS_justN_func(tau_UM,re_UM_temp*1e-6,'calc',[],ctt_UM,'N');
N_UM_low = N_UM;
N_UM_low(ihigh)=NaN;

%Restrict to specific region 
%[re_UM,iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges,iLAT,iLON,iT,iLAT2,iLON2,iT2] = restrict_array3D_to_region(re_UM_temp,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);
%ireg = find(gcm_Plat2D>=LAT_val(1) & gcm_Plat2D<=LAT_val(2) & gcm_Plon2D>=LON_val(1) & gcm_Plon2D<=LON_val(2));
ilat = find(gcm_Plat2D_UM(:,1)>=LAT_val(1) & gcm_Plat2D_UM(:,1)<=LAT_val(2));
ilon = find(gcm_Plon2D_UM(1,:)>=LON_val(1) & gcm_Plon2D_UM(1,:)<=LON_val(2));
re_UM_reg = re_UM(ilat,ilon);




