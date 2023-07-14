function []=RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,dat)
%LWP_PD_ALL, RWP_PD_ALL,...
%TLWP_PD_ALL,Conv_LWP_PD_ALL,Conv_RWP_PD_ALL,low_CF_PD_ALL,total_CF_PD_ALL,Nd_PD_ALL,max_cloud_height_in_cloud_LWC_IWC_PD)

%Convert all of the variables in the input structure into variables.
name_struc='dat'; %The name of the structure
names = eval(['fieldnames(' name_struc ');']);
for i=1:length(names)
    eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
    eval(eval_str);
end

%Set time tolerance within which they have to match the trajetory time
%time_tol_traj = 30/24; %in days
time_tol_traj = 3.5/24; %in days

%File below created using RyanE_traj_analysis.m - it's the times of Ryan's
%trajectories in Matlab time.
load_file='/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_matlab_dates.mat';
load(load_file);

%save_file='/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data.mat';


%Can use interpn because the grids are regular in lat, lon and time

%Example using all lat and lons offset by 0.2
x = gcm_Plat2D_UM(:,1);
y = gcm_Plon2D_UM(1,:)';
y(y<0)=y(y<0)+360; %ensures that it is monotonic (runs from 0 to 360).
[x3,y3,t3]=meshgrid(x,y,time_100); %Think can used ndgrid to avoid the need to permute
x3=permute(x3,[2 1 3]);
y3=permute(y3,[2 1 3]);
t3=permute(t3,[2 1 3]);


%traj_interp_temp = interpn(x3,y3,t3,dat,lats_traj(:),lons_traj(:),time_traj(:),'nearest');
%traj_interp_interpn = reshape(traj_interp_temp ,size(lats_traj));
%traj_interp = griddata(x3,y3,t3,Nd_PD_ALL(:,:,1:3),20,100,time_100(1)+1.1250);

%Test using own find nearest code since then get the difference between the
%nearest value and that requested.

%Match lat,lon and time separately
clear dlat_match ilat_match dlon_match ilon_match dtime_match itime_match time_match lat_match lon_match
for i=1:length(lats_traj(:))
    dlat = abs(lats_traj(i)-x);
    [dlat_match(i),ilat_match(i)]=min(dlat);
    lat_match(i) = x(ilat_match(i));
    
    dlon = abs(lons_traj(i)-y);
    [dlon_match(i),ilon_match(i)]=min(dlon);
    lon_match(i) = y(ilon_match(i));
    
    dtim = abs(time_traj(i)-time_100);
    [dtime_match(i),itime_match(i)]=min(dtim);
    time_match(i) = time_100(itime_match(i));
end

dtime_match = reshape(dtime_match,size(lats_traj)); %dtime is in days
dlon_match = reshape(dlon_match,size(lats_traj)); %dtime is in days
dlat_match = reshape(dlat_match,size(lats_traj)); %dtime is in days
time_match = reshape(time_match,size(lats_traj)); %dtime is in days
lon_match = reshape(lon_match,size(lats_traj)); %dtime is in days
lat_match = reshape(lat_match,size(lats_traj)); %dtime is in days

% Find the points that do not match within the required tolerance
inan=find(dtime_match>time_tol_traj | isnan(dtime_match)==1 );

%Now use the indices to get the required data at those points.
inds_match = sub2ind([length(x) length(y) length(time_100)],ilat_match,ilon_match,itime_match);


TLWP_PD_ALL = LWP_PD_ALL + RWP_PD_ALL; %LS LWP + RWP
TLWP_Conv_PD_ALL = TLWP_PD_ALL + Conv_LWP_PD_ALL + Conv_RWP_PD_ALL; %LS + conv TLWP

%Extract the data
traj_low_CF=reshape(low_CF_PD_ALL(inds_match),size(lats_traj)); traj_low_CF(inan)=NaN;
traj_total_CF=reshape(total_CF_PD_ALL(inds_match),size(lats_traj)); traj_total_CF(inan)=NaN;
traj_Nd=reshape(Nd_PD_ALL(inds_match),size(lats_traj)); traj_Nd(inan)=NaN;
traj_CTH=reshape(max_cloud_height_in_cloud_LWC_IWC_PD(inds_match),size(lats_traj)); traj_CTH(inan)=NaN;
traj_LWP=reshape(LWP_PD_ALL(inds_match),size(lats_traj)); traj_LWP(inan)=NaN;
traj_TLWP=reshape(TLWP_PD_ALL(inds_match),size(lats_traj)); traj_TLWP(inan)=NaN;
traj_TLWP_inc_Conv=reshape(TLWP_Conv_PD_ALL(inds_match),size(lats_traj)); traj_TLWP_inc_Conv(inan)=NaN;
traj_BL_height_combined_ALL=reshape(BL_height_combined_ALL(inds_match),size(lats_traj)); traj_BL_height_combined_ALL(inan)=NaN;
traj_BL_height_dthdz_ALL=reshape(BL_height_dthdz_ALL(inds_match),size(lats_traj)); traj_BL_height_dthdz_ALL(inan)=NaN;
traj_tot_surf_precip_rate_ALL=reshape(tot_surf_precip_rate_ALL(inds_match),size(lats_traj)); traj_tot_surf_precip_rate_ALL(inan)=NaN;
traj_liqCF_COSP_MODIS_PD_ALL=reshape(liqCF_COSP_MODIS_PD_ALL(inds_match),size(lats_traj)); traj_liqCF_COSP_MODIS_PD_ALL(inan)=NaN;
traj_totCF_COSP_MODIS_PD_ALL=reshape(totCF_COSP_MODIS_PD_ALL(inds_match),size(lats_traj)); traj_totCF_COSP_MODIS_PD_ALL(inan)=NaN;

traj_LTS=reshape(LTS_PD_ALL(inds_match),size(lats_traj)); traj_LTS(inan)=NaN;
traj_omega700=reshape(omega700_PD_ALL(inds_match),size(lats_traj)); traj_omega700(inan)=NaN;
traj_qv700=reshape(qv700_PD_ALL(inds_match),size(lats_traj)); traj_qv700(inan)=NaN;
traj_surface_temp=reshape(surface_temp_PD_ALL(inds_match),size(lats_traj)); traj_surface_temp(inan)=NaN;
traj_U10=reshape(U10_PD_ALL(inds_match),size(lats_traj)); traj_U10(inan)=NaN;
traj_V10=reshape(V10_PD_ALL(inds_match),size(lats_traj)); traj_V10(inan)=NaN;



%% Save the data for Ryan
notes.general=['Model values are interpolated to trajectory points using ''nearest'' interpolation method. They are save as traj_<var>. '...
    'Values for which the time difference between the model output and the trajectory is larger than +/-' num2str(24*time_tol_traj) ' hours are set to NaN. '...
    'Vaulues of CTH and Nd are also set to NaN when there is no cloud.'];

notes.dtime_match = 'abs difference between the model time point and the actual trajectory point in days';
notes.dlon_match = 'abs difference between the model lon point and the actual trajectory point in degrees';
notes.dlat_match = 'abs difference between the model lat point and the actual trajectory point in degrees';
notes.time_match = 'Closest mode time point to the trajectory point (Matlab time convention in days)';
notes.lon_match = 'Closest mode lon point to the trajectory point (degrees)';
notes.lat_match = 'Closest mode lat point to the trajectory point (degrees)';

notes.id_traj = 'ID number of trajectory';
notes.traj_low_CF = 'Low cloud (pressure > 680 hPa) fraction';
notes.total_CF = 'Total random max overlap cloud fraction';
notes.traj_Nd = 'Droplet number conc (m^{-3})';
notes.traj_CTH = 'Cloud top height (m) calculated as max height for which the in-cloud LWC+IWC is greater than 0.05 g/kg';
notes.traj_LWP = 'LWP (g/m2) from the model large scale cloud scheme';
notes.traj_TLWP = 'LWP + RWP (g/m2) from the large scale cloud scheme';
notes.traj_TLWP_inc_Conv = 'LWP+RWP (g/m2) from convection scheme and large scale scheme';
notes.traj_BL_height_combined_ALL = 'Top of boudary layer determined from BL diag 3-304 (turb. mixing height) for all BL types except 6 (cumulus) for which use 3-359 (top of parcel)';
notes.traj_BL_height_dthdz_ALL = 'Top of boundary layer determined from max of theta gradient below 7500m';
notes.traj_tot_surf_precip_rate_ALL = 'Total (large-scale + convective, rain+snow) surface precip rate (kg/m2/s)';
notes.traj_liqCF_COSP_MODIS_PD_ALL = 'COSP satellite simulator MODIS liquid cloud fraction (from UKMET model)';
notes.traj_totCF_COSP_MODIS_PD_ALL = 'COSP satellite simulator MODIS total cloud fraction (from UKMET model)';

notes.traj_LTS = 'LTS calculated as potemp at 700mb minus that at 1000mb (K)';
notes.traj_omega700 = 'Subsidence at 700 mb. Calculated as -rho700*9.81*w-component wind at 700 mb (Pa/s)';
notes.traj_qv700 = 'Specific humidity at 700 mb (kg/kg)';
notes.traj_surface_temp = 'Surface temperature (SST over ocean) (K)';
notes.traj_U10 = 'U-wind speed at 10m from surface (m/s)';
notes.traj_V10 = 'V-wind speed at 10m from surface (m/s)';

save(save_file,'notes','id_traj','dtime_match','dlon_match','dlat_match','time_match','lon_match','lat_match','traj_low_CF','traj_total_CF',...
    'traj_Nd','traj_CTH','traj_LWP','traj_TLWP','traj_TLWP_inc_Conv','traj_BL_height_combined_ALL', ...    
    'traj_BL_height_dthdz_ALL', 'traj_tot_surf_precip_rate_ALL', 'traj_liqCF_COSP_MODIS_PD_ALL', 'traj_totCF_COSP_MODIS_PD_ALL', ...
    'LWP_PD_ALL', 'RWP_PD_ALL','TLWP_PD_ALL','Conv_LWP_PD_ALL','Conv_RWP_PD_ALL','low_CF_PD_ALL','total_CF_PD_ALL','Nd_PD_ALL', ...
    'max_cloud_height_in_cloud_LWC_IWC_PD','BL_height_combined_ALL', 'BL_height_dthdz_ALL',...
    'tot_surf_precip_rate_ALL', 'liqCF_COSP_MODIS_PD_ALL', 'totCF_COSP_MODIS_PD_ALL', ...
    'traj_LTS','traj_omega700','traj_qv700','traj_surface_temp','traj_U10','traj_V10',...
    'LTS_PD_ALL','omega700_PD_ALL','qv700_PD_ALL','surface_temp_PD_ALL','U10_PD_ALL','V10_PD_ALL');


iplot_histo=0;
if iplot_histo==1

%% Calculate the number of complete trajectories (no NaNs in whole
%trajectory)
%size(traj_interp2)= [7 35087]

inan=isnan(traj_Nd);
%inan=isnan(traj_low_CF);
num_nans = sum(inan,1);
iwhole = find(num_nans==0);
Ntraj = length(iwhole);

qh=ndhistc_run([num_nans(:)],[-0.5:1:7.5]);
figure
plot([0:7],qh,'bx-');
xlabel('Number of NaNs in each trajectory');
ylabel('Frequency');

end




