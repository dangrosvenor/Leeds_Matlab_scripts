traj_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories'];
fid=fopen(traj_file,'rt');
traj_dat=fscanf(fid,'%f',[30 Inf]);
fclose(fid);

% Data Format: 
% column 1: trajectory reference #
% column 2: year at trajectory start
% columns 3-9: observation day (julian day 1-365, for times 0:12:72)
% columns 10-16: observation UTC (0-2400)
% columns 17-23: latitude of observation
% columns 24-30: longitude of observation

%um_cal_res = griddata(gcm_Plat2D_UM,gcm_Plon2D_UM,um_data,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);

%Constuct the 7 times of each trajetory (0,12,24,36,48,60,72 hours plus
%start time)
years_traj = repmat(traj_dat(2,:),[7 1]);
days_traj = traj_dat(3:9,:);
hours_traj = traj_dat(10:16,:);
[date,date_num] = date_from_day_of_year_func(days_traj(:),years_traj(:)); %Matlab time in days
date_num = reshape(date_num,size(hours_traj));
time_traj = date_num + hours_traj/24;

lats_traj = traj_dat(17:23,:);
lons_traj = traj_dat(24:30,:);
i180=find(lons_traj<0);
lons_traj(i180) = lons_traj(i180) + 360;

id_traj = traj_dat(1,:);

save_file='/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_matlab_dates.mat';
save(save_file,'time_traj','lats_traj','lons_traj','id_traj');


%Can use interpn because the grids are regular in lat, lon and time
run_simple_example=0;
if run_simple_example==1
    %Example using all lat and lons offset by 0.2
    x = gcm_Plat2D_UM(:,1);
    y = gcm_Plon2D_UM(1,:)';
    y(y<0)=y(y<0)+360; %ensures that it is monotonic (runs from 0 to 360).
    [x3,y3,t3]=meshgrid(x,y,time_out(1:3)); %think using ndgrid would aviod needing to permute here
    x3=permute(x3,[2 1 3]);
    y3=permute(y3,[2 1 3]);
    t3=permute(t3,[2 1 3]);
    [xi,yi,ti]=meshgrid(x+0.2,y+0.2,time_out(1)+1.1250/2);
    %[xi,yi,ti]=meshgrid(x,y,time_out(1));
    xi=permute(xi,[2 1 3]);
    yi=permute(yi,[2 1 3]);
    ti=permute(ti,[2 1 3]);
    %dat=inpaint_nans(Nd_PD_ALL(:,:,1:3)); %only works with 2d arrays...
    dat = Nd_PD_ALL(:,:,1:3);
    traj_interp = interpn(x3,y3,t3,dat,xi,yi,ti,'nearest');
    %traj_interp = griddata(x3,y3,t3,Nd_PD_ALL(:,:,1:3),20,100,time_out(1)+1.1250);
end



