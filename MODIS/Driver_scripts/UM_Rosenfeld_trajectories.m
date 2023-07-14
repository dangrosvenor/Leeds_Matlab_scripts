%Look at trajectories endig along 20S from model - using Ryan's
%trajectories

%Will look at ones ending at 75:5:100 W (20S). Require at least a certain
%time before this, though.
% Ryan's trajectories are every 12 hours with 7 points [0 12 24 36 48 60
% 72] hours from start of traj.

%Ryan's trajectories converted to Matlab time
load_file='/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_matlab_dates.mat';
%traj_dat = load(load_file);

% whos('-file',load_file)
%   Name           Size                 Bytes  Class     Attributes
% 
%   id_traj        1x35087             280696  double              
%   lats_traj      7x35087            1964872  double              
%   lons_traj      7x35087            1964872  double              
%   time_traj      7x35087            1964872  double  

%values in trajectories run from 0 to 360
lon_target = [-100:5:-75]+360; 
lat_target = -20 * ones(size(lon_target));
tol = 1; %+/- tolerance for the match
it_min = 1;

clear ii n_traj tot_lat tot_lon npoints
for i=1:length(lon_target)
    ii{i} = find(abs(traj_dat.lons_traj - lon_target(i)) <= tol & abs(traj_dat.lats_traj - lat_target(i)) <= tol ); 
    n_traj(i)=length(ii{i});
    [itime{i},itraj{i}] = ind2sub(size(traj_dat.lons_traj),ii{i});
    ivalid = find(itime{i}>it_min); %don't want ones that end at the target location
    nT = size(traj_dat.lons_traj,1);
    tot_lat{i} = zeros([nT 1]);
    tot_lon{i} = zeros(size(tot_lat{i}));
    npoints{i} = zeros(size(tot_lat{i}));    
    for i2=1:length(ivalid)
       it = ivalid(i2);
       istart = nT - itime{i}(it) + 1;       
       tot_lat{i}(istart:end) = tot_lat{i}(istart:end) +  traj_dat.lats_traj(1:itime{i}(it) , itraj{i}(it) );
       tot_lon{i}(istart:end) = tot_lon{i}(istart:end) +  traj_dat.lons_traj(1:itime{i}(it) , itraj{i}(it) );
       npoints{i}(istart:end) = npoints{i}(istart:end) + 1;
    end
    npoints{i}(npoints{i}==0) = NaN;
    mean_lat{i} = tot_lat{i} ./ npoints{i};
    mean_lon{i} = tot_lon{i} ./ npoints{i};
end

%% plot the mean trajectories

figure





