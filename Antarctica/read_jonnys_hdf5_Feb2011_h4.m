irecon=0;

file_name_h5='flight99\ice_size_dists_etc_flight99_Post_Jan2011_MAX_DOF';
%%%%file_name_h5='flight100\ice_size_dists_etc_flight100'; %NOT RIGHT - without max dof

file_name_h5='flight100\ice_size_dists_etc_flight100_Post_Jan2011_MAX_DOF';
%file_name_h5='flight101\ice_size_dists_etc_flight101_Post_Jan2011_MAX_DOF';
%file_name_h5='flight101\ice_size_dists_etc_flight101_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75';
%file_name_h5='flight102\ice_size_dists_etc_flight102_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75';
%file_name_h5='flight102\ice_size_dists_etc_flight102_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75_closed_particles_3';
%file_name_h5='flight102\ice_size_dists_etc_flight102_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75_RECON_SAMPLE_AREA';
%file_name_h5='flight117\ice_size_dists_etc';
%file_name_h5='flight104\ice_size_dists_etc_flight104_Post_Jan2011_MAX_DOF';
%file_name_h5='flight104\ice_size_dists_etc_flight104_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75';
%file_name_h5='flight104\ice_size_dists_etc_flight104_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75_TEST';
%file_name_h5='flight104\ice_size_dists_etc_flight104_Post_Jan2011_MAX_DOF_multRej_9_2ndary_area_0.75_TEST - CENTRE IN';irecon=1;
%file_name_h5='flight104\ice_size_dists_etc_flight104_Post_Jan2011_MAX_DOF_CENTRE-IN';irecon=1;
%file_name_h5='flight104\test';irecon=1;

%file_name_h5='flight105\ice_size_dists_etc_flight105';
%file_name_h5='flight108\ice_size_dists_etc_flight108';
%file_name_h5='flight113\ice_size_dists_etc_flight113';

%file_name_h5='flight117\ice_size_dists_etc_flight117';

%file_name_h5='flight120\ice_size_dists_etc_flight120';
%file_name_h5='flight122\ice_size_dists_etc_flight122';
%file_name_h5='flight123\ice_size_dists_etc_flight123';
%file_name_h5='flight129\ice_size_dists_etc_flight129';

if ~exist('comp')
    comp='uni';
    comp='laptop';
    comp='UWchallenger';
end
if ~exist('cut_off_size')
    cut_off_size=1;
end
if ~exist('CAS_LWC_cut_off_sizes')
    CAS_LWC_cut_off_sizes=[0 50];
end

switch comp
case 'uni'
    filedir='Y:\BAS_flights\';
    filename_h5=[filedir file_name_h5 '.h5'];
case 'laptop'
    filedir='C:\Documents and Settings\G\My Documents\logbook\Antarctica\Flights and instruments\';
    filename_h5=[filedir file_name_h5 '.h4'];
case 'UWchallenger'
    filedir='/home/disk/eos1/d.grosvenor/Ant_flights/';
    ifind=strfind(file_name_h5,'\');
    file_name_h5(ifind)='/';
    filename_h5=[filedir file_name_h5 '.h5'];    
end


switch comp
case {'uni','UWchallenger'}
    
    INFO = Hdf5INFO(filename_h5)
    
    INFO.GroupHierarchy.Datasets.Name %gives the info of all the dataset names
    
    %INFO.GroupHierarchy.Datasets(1)  %this gives all the info about dataset 1
    
    NVARS = length(INFO.GroupHierarchy.Datasets);
    %if strcmp(INFO.GroupHierarchy.Datasets(11).Name,'/Time_UTC')==1
    %    CIP_size_bins_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(10)); %this retrieves dataset 10
    %    CIP_time_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(11)); %this retrieves dataset 11
    %else
    CIP_size_bins_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(NVARS-1)); %this retrieves dataset 10
    CIP_time_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(NVARS)); %this retrieves dataset 11
    %end
    
    if NVARS>=17
        CIP_Naccept_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(10));    
        CIP_N_edge_Rej_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(11));
        CIP_N_IAT_Rej_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(12));
        CIP_N_mult_frames_Rej_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(13));
        CIP_N_secondary_image_area_Rej_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(14));    
        CIP_N_thresh_Rej_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(15));   
        CIP_N_nan_Rej_Jonny = hdf5read(INFO.GroupHierarchy.Datasets(16));   
        
        CIP_NRej_Jonny = CIP_N_edge_Rej_Jonny + CIP_N_IAT_Rej_Jonny + CIP_N_thresh_Rej_Jonny + ...
            CIP_N_mult_frames_Rej_Jonny + CIP_N_secondary_image_area_Rej_Jonny+CIP_N_nan_Rej_Jonny;
        
    end
    
    
    ice_PSD_read = hdf5read(INFO.GroupHierarchy.Datasets(1)); %this retrieves dataset 1
    round_PSD_read = hdf5read(INFO.GroupHierarchy.Datasets(4)); 
    small_PSD_read = hdf5read(INFO.GroupHierarchy.Datasets(7)); 
    
    ice_no_Jonny_read = hdf5read(INFO.GroupHierarchy.Datasets(2)); %
    round_no_Jonny_read = hdf5read(INFO.GroupHierarchy.Datasets(5)); %
    small_no_Jonny_read = hdf5read(INFO.GroupHierarchy.Datasets(8)); %
    ice_mass_Jonny_read = hdf5read(INFO.GroupHierarchy.Datasets(3)); %
    round_mass_Jonny_read = hdf5read(INFO.GroupHierarchy.Datasets(6)); %
    small_mass_Jonny_read = hdf5read(INFO.GroupHierarchy.Datasets(9)); %
    
    
case 'laptop'  %do a h4 read instead as Matlab R12 doesn't do h5 - need to convert from h5 to h4 using the h5toh4 utility


    SD_id = hdfsd('start',filename_h5,'read'); %open the file
    [NVARS,nglobal_attr,status] = hdfsd('fileinfo',SD_id); %get info on the data in this file
    NVARS=NVARS-1;
    
    CIP_size_bins_Jonny = read_h4_var(SD_id,NVARS-1);   
    CIP_time_Jonny = read_h4_var(SD_id,NVARS);
    
    if NVARS>=17
        CIP_Naccept_Jonny = read_h4_var(SD_id,9);
        CIP_N_edge_Rej_Jonny = read_h4_var(SD_id,10);
        CIP_N_IAT_Rej_Jonny = read_h4_var(SD_id,11);
        CIP_N_mult_frames_Rej_Jonny = read_h4_var(SD_id,12);
        CIP_N_secondary_image_area_Rej_Jonny = read_h4_var(SD_id,13);
        CIP_N_thresh_Rej_Jonny = read_h4_var(SD_id,14);
        CIP_N_nan_Rej_Jonny = read_h4_var(SD_id,15);
        
        CIP_NRej_Jonny = CIP_N_edge_Rej_Jonny + CIP_N_IAT_Rej_Jonny + CIP_N_thresh_Rej_Jonny + ...
            CIP_N_mult_frames_Rej_Jonny + CIP_N_secondary_image_area_Rej_Jonny+CIP_N_nan_Rej_Jonny;
        
    end
    
    
    ice_PSD_read = read_h4_var(SD_id,0);
    round_PSD_read = read_h4_var(SD_id,3);
    small_PSD_read = read_h4_var(SD_id,6);
    
    ice_no_Jonny_read = read_h4_var(SD_id,1);
    round_no_Jonny_read = read_h4_var(SD_id,4);
    small_no_Jonny_read = read_h4_var(SD_id,7);
    ice_mass_Jonny_read = read_h4_var(SD_id,2);
    round_mass_Jonny_read = read_h4_var(SD_id,5);
    small_mass_Jonny_read = read_h4_var(SD_id,8);
end





%NOTE - need to change the ice size distribution script in Jonny's code so that 
%it takes into account the airspeed as it alters the sample volume and so
%the number concentrations

%to a good approximation the CIP sample volume is linearly dependent on airspeed for all size bins
%e.g. try:
for i=1:62
    sVol2(i)=SampleVolumeCIP(100,i,1); %100 m/s
end
for i=1:62
    sVol(i)=SampleVolumeCIP(50,i,1); %50 m/s
end

%sVol./sVol2; %shows that the ratio is almost constant over the size bins for the two airspeeds

% function sVol=SampleVolumeCIP(tasVal,bin,sample_time)  
% %function sVolGetSampleVolumeCIP(tasVal,bin,sample_time)  
% %sVol is the sample volume per second of sample time in cm^3
% %tasVal is the true airspeed (m/s)
% %bin is the width of the particles in number of bins minus 1
% %sample_time is the length of time that the sample was taken for in seconds

%so should be able to just scale numbers by the true airspeed/Jonny's airspeed
%Jonny assumed an airspeed of 60 m/s? Actually, think he may have used the CIP airspeed...
%as TAS_Source is set to one, which I think means it uses the CIP airspeed
%actually it was using a constant speed of 60 m/s as that bit of code wasn't being used

%Jonny added an output variable, which is present for some of the files, but not others
%CIP size bins and the time are always the last two, though

air_speed_type='CIP probe';

%  [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number{icas_count}...
%                 ,CIP_total_number{icas_count},LWC_dist_cas,LWC_dist_cip...
%                 ,CAS_mode_diameter,CAS_mean_diameter,LWC_dist_cas_cutoff]...
%                 =cas_sample_volume_and_stats(dat_flt,CAS_time_all,...
%                 CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all...
%                 ,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes);  

[sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
        ,CAS_total_number_cutoff ...
    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
        ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
    =cas_sample_volume_and_stats2...
    (dat_flt,CAS_time_all,...
    CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,0); 


number_factor = 60./(air_speed_1D/100); %multiply by this to get the corrected numbers
%sample volume is proportional to the airspeed so number is inverslely proportional
%number_factor is based on CAS time (it is interpolated from CIP_time_all to CAS_time_all in 
%cas_sample_volume_and_stats

number_f2 = interp1(CAS_time_all,number_factor,CIP_time_Jonny); %will have lots of NaNs


TAS_Jonny = interp1(CIP_time_all,TAS_all,CIP_time_Jonny); %airspeed on Jonny time

clear sVol_Jonny
for i=0:size(ice_PSD_read,1)-1   %61
    if irecon==1
        disp('**** WARNING - USING RECONSTRUCTED SAMPLE VOLUME *****');
        sVol_Jonny(i+1,:)=SampleVolumeCIP_recon(TAS_Jonny,i,1); %TAS_all is the airspeed in m/s        
    else
        sVol_Jonny(i+1,:)=SampleVolumeCIP(TAS_Jonny,i,1); %TAS_all is the airspeed in m/s
        %on the same time base as the CIP counts
    end
end



number_f2_rep = repmat(number_f2,[1 size(ice_PSD_read,1)])';
ice_PSD = ice_PSD_read.*number_f2_rep;
round_PSD = round_PSD_read.*number_f2_rep;
small_PSD = small_PSD_read.*number_f2_rep;
tot_PSD = ice_PSD + round_PSD + small_PSD;






ice_no_Jonny = number_f2(1:end-1).*ice_no_Jonny_read; %
round_no_Jonny = number_f2(1:end-1).*round_no_Jonny_read; %
small_no_Jonny = number_f2(1:end-1).*small_no_Jonny_read; %
ice_mass_Jonny = number_f2(1:end-1).*ice_mass_Jonny_read; %
round_mass_Jonny = number_f2(1:end-1).*round_mass_Jonny_read; %
small_mass_Jonny = number_f2(1:end-1).*small_mass_Jonny_read; %

ice_no_particles_PSD = ice_PSD .* sVol_Jonny;
round_no_particles_PSD = round_PSD .* sVol_Jonny;
small_no_particles_PSD = small_PSD .* sVol_Jonny;
tot_no_particles_PSD = ice_no_particles_PSD + round_no_particles_PSD + small_no_particles_PSD;

ice_no_total_Jonny = ice_no_Jonny + round_no_Jonny + small_no_Jonny;

total_mass_Jonny = ice_mass_Jonny + round_mass_Jonny + small_mass_Jonny;

CIP_mid = 0.5*(CIP_size_bins_Jonny(1:end-1) + CIP_size_bins_Jonny(2:end) );
CIP_mid_rep = repmat(CIP_mid,[1 size(ice_PSD,2)]);
ice_no_CIP_Dan = sum(ice_PSD,1);

mean_ice_size = sum(ice_PSD.*CIP_mid_rep,1) ./ ice_no_CIP_Dan;
mean_ice_size(ice_no_CIP_Dan==0)=0;



ice_no_tot_CIP_Dan = sum(ice_PSD+round_PSD+small_PSD,1);

CIP_time_Jonny2=CIP_time_Jonny(1:end-1); %tested out which times correspond to the correct numbers
%in ice_no_CIP_Dan (as calculated from the PSDs) and it is 1:end-1 rather than 2:end, at least for flt117


disp('Done read Jonny');

%1 /CIP_Ice_PSD_cm3  [size * time]
%2 /CIP_Ice_cm3      [time]
%3 /CIP_Ice_gm3
%4 /CIP_Round_PSD_cm3
%5 /CIP_Round_cm3
%6 /CIP_Round_gm3
%7 /CIP_Small_PSD_cm3
%8 /CIP_Small_cm3
%9 /CIP_Small_gm3
%10 /Size_bins_um
%11 /Time_UTC