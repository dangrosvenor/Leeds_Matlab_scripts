%load in the cross sections taken along the streamlines

file_dir2 = 'C:\Documents and Settings\dan\My Documents\WRF\';

file_case=2;
switch file_case
    case 1 %LAT=67.5, most northern streamline
        file_dir=[file_dir2 'ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\'];
        file='pot_slice_streamline_LAT=67.5_z0=950_str_succ=129.mat';
    case 2 %LAT=68.1
        file_dir=[file_dir2 'ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\'];
        file='pot_slice_streamline_LAT=68.1_z0=950_str_succ=149.mat';
    case 3 %LAT=68.6
        file_dir=[file_dir2 'ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\'];
        file='pot_slice_streamline_LAT=68.6_z0=1200_str_succ=145.mat';
    case 4
        file_dir = 'Y:\WRF\ecmwf_ml_0.5_nudging\figures_for_paper_Aug2010\';
        file='pot_slice_streamline_LAT=67.5_z0=950_str_succ=129_0-5km.mat';
end
%N.B. the str_succ numbers above will only apply if use the same
%nstream=200 value in streamlines_threeD_draw.m

filename=[file_dir file];
whos('-file',filename)
load(filename);

%calculate the terrain baesed on the locations of the NaNs in the cross section
for i=1:size(pot_cross_6th_Jan_6UTC,2)
    inan=isnan(pot_cross_6th_Jan_6UTC(:,i));
    iterr=find(inan==0);
    terr_cross(i)=XY_pot_cross_data.Y_cross(min(iterr));
end

%subtract the miniumum - not sure if this is because of the terrain being above sea-level?
terr_cross2 =terr_cross - min(terr_cross);
disp('Loaded cross section data');
