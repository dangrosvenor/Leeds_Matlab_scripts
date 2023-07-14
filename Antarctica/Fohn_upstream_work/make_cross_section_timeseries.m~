%% read in the filenames of the cross sections that we want to go into the
%% timeseries
clear dir U_mean_timser N_Brunt_timser F0_timser Hbar_timser

direc = '/home/disk/eos1/d.grosvenor/matlab/work/Antarctica/Fohn_upstream_work/';
list_of_pot_files = dir([direc 'Pot*20130808*']);
list_of_U_files = dir([direc 'Component*20130808*']);
list_of_WindDir_files = dir([direc 'Wind direction*20130808*']);

%% load them in one at a time and extract the Froude number and U and N
%% from the first (far left) column in the cross section

for itimser_F0=1:length(list_of_pot_files)

override_load_cross_section=1;
pot_filename = [direc list_of_pot_files(itimser_F0).name];
U_filename = [direc list_of_U_files(itimser_F0).name];
WindDir_filename = [direc list_of_WindDir_files(itimser_F0).name];

day(itimser_F0) = str2num(list_of_U_files(itimser_F0).name(41:42));
hour(itimser_F0) = str2num(list_of_U_files(itimser_F0).name(48:49));
min(itimser_F0) = str2num(list_of_U_files(itimser_F0).name(50:51));



load_cross_section_data_July2013

U_mean_timser(itimser_F0) = U_mean;
N_Brunt_timser(itimser_F0) = N_Brunt;
F0_timser(itimser_F0) = F0;
Hbar_timser(itimser_F0) = Hbar;



end

matlab_time = datenum(2006,1,day,hour,min,0);

