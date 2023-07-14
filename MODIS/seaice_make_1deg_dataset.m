% Also see 
%   seaice_match_times.m
% for the script that then matches the times to a MODIS dataset

%Automatically saves to filename stored in seaice_savefile
% This makes a 1x1 deg sea-ice dataset with both daily values and also one with
% the max value for a moving 2 week window (seaice_array_1deg_max)
% Looks like it just needs MLAT and MLON from the MODIS L3 data

%To load data ready for processing *** run read_sea_ice_multi ***, which will give seaice_array and
%seaice_array_datenum, and seaice_lat and seaice_lon
% seaice_lon,seaice_lat,seaice_array

%% Data produced from read_sea_ice_multi
% I.e. load the data specified here instead of running read_sea_ice_multi

% seaice_fileload_ALL = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_midDec2006_to_midJan2008_20140424T071237.mat'
% load(seaice_fileload_ALL)
%  (S. hemisphere for mid Dec 2006 to mid Jan 2008)

% S. hemisphere For mid Nov 2005 - mid Jan 2008 :-
%seaice_fileload_ALL = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_midNov2005_to_midJan2008_20140508T071821.mat';
%load(seaice_fileload_ALL);

%SHemisphere_midMay2008_to_midSep2008
%seaice_savefile = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/daily//saved_seaice_SHemisphere_midMay2008_to_midSep2008_20161221T091809.mat';


%% Data produced from this routine (seaice_make_1deg_dataset)
%OR load in the data for specific regions (and skip this routine):-
% 60-160 W, 60-70S :-
%     seaice_fileload_SO = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2007_60-140W_60-70S_20140425T051130.mat';
%     load(seaice_fileload_SO);

% Global, (SH only) 2007 :-
%     seaice_fileload_SO = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2007_global_20140508T014643.mat'
%      load(seaice_fileload_SO);

% Global (SH only) mid Nov 2005 - mid Jan 2008 :-
%   seaice_fileload_1deg_global ='/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2007_global_20140508T094510.mat';
%   load(seaice_fileload_1deg_global);

% NH (all lons) - 1st Nov 2006 to 31st Dec 2007
% seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2007_all_NHemisphere_20150529T103046.mat'
% load(seaice_fileload_1deg_NH)

% SH, JJA 2008
% /home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2008_JJA_all_NHemisphere_20161222T222241.mat
% NH, JJA 2008
% /home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2008_JJA_all_SHemisphere_20161221T095015.mat
% NH all 2000-2015
% /home/disk/eos8/d.grosvenor/sea_ice_data/north/daily//saved_seaice_NHemisphere_all_2000-2015_20161223T012957.mat



% N.B. will need to run 
%   seaice_match_times 
% afterwards to match the times to the MODIS data

seaice_savedir = '/home/disk/eos8/d.grosvenor/sea_ice_data/south';
%seaice_savedir = '/home/disk/eos8/d.grosvenor/sea_ice_data/north';

isave_seaice = 1;
tag = '2007_60-140W_60-70S';
tag = '2007_global';
tag = '2007_all_NHemisphere';
tag = '2007_all_NHemisphere';
tag = '2008_JJA_all_SHemisphere';
tag = '2008_JJA_all_NHemisphere';
tag = '2000-2015_all_NHemisphere';
tag = '2000-2015_all_SHemisphere';

MLAT = [-89.5 : 89.5];
MLON = [-179.5 : 179.5];

%% N.B. - dates don't have to match e.g. MODIS data at this stage - can
%% just process all of the data that has been loaded in and then cut down
%% later using seaice_match_times

% For these dates put the actual says that want (i.e. don't take worry about the 2 week window thing - just need to be sure have the data for 7 or 8 days either side)
%        start_date = datenum('01-Jan-2007');

%        start_date = datenum('01-Nov-2006');        
%        end_date = datenum('31-Dec-2007');        

        start_date = datenum('01-Jun-2008');
        end_date = datenum('31-Aug-2008'); 
        
        start_date = datenum('01-Jan-2001');
        end_date = datenum('31-Dec-2014');         
        
        ndays = end_date - start_date + 1;
        
        period = 15; %make this an odd number

%         seaice_array_max = NaN*ones([size(seaice_data) ndays]);
         seaice_array_1deg_datenum = NaN*ones([1 ndays]);
         
         seaice_array_1deg = NaN*ones([length(MLAT) length(MLON) ndays]);         
         seaice_array_1deg_max = NaN*ones([length(MLAT) length(MLON) ndays]);
        
        for isea2=1:ndays  %length(files)
            fprintf(1,'Processing day %i out of %i\n',isea2,ndays);
            
            day = start_date + isea2 - 1;
            day_minus = day - (period-1)/2;
            day_plus = day + (period-1)/2;
           
            iday = find(seaice_array_datenum == day);
            iminus = find(seaice_array_datenum == day_minus);
            iplus = find(seaice_array_datenum == day_plus);                        
            
            seaice_array_max = max( seaice_array(:,:,iminus:iplus) ,[],3 );
            seaice_array_1deg_datenum(isea2) = day;
            
            seaice_array_1deg(:,:,isea2) = griddata(seaice_lon,seaice_lat,seaice_array(:,:,iday),MLON,MLAT');
            seaice_array_1deg_max(:,:,isea2) = griddata(seaice_lon,seaice_lat,seaice_array_max,MLON,MLAT');


        end



if isave_seaice==1
   seaice_savefile = [seaice_savedir '/saved_seaice_1deg_' tag '_' datestr(now,30) '.mat']; 
   
   MLON_seaice = MLON;
   MLAT_seaice = MLAT;
 
   save(seaice_savefile,'MLON_seaice','MLAT_seaice','seaice_array_1deg','seaice_array_1deg_max','seaice_array_1deg_datenum','-V7.3');
 
end


