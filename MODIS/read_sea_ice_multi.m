%Also see 
% seaice_make_1deg_dataset.m 
%     for the script that makes the 1x1 deg gridded dataset

%Have some data saved following the running of this script to save loading
%in lots of data again:-

%N. Hemisphere mid Nov 2005 to mid Jan 2008
%file_seaice_load='/home/disk/eos8/d.grosvenor/sea_ice_data/north/daily/saved_seaice_NHemisphere_midNov2005_to_midJan2008_20150602T040143.mat';

%S. Hemisphere mid Dec 2006 to mid Jan 2008
% seaice_fileload_ALL = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_midDec2006_to_midJan2008_20140424T071237.mat'
% load(seaice_fileload_ALL)

%S. Hemisphere mid Nov 2005 - mid Jan 2008
%seaice_fileload_ALL = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_midNov2005_to_midJan2008_20140508T071821.mat';
%load(seaice_fileload_ALL);

%SHemisphere_midMay2008_to_midSep2008
%seaice_savefile = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/daily//saved_seaice_SHemisphere_midMay2008_to_midSep2008_20161221T091809.mat';

pole='Antarctic monthly';
%pole='Antarctic fortnightly';
pole='Antarctic daily';
%pole='Arctic';
%pole='Arctic daily';

isave_seaice = 1; %whether to save the seaice_array etc. after loading in - change the tag below
tag = 'NHemisphere_midNov2005_to_midJan2008';
tag = 'SHemisphere_midMay2008_to_midSep2008'; %i.e. JJA of 2008 with some extra for 2 week window
tag = 'NHemisphere_midMay2008_to_midSep2008'; %i.e. JJA of 2008 with some extra for 2 week window
tag = 'NHemisphere_all_2000-2015'; %all days for 2000-2015
tag = 'SHemisphere_all_2000-2015'; %all days for 2000-2015

isea=1;
clear seaice_file2

switch pole
    case 'Arctic'
        seaice_dir = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/daily/';
        
        seaice_file2{isea} = '2006/nt_20061201_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070101_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070201_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070301_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070401_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070501_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070601_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070701_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070801_f13_v01_n.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070901_f13_v01_n.bin'; isea=isea+1;        
        
    case 'Antarctic monthly'
        seaice_dir = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/daily/';
        
%        seaice_file2{isea} = '2006/nt_20061201_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070101_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070201_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070301_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070401_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070501_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070601_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20070701_f13_v01_s.bin'; isea=isea+1;        
        seaice_file2{isea} = '2007/nt_20070801_f13_v01_s.bin'; isea=isea+1;        
        seaice_file2{isea} = '2007/nt_20070901_f13_v01_s.bin'; isea=isea+1;                
        seaice_file2{isea} = '2007/nt_20071001_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20071001_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20071101_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20071201_f13_v01_s.bin'; isea=isea+1;
        seaice_file2{isea} = '2007/nt_20071231_f13_v01_s.bin'; isea=isea+1;        


    
    case 'Antarctic fortnightly'
        seaice_dir = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/daily/';
        for isea2=1:26  %length(files)
            isea = (isea2-1)*14 + 1;
            file_str = datestr( datenum('01-Jan-2007')+isea-1 ,'mmdd');
            seaice_file2{isea2} = ['2007/nt_2007' file_str '_f13_v01_s.bin'];
        end
            seaice_file2{isea2} = ['2007/nt_20071231_f13_v01_s.bin']; %also load in 31st Dec
            
    case {'Antarctic daily','Arctic daily'}
        
        switch pole
            case 'Antarctic daily'
                seaice_dir = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/daily/';
                pole_str='s';
            case 'Arctic daily'
                seaice_dir = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/daily/';
                pole_str='n';
        end

        start_date = datenum('15-Nov-2005');
        end_date = datenum('15-Jan-2008');
        
        start_date = datenum('15-May-2008');
        end_date = datenum('15-Sep-2008');    
        
        start_date = datenum('01-Jan-2000');
        end_date = datenum('31-Dec-2015');       
        
        period = 14;
        period = 1;
        ndays = end_date - start_date + 1;
        nfiles = floor(ndays/period);
        
        seaice_array_datenum = NaN*ones([1 nfiles]);
        
        
        for isea2=1:nfiles  %length(files)
            isea = (isea2-1)*period + 1;
            
            file_date = start_date +isea-1;

            file_str = datestr( file_date ,'mmdd');
            year_str = datestr( file_date ,'yyyy');
            if str2num(year_str) <= 2006
                    fstr = 'f13';
            else
                    fstr = 'f17';
            end
            if str2num(year_str) <= 2010
                    vstr = 'v01';
            else
                    vstr = 'v1.1';
            end
            
            
            seaice_file2{isea2} = [year_str '/nt_' year_str file_str '_' fstr '_' vstr '_' pole_str '.bin'];
            seaice_array_datenum(isea2) = file_date;
        end
%        seaice_file2{isea2} = ['2007/nt_20071231_f13_v01_s.bin']; %also load in 31st Dec

end

for isea=1:length(seaice_file2)
    fprintf(1,'\nLoading day %i out of %i',isea,length(seaice_file2));
    seaice_file = seaice_file2{isea};
    sea_str = seaice_file2{isea}(13:16);

    ioverride_seaice=1;
    read_sea_ice_nimbus
    close(gcf);

    eval(['seaice_data_' sea_str '=seaice_data;']);

    %create the array first time through
    if isea==1
        seaice_array = NaN*ones([size(seaice_data) length(seaice_file2)]);
    end
    
    seaice_array(:,:,isea) = seaice_data;
    
    
end

fprintf(1,'\nDone read_sea_ice_multi\n');


if isave_seaice==1
   seaice_savefile = [seaice_dir '/saved_seaice_' tag '_' ...
       datestr(now,30) '.mat']; 
 
   save(seaice_savefile,'seaice_lon','seaice_lat','seaice_array','seaice_array_datenum','-V7.3');
 
end
    
    
   