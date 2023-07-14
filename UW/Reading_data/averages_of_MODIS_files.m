%read and average multiple modis files
% The script to specify what variables to load depends on the Collection -
% search for read_MODIS_variables_01 below
% Can specify multiple dataset_modis_list values and it will loop through
% each one and make a .mat file in 
% savedir_var='~/modis_work/saved_data_L3/'  with filename:-
% savevarname = [savedir_var 'timeseries3_data_' dataset_modis '_' datestr(now,30)];  
% Can be used just to read in the data and make a .mat file for each day,
% as opposed to a the daily .hdf files that are origninally present that
% are harder to deal with.Can also just extract the key variables.
% But it also runs make_timeseries_MODIS that concatenates the daily arrays
% together.

collection2='C5'; 
collection2='C6.1'; 

idataset=1; clear dataset_modis_list;

switch collection2
    case 'C5'   
        collection='C5';
        
        %dataset_modis_list{idataset} = 'y2000_TERRA'; idataset=idataset+1;  %day
        %55-366 (launch year) %117-188 & 219-230 not available form LAADS
        %dataset_modis_list{idataset} = 'y2001_TERRA'; idataset=idataset+1; %219-230 not available from LAADS
        
        %dataset_modis_list{idataset} = 'y2002_TERRA';  idataset=idataset+1; %79-88 missing from LAADS web
        %but Rob has 87 & 88 in his dir??  199-203 were missing, but have
        %downloaded. 105 missing form LAADS web.
        
        dataset_modis_list{idataset} = 'y2002_AQUA'; idataset=idataset+1; %184-365 as Aqua launched 2002
        %186 was labelled with a .1 at the end-? 211-217 missing from LADS web. Downloaded 300, which was
        %missing
        
        %dataset_modis_list{idataset} = 'y2003_TERRA'; idataset=idataset+1; %351-358 missing form LAADS web - downloaded
        dataset_modis_list{idataset} = 'y2003_AQUA'; idataset=idataset+1;  %
        
        
        %dataset_modis_list{idataset} = 'y2004_TERRA';  idataset=idataset+1; %
        dataset_modis_list{idataset} = 'y2004_AQUA';  idataset=idataset+1; %
        
        
        %dataset_modis_list{idataset} = 'y2005_TERRA'; idataset=idataset+1;
        dataset_modis_list{idataset} = 'y2005_AQUA';  idataset=idataset+1; %full
        %files
        
        %dataset_modis_list{idataset} = 'y2006_TERRA'; idataset=idataset+1;
        %dataset_modis_list{idataset} = 'y2006_AQUA'; idataset=idataset+1;
        dataset_modis_list{idataset} = 'y2006_AQUA_updated'; idataset=idataset+1;
        
        %dataset_modis_list{idataset} = 'y2007_TERRA'; idataset=idataset+1;
        dataset_modis_list{idataset} = 'y2007_AQUA_51'; idataset=idataset+1; %all coll 5.1
        
        
        %dataset_modis_list{idataset} = 'y2008_TERRA'; idataset=idataset+1; %coll 5
        dataset_modis_list{idataset} = 'y2008_AQUA'; idataset=idataset+1; %coll 5    -
        
        %dataset_modis_list{idataset} = 'y2009_TERRA'; idataset=idataset+1; %Rob's coll 5
        dataset_modis_list{idataset} = 'y2009_AQUA'; idataset=idataset+1; %coll 5.1- downloaded, only have 1-99
        
        %dataset_modis_list{idataset} = 'y2010_TERRA'; idataset=idataset+1; %
        dataset_modis_list{idataset} = 'y2010_AQUA'; idataset=idataset+1; %
        
        %dataset_modis_list{idataset} = 'y2011_TERRA'; idataset=idataset+1; %
        dataset_modis_list{idataset} = 'y2011_AQUA'; idataset=idataset+1; %
        
        %dataset_modis_list{idataset} = 'y2012_TERRA'; idataset=idataset+1; %
        dataset_modis_list{idataset} = 'y2012_AQUA'; idataset=idataset+1; %
        
        %dataset_modis_list{idataset} = 'y2013_TERRA'; idataset=idataset+1; %
        dataset_modis_list{idataset} = 'y2013_AQUA'; idataset=idataset+1; %
        
        %dataset_modis_list{idataset} = 'y2014_TERRA'; idataset=idataset+1; %
        dataset_modis_list{idataset} = 'y2014_AQUA'; idataset=idataset+1; %
        
        
        
        %%added days 363-365 (but are collection 5.1, others are 5)
        
    case 'C6'
        
        %% Collection 6
        
        
        
        %dataset_modis_list{idataset} = 'y2005_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2006_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2007_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2008_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2009_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2010_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2011_AQUA_C6'; idataset=idataset+1;  %
        dataset_modis_list{idataset} = 'y2012_AQUA_C6'; idataset=idataset+1;  %
        %dataset_modis_list{idataset} = 'y2013_AQUA_C6'; idataset=idataset+1;  %
        %dataset_modis_list{idataset} = 'y2014_AQUA_C6'; idataset=idataset+1;  %
        
    case 'C6.1'
        %dataset_modis_list{idataset} = 'y2009_AQUA_C61_Monthly'; idataset=idataset+1; %E.g., monthly AOD
        %for year=2004:2022
        for year=2021:2021
            years_dataset{idataset}=year; dataset_modis_list{idataset} = 'various_years_AQUA_C61_Monthly'; idataset=idataset+1; %E.g., monthly AOD
        end
        %dataset_modis_list{idataset} = 'various_years_AQUA_C61_Monthly'; idataset=idataset+1; %E.g., monthly AOD
        %dataset_modis_list{idataset} = 'y2016_AQUA_C61'; idataset=idataset+1; %for Leighton
        %dataset_modis_list{idataset} = 'y2017_AQUA_C61'; idataset=idataset+1; %for Leighton
        
end


for idataset=1:length(dataset_modis_list)
    dataset_modis = dataset_modis_list{idataset};
    
    %Set default filedir - may be changed below
    filedir='/home/disk/eos10/robwood/';
    isave_MODIS_data=1;  %flag to say whether to save all of the data gathered or not

   
    clear modis_dir

    switch dataset_modis
        case 'y2000_TERRA'
            modis_dir='MODIS/D3/y2000/';   %C5 (not C5.1)
            %        files_to_read=[1 3:299];  %no. 300 is just a checksum and number 2 is a data list or something
        case 'y2001_TERRA'   %C5 (not C5.1)
            modis_dir='MODIS/D3/y2001/';    
        case 'y2002_AQUA' %only 184-365 as Aqua launched in 2002. Mixture of C5 and C5.1
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2002_a/';
        case 'y2002_TERRA'  %C5 (not C5.1)
            modis_dir='MODIS/D3/y2002/';
        case 'y2003_TERRA'  %C5 (not C5.1), except for a few days that were originally missing.
            modis_dir='MODIS/D3/y2003/';
        case 'y2003_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2003_a/';  %C5 (not C5.1)

            %         modis_dir='MODIS/D3/y2003_a/';
        case 'y2005_TERRA'
            modis_dir='MODIS/D3/y2005/';   %C5 (not C5.1)
            %        files_to_read=[1 3:299];  %no. 300 is just a checksum and number 2 is a data list or something
        case 'y2004_TERRA'
            %        modis_dir='MODIS/D3/y2004/';
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2004/';   %Collecction 5.1
        case 'y2004_AQUA'
            modis_dir='MODIS/D3/y2004_a/'; %C5 (not C5.1)
        case 'y2005_AQUA'
            modis_dir='MODIS/D3/y2005_a/'; %C5 (not C5.1)
        case 'y2006_TERRA'
            modis_dir='MODIS/D3/y2006/'; %C5 (not C5.1)
        case 'y2006_AQUA'
            modis_dir='MODIS/D3/y2006_a/';  %C5 (not C5.1), but with a few days 5.1
        case 'y2006_AQUA_updated'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2006_a/'; %C5 (not C5.1), but with a few days 5.1
        case 'y2007_TERRA'
            modis_dir='MODIS/D3/y2007/';
        case 'y2007_AQUA_51'  %this wasn't present so have downloaded it (collection 5.1)
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2007_a/';
        case 'y2008_TERRA' %collection 5 (Rob's)
            modis_dir='MODIS/D3/y2008/';
         case 'y2008_AQUA' %collection 5 (Rob's)
            modis_dir='MODIS/D3/y2008_a/';    
        case 'y2009_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2009_a/';  %Collecction 5.1
        case 'y2009_TERRA' %collection 5 (Rob's)
            modis_dir='MODIS/D3/y2009/';            
        case 'y2010_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2010_a/'; %Collecction 5.1
        case 'y2010_TERRA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2010/'; %Collecction 5.1
        case 'y2011_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2011_a/'; %Collecction 5.1
        case 'y2011_TERRA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2011/'; %Collecction 5.1
        case 'y2012_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2012_a/'; %Collecction 5.1
        case 'y2012_TERRA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2012/'; %Collecction 5.1            
        case 'y2013_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2013_a/'; %Collecction 5.1
        case 'y2013_TERRA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2013/'; %Collecction 5.1
        case 'y2014_AQUA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2014_a/'; %Collecction 5.1
        case 'y2014_TERRA'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/y2014/'; %Collecction 5.1  
            
%% Colleciton 6  
        case 'y2005_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2005_a/'; %Collecction 6 
            collection='C6';   
        case 'y2006_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2006_a/'; %Collecction 6 
            collection='C6';   
        case 'y2007_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2007_a/'; %Collecction 6 
            collection='C6';   
        case 'y2008_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2008_a/'; %Collecction 6 
            collection='C6';   
        case 'y2009_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2009_a/'; %Collecction 6 
            collection='C6';               
        case 'y2010_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2010_a/'; %Collecction 6 
            collection='C6';   
        case 'y2011_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2011_a/'; %Collecction 6 
            collection='C6';               
        case 'y2012_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2012_a/'; %Collecction 6 
            collection='C6';   
        case 'y2013_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2013_a/'; %Collecction 6 
            collection='C6';               
        case 'y2014_AQUA_C6'
            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='MODIS_L3_data/C6/y2014_a/'; %Collecction 6 
            collection='C6'; 
        case 'y2009_AQUA_C61_Monthly'
            filedir='/home/disk/eos15/d.grosvenor/';
            modis_dir='MODIS_L3/aqua_monthly/'; 
            collection='C61_Monthly_AOD'; 
            %Contains monthly data from day 60 of 2009 to day 91 of 2010
            %Day 60 of 2009 is 1st March 2009
            %Day 91 of 2010 is 1st April 2010  
            
        case 'various_years_AQUA_C61_Monthly'
            filedir='/home/disk/eos15/d.grosvenor/';
            %year='2003';
            year = num2str(years_dataset{idataset});
            modis_dir=['MODIS_L3/aqua_monthly/MYD08_M3/' year '/']; 
            collection='C61_Monthly_AOD'; 
            dataset_modis = ['y' year '_AQUA_C61_Monthly_AOD'];
            
        case 'y2016_AQUA_C61' %for Leighton - Dec 2016 and all of 2017.
            filedir='/home/disk/eos13/d.grosvenor/';
            modis_dir='MODIS_L3_data/C61/aqua/2016/'; %Collecction 6 
            collection='C6'; %For the vars in the files, etc.
            %Contains daily data from day 336 to 366 in 2016 (Dec 2016)
            
       case 'y2017_AQUA_C61' %for Leighton - Dec 2016 and all of 2017.
            filedir='/home/disk/eos13/d.grosvenor/';
            modis_dir='MODIS_L3_data/C61/aqua/2017/'; %Collecction 6 
            collection='C6'; %For the vars in the files, etc.
            %Contains daily data for all of 2017     
            
    end

    disp(dataset_modis);

    files_mod = dir([filedir modis_dir]);
    %files_mod(1:2)=[]; %these are the . and .. listings
    %these are removed below, along with any other files that don't end in .hdf

    ndir_files = length(files_mod);

    savedir_var='~/modis_work/saved_data_L3/'
%    savevarname = [savedir_var 'timeseries3_TProfiles_CTT_CTP_data_' dataset_modis '_' datestr(now,30)];
    savevarname = [savedir_var 'timeseries3_data_' dataset_modis '_' datestr(now,30)];    

    if exist('SD_id')
        status = hdfsd('end',SD_id);
    end

    %first of all go through all of the filenames to decide the ones that are
    %proper files (not checksums, etc.)
    clear files_to_read days_to_read
    imod_read=0;
    for imr=1:ndir_files
        %test to see if we want to process this file
        if length(files_mod(imr).name)>=4 & strcmp(files_mod(imr).name(end-3:end),'.hdf')==1
            imod_read=imod_read+1;
            files_to_read(imod_read)=imr;
            days_to_read(imod_read) = str2num(files_mod(imr).name(15:17));
        end

    end

    %files_to_read=[1:2];
    
    %days_to_read are the actual day numbers

    nMOD_av=length(files_to_read);

    for imr=1:nMOD_av
        fprintf(1,'\n%d of %d ',imr,nMOD_av);

        imod_file=files_to_read(imr);

        file_name_h5 = [modis_dir files_mod(imod_file).name];

        imodis_file_override=1;
        open_MODIS_file_01

        iaverage_modis=1;
        itimeseries_MODIS=1;
        
        switch collection
            case 'C5'
                read_MODIS_variables_01; %this will read the requested variables and average them
            case {'C6'}
                read_MODIS_variables_01_C6; %this will read the requested variables and average them
            case {'C61_Monthly_AOD'}
                read_MODIS_variables_01_C61_monthly_AOD
        end
        
        status = hdfsd('end',SD_id); %end access to the file - think otherwise things start to go wrong

    end


    if isave_MODIS_data==1
        save_MODIS_L3_data
    end
        

end

