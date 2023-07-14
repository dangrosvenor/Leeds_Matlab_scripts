%read and average multiple modis files
%read_MODIS_variables_01 specifies what variables to load
% Can specify multiple dataset_modis_list values and it will loop through
% each one and make a .mat file in 
% savedir_var='~/modis_work/saved_data_L3/'  with filename:-
% savevarname = [savedir_var 'timeseries3_data_' dataset_modis '_' datestr(now,30)];  
% Can be used just to read in the data and make a .mat file for each day,
% as opposed to a the daily .hdf files that are origninally present that
% are harder to deal with.Can also just extract the key variables.
% But it also runs make_timeseries_MODIS that concatenates teh daily arrays
% together.

isave_MODIS_data=1;

idataset=1; clear dataset_modis_list
    
%dataset_modis_list{idataset} = 'y2013_asc'; idataset=idataset+1;  %
dataset_modis_list{idataset} = 'y2013_desc'; idataset=idataset+1;  %



for idataset=1:length(dataset_modis_list)
    dataset_modis = dataset_modis_list{idataset};
    
    %Set default filedir of where to find the data to read in
    filedir='/home/disk/eos8/d.grosvenor/SMOS/SMOS_L3_CPDC/';
    isave_MODIS_data=1;  %flag to say whether to save all of the data gathered or not
    
    %Dir where to save .mat files
    savedir_var='/home/disk/eos8/d.grosvenor/SMOS/SMOS_L3_CPDC/saved_SMOS_L3/';
%    savevarname = [savedir_var 'timeseries3_TProfiles_CTT_CTP_data_' dataset_modis '_' datestr(now,30)];
    savevarname = [savedir_var 'timeseries3_data_' dataset_modis '_' datestr(now,30)];      

   
    clear modis_dir

    switch dataset_modis
        case 'y2013_asc'
            %            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='y2013/'; %Ascending
        case 'y2013_desc'
            %            filedir='/home/disk/eos8/d.grosvenor/';
            modis_dir='y2013_desc/'; %Ascending
    end

%    disp(dataset_modis);cd
% Need the above??

%    files_mod = dir([filedir modis_dir]);
    
    %files_mod(1:2)=[]; %these are the . and .. listings
    %these are removed below, along with any other files that don't end in .hdf
    
    files_mod=[];
    
    for i=1:12
        mon_str = num2str(i,'%02g');    
%        files_mod = eval(['!find ' filedir modis_dir ' -name "*.nc"']);    
        files_mod2 = dir([filedir modis_dir '/' mon_str '/']);
        for ifile=1:length(files_mod2)
            files_mod2(ifile).name = [mon_str '/' files_mod2(ifile).name];
        end
        files_mod = cat(1,files_mod,files_mod2);
    end

    ndir_files = length(files_mod);

  


    %first of all go through all of the filenames to decide the ones that are
    %proper files (not checksums, etc.)
    clear files_to_read days_to_read years_to_read
    imod_read=0;
    for imr=1:ndir_files
        %test to see if we want to process this file
        if length(files_mod(imr).name)>=4 & strcmp(files_mod(imr).name(end-2:end),'.nc')==1
            imod_read=imod_read+1;
            files_to_read(imod_read)=imr;
            year_file = str2num(files_mod(imr).name(23:26));            
            mon_file = str2num(files_mod(imr).name(27:28));
            day_file = str2num(files_mod(imr).name(29:30));
            day_of_year = day_of_year_from_date_func(datenum(year_file,mon_file,day_file));
            days_to_read(imod_read) = day_of_year;
            years_to_read(imod_read) = year_file;            
%            files_mod(imr).name
        end

    end

    %files_to_read=[1:2];
    
    %days_to_read are the actual day numbers

    nMOD_av=length(files_to_read);

    for imr=1:nMOD_av
        fprintf(1,'\n%d of %d ',imr,nMOD_av);

        modis_day_str = num2str(days_to_read(imr),'%03g');
        modis_year_str = num2str(years_to_read(imr));
        
        imod_file=files_to_read(imr);
        

        
        file_name_nc = [filedir modis_dir files_mod(imod_file).name];

        imodis_file_override=1;
%        open_MODIS_file_01

        nc = netcdf(file_name_nc);

        iaverage_modis=1;
        itimeseries_MODIS=1;
        read_SMOS_variables_01; %this will read the requested variables and average them

%        status = hdfsd('end',SD_id); %end access to the file - think otherwise things start to go wrong

        clear nc

    end
    
    [gcm_Plon2D_SMOS,gcm_Plat2D_SMOS]=meshgrid(MLON,MLAT);
    [gcm_Plat2D_edges_SMOS,gcm_Plon2D_edges_SMOS] = get_edges_lat_lon(gcm_Plat2D_SMOS,gcm_Plon2D_SMOS);
    gcm_str = 'SMOS';
    daynum_timeseries3_SMOS = days_to_read;
    gcm_time_UTC_SMOS = 0;
    modisyears_str = dataset_modis_list{1};


    if isave_MODIS_data==1
        save_MODIS_L3_data
    end
        

end

