%read monthly IPSL GOCCPv2.1 CALIPSO (no CloudSat included) data - see
%Chepfer, JGR, 2010
% N.B. - have saved data in mat files, so don't
%need to run the more lengthy 'process' action each time
%Be aware that daytime means daylight in this dataset. So e.g. at high lats
%(e.g. Arctic summer) there will be no nighttime data

%Makes data arrays of size [nyears*12 nalts nlat nlon]

%
calipso_daynight_label = 'nighttime';
calipso_daynight_label = 'daytime';


%whether we want the CF from CFADS or from the 3D CF product
data_type_cf_cal = 'CFAD';
data_type_cf_cal = 'CF3D';


switch data_type_cf_cal
    case 'CFAD'
        %for saved data
        savedir_data = '/home/disk/eos8/d.grosvenor/CPT/CFADS/processed_CFADs_saved_matlab_data/';
    case 'CF3D'
        savedir_data = '/home/disk/eos8/d.grosvenor/CPT/CFADS/processed_CF3D_saved_matlab_data/';       
end

%for plots
savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';



isavemem=1;

sr_thresh = 7;

action_flag = 'process'; %read in data, process to calculate the CF using a certain sr
%threshold and then save
action_flag = 'load'; %Load in previously processed data


switch action_flag

    case 'load'
        switch data_type_cf_cal
            case 'CFAD'
                switch calipso_daynight_label
                    case 'nighttime'
                        savefilename2='CFAD_CALIPSO_cf_thresh_1.2_nighttime.mat';
                        %        savefilename2='CFAD_CALIPSO_cf_thresh_3_nighttime.mat';
                        savefilename2='CFAD_CALIPSO_cf_thresh_7_nighttime.mat';
                    case 'daytime'
                        savefilename2='CFAD_CALIPSO_cf_thresh_7_daytime.mat';
                end

            case 'CF3D'
                switch calipso_daynight_label
                    case 'nighttime'
                        savefilename2='CF3D_CALIPSO_cf_nighttime.mat';
                    case 'daytime'
                        savefilename2='CF3D_CALIPSO_cf_daytime.mat';
                end
                                          
        end
        
        savefilename=[savedir_data savefilename2];

        load(savefilename);


    case 'process'

        clear field_names
        switch data_type_cf_cal
            case 'CFAD'
                switch calipso_daynight_label
                    case 'daytime';
                        nc_dir='/home/disk/eos8/d.grosvenor/CPT/CFADS/CALIPSO/daytime/';
                        ifield=1;
                        field_names{ifield}='cfadLidarsr532'; ifield=ifield+1; %
                        %        field_names{ifield}='cfad_lidarsr532_Occ2'; ifield=ifield+1; %
                    case 'nighttime'
                        nc_dir='/home/disk/eos8/d.grosvenor/CPT/CFADS/CALIPSO/nighttime/';
                        ifield=1;
                        field_names{ifield}='cfad_lidarsr532_Occ'; ifield=ifield+1; %
                        field_names{ifield}='cfad_lidarsr532_Occ2'; ifield=ifield+1; %
                end

            case 'CF3D'
                switch calipso_daynight_label
                    case 'daytime';
                        nc_dir='/home/disk/eos8/d.grosvenor/CPT/cf_3d/';
                        daynight_dirstr = 'day';
                        ifield=1;
                        field_names{ifield}='clcalipso'; ifield=ifield+1; %
                    case 'nighttime'
                        nc_dir='/home/disk/eos8/d.grosvenor/CPT/cf_3d/';
                        daynight_dirstr = 'night';
                        ifield=1;
                        field_names{ifield}='clcalipso'; ifield=ifield+1; %
                end

        end

        gcm_time_of_day_select=0;








        years_requested = [2007:2010]; %will leave out 2006 as it is only a partial year
        %years_requested = [2007]; %will leave out 2006 as it is only a
        %partial year (starts in June since the satellite only went up shortly before that)
        years_calipso_str='';
        for i=1:length(years_requested)
            years_calipso_str=[years_calipso_str ' ' num2str(years_requested(i))];
        end



        %% Search through all of the files in the directory to find the ones we
        %% want to process
        
        filesdir=[];
        switch data_type_cf_cal
            case 'CFAD'
                %files = dir([nc_dir field_names{ifield} '*']);
                files = dir([nc_dir '*']);
                for i=1:length(files)
                    filesdir(i).name = '';
                end
            case 'CF3D'
                ifiles=1;
                files=[];
               
                for i=1:length(years_requested)
                    filedirstr = [num2str(years_requested(i)) '/' daynight_dirstr '/'];
                    filenames = dir([nc_dir filedirstr '/*']);
                    nfiles = length(filenames);
                    files = [files; filenames];
                    for i=ifiles:ifiles+nfiles-1
                        filesdir(i).name = filedirstr;
                    end
                    ifiles = ifiles + nfiles;
                end

        end
        
        
        clear filenames_csat years_csat months_csat iyear_csat
        j=0;
        iyear=1;
        for i=1:length(files)

            if length(strfind(files(i).name,'.nc'))>0
                
                switch data_type_cf_cal
                    case 'CFAD'
                        switch calipso_daynight_label
                            case 'daytime';
                                year_of_file = str2num(files(i).name(55:58));
                            case 'nighttime'
                                year_of_file = str2num(files(i).name(14:17));
                        end

                    case 'CF3D'
                        switch calipso_daynight_label
                            case 'daytime';
                                year_of_file = str2num(files(i).name(22:25));
                            case 'nighttime'
                                year_of_file = str2num(files(i).name(22:25));
                        end

                end




                if length(find(year_of_file==years_requested))>0
                    j=j+1;
                    filenames_csat{j} = files(i).name;
                    dirnames_csat{j} = filesdir(i).name;    
                    
                    switch data_type_cf_cal
                        case 'CFAD'
                            switch calipso_daynight_label
                                case 'daytime';
                                    months_csat(j) = str2num(files(i).name(59:60));
                                case 'nighttime'
                                    months_csat(j) = str2num(files(i).name(15:16));
                            end

                        case 'CF3D'
                            switch calipso_daynight_label
                                case 'daytime';
                                    months_csat(j) = str2num(files(i).name(26:27));
                                case 'nighttime'
                                    months_csat(j) = str2num(files(i).name(26:27));
                            end

                    end



                    years_csat(j) = year_of_file;
                    %on first pass set to current year
                    if j==1
                        year_old = years_csat(j);
                    end

                    if years_csat(j) ~= year_old
                        iyear=iyear+1;
                    end
                    iyear_csat(j)=iyear;
                    year_old = years_csat(j);
                end

            end

        end

        %% Now process the files we identified

        %first just process the first file to check the size of the arrays

        nc = netcdf([nc_dir dirnames_csat{1} filenames_csat{1}]);
        dat = nc{field_names{1}}(:,:,:);
        %lat_csat = nc{'LAT'}(:);
        %lon_csat = nc{'LON'}(:);
        %nasc_csat = nc{'ASCDES'}(:);

        %for CF3D :-
        
%          longitude:lon_name = "Longitude" ;
%                 longitude:units = "degrees_east" ;
%                 longitude:axis = "X" ;
%         float latitude(latitude) ;
%                 latitude:lon_name = "Latitude" ;
%                 latitude:units = "degrees_north" ;
%                 latitude:axis = "Y" ;
%         float alt_mid(altitude) ;
%                 alt_mid:lon_name = "Middle of the altitude bin" ;
%                 alt_mid:units = "kilometer" ;
%                 alt_mid:positive = "up" ;
%                 alt_mid:axis = "Z" ;
%         float alt_bound(nv, altitude) ;
%                 alt_bound:lon_name = "Boundaries of the altitude bin" ;
%                 alt_bound:units = "kilometer" ;
%         float time(time) ;
%                 time:lon_name = "Time" ;
%                 time:units = "days since 2000-01-01 00:00:00" ;
%                 time:axis = "T" ;
%                 time:comment =

switch data_type_cf_cal
    case 'CF3D'

        switch calipso_daynight_label
            case {'daytime','nighttime'}
                gcm_lat = nc{'latitude'}(:);
                gcm_lon = nc{'longitude'}(:);

                cfad_alts_centre_CALIPSO = 1e3*nc{'alt_mid'}(:); %metres
                cfad_alts_edges_CALIPSO = 1e3*nc{'alt_bound'}(:);
%                cfad_alts_edges_CALIPSO = cfad_alts_edges_CALIPSO';
                cfad_alts_edges_CALIPSO = [cfad_alts_edges_CALIPSO(1,:) cfad_alts_edges_CALIPSO(2,end)];

                %                 cfad_cal_srs_centre_CALIPSO = nc{'srbox'}(:)/1e3; %seem to be 1e3 bigger than the nightime ones...?
                %                 cfad_cal_srs_edges_CALIPSO = nc{'srbox_bnds'}(:)/1e3;
                %                 cfad_cal_srs_edges_CALIPSO = cfad_cal_srs_edges_CALIPSO';
                %                 cfad_cal_srs_edges_CALIPSO = [cfad_cal_srs_edges_CALIPSO(1,:) cfad_cal_srs_edges_CALIPSO(2,end)];
                % cfad_cal_srs_centre_CALIPSO2 = nc{'srbox_mid2'}(:); %these are for the bleow  surface elev, rejected & noisy pixels
                % cfad_cal_srs_edges_CALIPSO2 = nc{'srbox_bound2'}(:);
                % cfad_cal_srs_edges_CALIPSO2 = [cfad_cal_srs_edges_CALIPSO2(1,:) cfad_cal_srs_edges_CALIPSO2(2,end)];

                %             case 'nighttime'
                %                 gcm_lat = nc{'latitude'}(:);
                %                 gcm_lon = nc{'longitude'}(:);
                %
                %                 cfad_alts_centre_CALIPSO = 1e3*nc{'alt_mid'}(:);
                %                 cfad_alts_edges_CALIPSO = 1e3*nc{'alt_bound'}(:);
                %                 cfad_alts_edges_CALIPSO = [cfad_alts_edges_CALIPSO(1,:) cfad_alts_edges_CALIPSO(2,end)];
                %                 cfad_cal_srs_centre_CALIPSO = nc{'srbox_mid'}(:);
                %                 cfad_cal_srs_edges_CALIPSO = nc{'srbox_bound'}(:);
                %                 cfad_cal_srs_edges_CALIPSO = [cfad_cal_srs_edges_CALIPSO(1,:) cfad_cal_srs_edges_CALIPSO(2,end)];
                %                 cfad_cal_srs_centre_CALIPSO2 = nc{'srbox_mid2'}(:); %these are for the bleow  surface elev, rejected & noisy pixels
                %                 cfad_cal_srs_edges_CALIPSO2 = nc{'srbox_bound2'}(:);
                %                 cfad_cal_srs_edges_CALIPSO2 = [cfad_cal_srs_edges_CALIPSO2(1,:) cfad_cal_srs_edges_CALIPSO2(2,end)];

        end

    case 'CFAD'
        switch calipso_daynight_label
            case 'daytime';
                gcm_lat = nc{'lat'}(:);
                gcm_lon = nc{'lon'}(:) - 180;

                cfad_alts_centre_CALIPSO = nc{'alt40'}(:);
                cfad_alts_edges_CALIPSO = nc{'alt40_bnds'}(:);
                cfad_alts_edges_CALIPSO =cfad_alts_edges_CALIPSO';
                cfad_alts_edges_CALIPSO = [cfad_alts_edges_CALIPSO(1,:) cfad_alts_edges_CALIPSO(2,end)];

                cfad_cal_srs_centre_CALIPSO = nc{'srbox'}(:)/1e3; %seem to be 1e3 bigger than the nightime ones...?
                cfad_cal_srs_edges_CALIPSO = nc{'srbox_bnds'}(:)/1e3;
                cfad_cal_srs_edges_CALIPSO = cfad_cal_srs_edges_CALIPSO';
                cfad_cal_srs_edges_CALIPSO = [cfad_cal_srs_edges_CALIPSO(1,:) cfad_cal_srs_edges_CALIPSO(2,end)];
                % cfad_cal_srs_centre_CALIPSO2 = nc{'srbox_mid2'}(:); %these are for the bleow  surface elev, rejected & noisy pixels
                % cfad_cal_srs_edges_CALIPSO2 = nc{'srbox_bound2'}(:);
                % cfad_cal_srs_edges_CALIPSO2 = [cfad_cal_srs_edges_CALIPSO2(1,:) cfad_cal_srs_edges_CALIPSO2(2,end)];

            case 'nighttime'
                gcm_lat = nc{'latitude'}(:);
                gcm_lon = nc{'longitude'}(:);

                cfad_alts_centre_CALIPSO = 1e3*nc{'alt_mid'}(:);
                cfad_alts_edges_CALIPSO = 1e3*nc{'alt_bound'}(:);
                cfad_alts_edges_CALIPSO = [cfad_alts_edges_CALIPSO(1,:) cfad_alts_edges_CALIPSO(2,end)];
                cfad_cal_srs_centre_CALIPSO = nc{'srbox_mid'}(:);
                cfad_cal_srs_edges_CALIPSO = nc{'srbox_bound'}(:);
                cfad_cal_srs_edges_CALIPSO = [cfad_cal_srs_edges_CALIPSO(1,:) cfad_cal_srs_edges_CALIPSO(2,end)];
                cfad_cal_srs_centre_CALIPSO2 = nc{'srbox_mid2'}(:); %these are for the bleow  surface elev, rejected & noisy pixels
                cfad_cal_srs_edges_CALIPSO2 = nc{'srbox_bound2'}(:);
                cfad_cal_srs_edges_CALIPSO2 = [cfad_cal_srs_edges_CALIPSO2(1,:) cfad_cal_srs_edges_CALIPSO2(2,end)];

        end



end


        dlat=mean(diff(gcm_lat));
        dlon=mean(diff(gcm_lon));

        gcm_slat = [gcm_lat-dlat/2; gcm_lat(end)+dlat/2];
        gcm_slon = [gcm_lon-dlon/2; gcm_lon(end)+dlon/2];


        nlat=length(gcm_lat);
        nlon=length(gcm_lon);


        nalts=length(cfad_alts_centre_CALIPSO);
        
        
        years_csat_unique = unique(years_csat);

        %nyears=length(years_csat_unique);
        nyears=length(years_requested);

        switch data_type_cf_cal
            case 'CFAD'

                nsrs=length(cfad_cal_srs_centre_CALIPSO);
                if exist('cfad_cal_srs_centre_CALIPSO2')
                    nsrs2=length(cfad_cal_srs_centre_CALIPSO2);
                end

                varname_cfad = ['cfad_lidarsr532_monthly_' calipso_daynight_label];

                clear field_name_out

                if isavemem==0

                    for ifield=1:length(field_names)
                        switch field_names{ifield}
                            case {'cfad_lidarsr532_Occ','cfadLidarsr532'}
                                eval_str = [varname_cfad '=NaN*ones([nyears*12 nsrs nalts nlat nlon]);'];
                                field_name_out{ifield} = varname_cfad;
                            case 'cfad_lidarsr532_Occ2'
                                eval_str = ['cfad_lidarsr532_monthly2_' calipso_daynight_label '=NaN*ones([nyears*12 nsrs2 nalts nlat nlon]);'];
                                field_name_out{ifield} = ['cfad_lidarsr532_monthly2_' calipso_daynight_label];
                        end
                        eval(eval_str);
                    end

                end

                eval_str = ['Ntot=NaN*ones([nyears*12 nalts nlat nlon]);'];
                eval(eval_str);

                istart=find(cfad_cal_srs_edges_CALIPSO>=sr_thresh);


        end

        eval_str = ['cf_CFAD_sr_CALIPSO=NaN*ones([nyears*12 nalts nlat nlon]);'];
        eval(eval_str);

        
        clear month_calipso_cf year_calipso_cf

        for j=1:length(filenames_csat)

            fprintf(1,'\nProcessing file %d out of %d',j,length(filenames_csat));

            ind = months_csat(j)+(years_csat(j)-years_csat(1))*12;
            %        low_CF_calipso = read_calipso_cmor(nc_dir,nc_inst_file,'cllcalipso',ilat,ilon,'(:,ilat,ilon)');

            for ifield=1:length(field_names)
                field_name = field_names{ifield};
                switch data_type_cf_cal
                    case 'CFAD'
                        eval(['dat = read_calipso_ipsl_night(nc_dir,[dirnames_csat{j} filenames_csat{j}],''' field_name ''',0,0,''(:,:,:,:)'',calipso_daynight_label);']);
                    case 'CF3D'
                        eval(['dat = read_calipso_ipsl_night(nc_dir,[dirnames_csat{j} filenames_csat{j}],''' field_name ''',0,0,''(:,:,:)'',''nighttime'');']); %use the nighttime flag
                        %for both day and night since the CF3D data is like
                        %the CFAD nighttime (fill value of -9999, etc.)
                end
                if isavemem==0
                    eval([field_name_out{ifield} '(ind,:,:,:,:) = 100*dat(:,:,:,:);']); %convert to % to be consistent with the other CALIPSO files
                end
                switch field_names{ifield}
                    case {'cfad_lidarsr532_Occ','cfadLidarsr532'}
                        Ntot(ind,:,:,:) = meanNoNan(dat,1,'sum');
                        %basically am separating the sr bins into cloudy
                        %and clear for each location, height and time. Cloud Fraction is the 
                        %fraction of obs that were cloudy. So is
                        %essentially a cloud frequency taken 
                        %over time (Ncloudy_obs/Ntotal_obs)
                        [cf_CFAD_sr_CALIPSO2,N2] = meanNoNan(dat(istart(1):end,:,:,:),1,'sum');
                        cf_CFAD_sr_CALIPSO(ind,:,:,:) = cf_CFAD_sr_CALIPSO2 ./ squeeze(Ntot(ind,:,:,:));
                    case 'clcalipso'
                        cf_CFAD_sr_CALIPSO(ind,:,:,:) = dat;
                end

            end



            %        eval([lower(field_name) '(ind+1,:,:) = dat(2,:,:);']); %descending
            %        eval([lower(field_name) '_daily(j,:,:) = 0.5*(dat(1,:,:) + dat(2,:,:));']); %daily average
            %        asc_calipso_matt(ind) = 1;
            %        asc_calipso_matt(ind+1) = 2;

            month_calipso_cf(ind) = months_csat(j);
            year_calipso_cf(ind) = years_csat(j);

            %        month_calipso_matt_daily(j) = months_csat(j);
            %        year_calipso_matt_daily(j) = years_csat(j);


            %    end

        end








        Plat=gcm_slat;
        Plon=gcm_slon;


        i180=find(Plon>180);
        Plon2=Plon;
        Plon2(i180)=Plon2(i180)-360;

        %Plon=[Plon2(1:end); Plon2(1)];
        Plon=Plon2;



        [gcm_Plon2D_edges_CALIPSO_monthly,gcm_Plat2D_edges_CALIPSO_monthly]=meshgrid(Plon,Plat);


        Plat=gcm_lat;
        Plon=gcm_lon;

        i180=find(Plon>180);
        Plon(i180)=Plon(i180)-360;

        [gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly]=meshgrid(Plon,Plat);





        sr_thresh_CALIPSO = sr_thresh;

        switch data_type_cf_cal
            case 'CFAD'
                savefilename2 = ['CFAD_CALIPSO_cf_thresh_' num2str(sr_thresh) '_' calipso_daynight_label '.mat']
                savefilename=[savedir_data savefilename2];
                save(savefilename,'cf_CFAD_sr_CALIPSO','Ntot','sr_thresh_CALIPSO','calipso_daynight_label','years_calipso_str','gcm_Plon2D_edges_CALIPSO_monthly','gcm_Plat2D_edges_CALIPSO_monthly','gcm_Plon2D_CALIPSO_monthly','gcm_Plat2D_CALIPSO_monthly','cfad_alts_edges_CALIPSO','cfad_alts_centre_CALIPSO','-V7.3');
            case 'CF3D'
                savefilename2 = ['CF3D_CALIPSO_cf_' calipso_daynight_label '.mat']
                savefilename=[savedir_data savefilename2];
                save(savefilename,'cf_CFAD_sr_CALIPSO','calipso_daynight_label','years_calipso_str','gcm_Plon2D_edges_CALIPSO_monthly','gcm_Plat2D_edges_CALIPSO_monthly','gcm_Plon2D_CALIPSO_monthly','gcm_Plat2D_CALIPSO_monthly','cfad_alts_edges_CALIPSO','cfad_alts_centre_CALIPSO','-V7.3');
        end


end

gcm_years_loaded_str = [years_calipso_str ' OBS ' calipso_daynight_label ' '];

daynum_timeseries3_CALIPSO = 1:size(cf_CFAD_sr_CALIPSO,1);

daynum_timeseries3_CALIPSO_monthly = 1:size(cf_CFAD_sr_CALIPSO,1);
gcm_time_UTC_CALIPSO_monthly = zeros([1 length(daynum_timeseries3_CALIPSO_monthly)]);

fprintf(1,'\n Done read calipso monthly\n');






