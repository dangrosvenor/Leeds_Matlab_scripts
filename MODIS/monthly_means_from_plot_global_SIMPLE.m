try
    
clear season_mean_ALL season_Ndatap_ALL season_std_ALL

if ~exist('ioverride_monthly_options') | ioverride_monthly_options==0

    multi_case = 'CTP';
    multi_case = 'Monthly';
%    multi_case = 'Fortnightly';
    %multi_case = 'Bimonthly';
    %multi_case = 'Annual';
%    multi_case = 'Seasonal';
    %multi_case = 'Regional';
    %multi_case = '20th of each month, vs lat';

   
    
    isave=0;   
    
    data_select='specific_modis_data';
    years_required_for_mean = unique(modisyear_timeseries3);
    screen_type ='none';
   
%    proj_type = 'global oval';
    ifilter_ndays=0;       

end

clear season_mean season_Ndatap season_std season_datatype season_vals_timelabs season_LON_str season_LAT_str season_zonal
clear season_mean_vs_LAT_2D  season_Ndatap_vs_LAT_2D season_max_vs_LAT_2D season_min_vs_LAT_2D season_std_vs_LAT_2D
clear season_timestd_vs_LAT_2D season_meanoverall_vs_LAT_2D season_Noverall_vs_LAT_2D season_stdoverall_vs_LAT_2D


%for '20th of each month, vs lat' will do global plots - so need to switch
%off restrict_domain
if ~exist('ioverride_location_selection2') | ioverride_location_selection2==0

    %    LAT_vals = [-22 -18; -22 -18; -22 -18; -22 -18]; LON_vals = [-76.5 -72; -81 -76.5; -85.5 -81; -90 -85.5];
    LAT_vals = [-22 -18; -22 -18; -22 -18; -20 -10]; LON_vals = [-75 -70; -80 -75; -86 -80; -90 -80];
    LAT_vals = [-22 -18; -22 -18; -22 -18; -20 -10]; LON_vals = [-75 -70; -80 -75; -86 -80; -90 -80];
    LAT_vals = [-22.74 -18; -22.74 -18; -22.74 -18; -20 -10]; LON_val = [-76.25 -71.25; -81.25 -76.25; -87.25 -81.25; -91.25 -81.25];


    %    LAT_vals = [-60 -40;]; LON_vals = [-179.5 179.5;]; %All lons
    %    LAT_vals = [-60 -40;]; LON_vals = [-179.5 179.5;]; %All lons
    LAT_vals = [-60 -40;]; LON_vals = [50 100;]; %JonnyW box
%    LAT_vals = [-60 -45;]; LON_vals = [50 150;]; %JonnyW box
    %    LAT_vals = [-60 -40;]; LON_vals = [-70 50;];
    %    LAT_vals = [-60 -40;]; LON_vals = [-179.5 -100;];

    %    thresh_LAT = [0 60]; thresh_LON = [90 180]; %China
    
    
    LAT_vals = [-55 -35;]; LON_vals = [-179.5 179.5;]; %Jane M./ Dan McCoy SO box , all lons    

    thresh_ndays=1;

else
    ioverride_location_selection=1;
end

thresh_LAT = LAT_vals; thresh_LON = LON_vals;


CTPs=[900:-100:200];


clear monthly_mean monthly_mean2 seasonal_mean seasonal_mean2 season_vals_save seasonal_NX seasonal_NY

switch multi_case
    case 'Monthly_old'

        for month_no=1:12
            
            % --- Have now created the days_of_month function to give
            % e.g. 32:60 for Feb --

            switch month_no
                case 1
                    days_required_for_mean = [1:31]; time_mean_str = 'Jan';
                case 2
                    days_required_for_mean = [32:60]; time_mean_str = 'Feb';
                case 3
                    days_required_for_mean = [61:91]; time_mean_str = 'Mar';
                case 4
                    days_required_for_mean = [92:121]; time_mean_str = 'Apr';
                case 5
                    days_required_for_mean = [122:152]; time_mean_str = 'May';
                case 6
                    days_required_for_mean = [153:182]; time_mean_str = 'Jun';
                case 7
                    days_required_for_mean = [183:213]; time_mean_str = 'Jul';
                case 8
                    days_required_for_mean = [214:244]; time_mean_str = 'Aug';
                case 9
                    days_required_for_mean = [245:274]; time_mean_str = 'Sep';
                case 10
                    days_required_for_mean = [275:305]; time_mean_str = 'Oct';
                case 11
                    days_required_for_mean = [306:335]; time_mean_str = 'Nov';
                case 12
                    days_required_for_mean = [336:366]; time_mean_str = 'Dec';
            end



            ioverride_time_selection=1;
            eval(plot_script);

            switch plot_script
                case 'plot_global_maps'
                    monthly_mean(month_no)=Pmean;
                case 'plotTimeHeightVap3'
                    monthly_mean(month_no)= X_mean_overall;
                    monthly_mean2(month_no)= Y_mean_overall;
            end


            saveas_ps_fig_emf(gcf,[savename]);

        end

        seasons(1)=mean(monthly_mean([12 1:2]));  sea_month(1)=1;
        seasons(2)=mean(monthly_mean([3:5]));  sea_month(2)=4;
        seasons(3)=mean(monthly_mean([6:8]));  sea_month(3)=7;
        seasons(4)=mean(monthly_mean([9:11]));  sea_month(4)=10;



    case {'Seasonal','Regional','Monthly','Bimonthly','Fortnightly','Annual','20th of each month, vs lat'}
        switch multi_case
            case 'Seasonal'
                Npdfs = 4;
            case 'Regional'
                Npdfs = size(LAT_vals,1);
            case {'Monthly','20th of each month, vs lat'}
                Npdfs=12;
            case {'Fortnightly'}
                Npdfs=26;                
            case 'Bimonthly'
                Npdfs=6;
            case 'Annual'
                Npdfs=1;
        end

        for iseason=1:Npdfs


            switch multi_case
                case {'Monthly','20th of each month, vs lat','Fortnightly'}

                    ioverride_time_selection=1;
                    xlab_str = 'Month';

                    %for CALIPSO
                    months_required_for_mean = iseason;


                    switch multi_case
                        case '20th of each month, vs lat'
                            LAT_str = ['Global'];
                            LON_str = [''];

                            %N.B. - based these on 2007, so may be different for
                            %other years
                            switch iseason
                                case 1
                                    %changed this to 19 since 20th Jan 2007
                                    %is missing
                                    days_required_for_mean = [19]; time_mean_str = '20-Jan';
                                case 2
                                    days_required_for_mean = [51]; time_mean_str = '20-Feb';
                                case 3
                                    days_required_for_mean = [79]; time_mean_str = '20-Mar';
                                case 4
                                    days_required_for_mean = [110]; time_mean_str = '20-Apr';
                                case 5
                                    days_required_for_mean = [140]; time_mean_str = '20-May';
                                case 6
                                    days_required_for_mean = [171]; time_mean_str = '20-Jun';
                                case 7
                                    days_required_for_mean = [201]; time_mean_str = '20-Jul';
                                case 8
                                    days_required_for_mean = [232]; time_mean_str = '20-Aug';
                                case 9
                                    days_required_for_mean = [263]; time_mean_str = '20-Sep';
                                case 10
                                    days_required_for_mean = [293]; time_mean_str = '20-Oct';
                                case 11
                                    days_required_for_mean = [324]; time_mean_str = '20-Nov';
                                case 12
                                    days_required_for_mean = [354]; time_mean_str = '20-Dec';
                            end



                        case 'Monthly'
                            LAT_str = [num2str(thresh_LAT(1)) ' to ' num2str(thresh_LAT(2))];
                            LON_str = [num2str(thresh_LON(1)) ' to ' num2str(thresh_LON(2))];



                            % --- Have now created the days_of_month function to give
                            % e.g. 32:60 for Feb --
                            %for MODIS
                            switch iseason
                                case 1
                                    days_required_for_mean = [1:31]; time_mean_str = 'Jan';
                                case 2
                                    days_required_for_mean = [32:60]; time_mean_str = 'Feb';
                                case 3
                                    days_required_for_mean = [61:91]; time_mean_str = 'Mar';
                                case 4
                                    days_required_for_mean = [92:121]; time_mean_str = 'Apr';
                                case 5
                                    days_required_for_mean = [122:152]; time_mean_str = 'May';
                                case 6
                                    days_required_for_mean = [153:182]; time_mean_str = 'Jun';
                                case 7
                                    days_required_for_mean = [183:213]; time_mean_str = 'Jul';
                                case 8
                                    days_required_for_mean = [214:244]; time_mean_str = 'Aug';
                                case 9
                                    days_required_for_mean = [245:274]; time_mean_str = 'Sep';
                                case 10
                                    days_required_for_mean = [275:305]; time_mean_str = 'Oct';
                                case 11
                                    days_required_for_mean = [306:335]; time_mean_str = 'Nov';
                                case 12
                                    days_required_for_mean = [336:366]; time_mean_str = 'Dec';
                            end
                            
                        case 'Fortnightly'
                            LAT_str = [num2str(thresh_LAT(1)) ' to ' num2str(thresh_LAT(2))];
                            LON_str = [num2str(thresh_LON(1)) ' to ' num2str(thresh_LON(2))];



                            %for MODIS

                                %work out the day based on iseason (number
                                %of fortnights into the year)
                                    day_current = (iseason-1)*14 + 1;
                                    day_start_str = datestr(datenum('01-Jan')+day_current-1,'dd mmm');
                                    day_end_str = datestr(datenum('01-Jan')+day_current+13-1,'dd mmm');                                    
                                    days_required_for_mean = [day_current:day_current+13]; time_mean_str = day_start_str;                                  

                                    seaice_start_str = datestr(datenum(day_start_str),'mmdd');
                                    seaice_end_str = datestr(datenum(day_end_str)-1,'mmdd'); %make the end date that of the starth
                                    %of the next fortnight to save number
                                    %of files that need to be loaded
                                    
                                   %except for last fortnight of the year
                                   if iseason==26
                                       seaice_end_str = '1231';
                                   end

                            

                    end

                    pdflab = time_mean_str;

                case 'Annual'
                    ioverride_time_selection=1;

                    xlab_str = 'Month';

                    switch iseason
                        case 1
                            days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';
                    end

                    pdflab = time_mean_str;

                case 'Bimonthly'
                    ioverride_time_selection=1;

                    xlab_str = 'Month';

                    switch iseason
                        case 1
                            days_required_for_mean = [1:60]; time_mean_str = 'JanFeb';
                        case 2
                            days_required_for_mean = [61:121]; time_mean_str = 'MarApr';
                        case 3
                            days_required_for_mean = [122:182]; time_mean_str = 'MayJun';
                        case 4
                            days_required_for_mean = [183:244]; time_mean_str = 'JulAug';
                        case 5
                            days_required_for_mean = [245:305]; time_mean_str = 'SepOct';
                        case 6
                            days_required_for_mean = [306:366]; time_mean_str = 'NovDec';
                    end

                    pdflab = time_mean_str;

                case 'Seasonal'
                    xlab_str = 'Season';
                    ioverride_time_selection=1;

                    switch iseason
                        case 1
                            %straight DJF
                            days_required_for_mean = [336:366 1:60]; time_mean_str = 'DJF';
                        case 2
                            %    %straight MAM
                            days_required_for_mean = [61:152]; time_mean_str = 'MAM';
                        case 3
                            %    %straight JJA
                            days_required_for_mean = [153:244]; time_mean_str = 'JJA';
                        case 4
                            %    %straight SON
                            days_required_for_mean = [245:335]; time_mean_str = 'SON';
                    end

                    pdflab = time_mean_str;

                case 'Regional'
                    xlab_str = 'Region';
                    ioverride_location_selection=1;

                    LAT_val = LAT_vals(iseason,:);
                    LON_val = LON_vals(iseason,:);



            end


            x_lab_seasonal(iseason).lab=xlab_str;

            switch plot_script
                case 'plot_global_maps'
                    ioverride_plotglobal_loc=1;
                    ioverride_plotglobal_thresh=1;  %comment out if want to use screenings set in plot_global
                    %time override should aready be set (ioverride_time_selection)
                    ioverride_years_time_screen=1; %required to specify the different years
                    inew_cticks=1;  %colorbar is of the non-linear type
               %Run the plot global script
                    plot_global_maps
                    if exist('isave') & isave==1
                        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',save_notes_filepath);
                        %                    saveas_ps_fig_emf(gcf,savename,'',0);
                    end
                    
                    %Store the monthly data for all the locations.
                    season_mean_ALL(:,:,iseason) = P_save;
                    season_Ndatap_ALL(:,:,iseason) = Npoints;
                    season_std_ALL(:,:,iseason) = P_std_dev;
                    
                    if length(Pmean)>0
                        season_mean(iseason) = Pmean;
                        season_Ndatap(iseason) = Pnpoints;
                        season_Ndatap_tot(iseason) = Pntotal_avail;
                        season_std(iseason) = Pstd;
                        [mean_sig2,Nspatial] = meanNoNan(Pstd(:).^2,1);
                        season_timestd(iseason) = sqrt(mean_sig2./Nspatial);
                        
                        %average zonally (over lons)
%                        season_zonal(iseason,:) = meanNoNan(P_save,2);

                        season_datatype{iseason} = modis_data_plot;
                        season_vals_timelabs{iseason} = time_mean_str;

                        if exist('LON_str')
                            season_LON_str{iseason} = LON_str;
                            season_LAT_str{iseason} = LAT_str;
                        end
                        
                        [tmp_mean,tmp_N,tmp_std] = meanNoNan(P_save,2);
                        season_mean_vs_LAT_2D(1:length(MLAT),iseason) = tmp_mean;
                        season_Ndatap_vs_LAT_2D(1:length(MLAT),iseason) = tmp_N;

                        switch multi_case
                            case '20th of each month, vs lat'



                                season_max_vs_LAT_2D(1:length(MLAT),iseason) = max(P,[],2);
                                season_min_vs_LAT_2D(1:length(MLAT),iseason) = min(P,[],2);

                                %spatial std for all of the time means
                                season_std_vs_LAT_2D(1:length(MLAT),iseason) = tmp_std;
                                [tmp_mean_sig2,tmp_Nspatial] = meanNoNan(tmp_std.^2,1);
                                season_timestd_vs_LAT_2D(1:length(MLAT),iseason) = sqrt(tmp_mean_sig2./tmp_Nspatial);
                                if exist('dat_modis2')
                                    [tmp_meanoverall,tmp_Noverall(iseason),tmp_stdoverall(iseason)] = meanNoNan(dat_modis2(:),1);
                                    season_meanoverall_vs_LAT_2D(1:length(MLAT),iseason) = tmp_meanoverall;
                                    season_Noverall_vs_LAT_2D(1:length(MLAT),iseason) = tmp_Noverall;
                                    season_stdoverall_vs_LAT_2D(1:length(MLAT),iseason) = tmp_stdoverall;
                                end

                        end

                        if exist('dat_modis2')
                            [tmp_meanoverall,tmp_Noverall(iseason),tmp_stdoverall(iseason)] = meanNoNan(dat_modis2(:),1);
                            season_meanoverall(iseason) = tmp_meanoverall;
                            season_Noverall(iseason) = tmp_Noverall(iseason);
                            season_stdoverall_(iseason) = tmp_stdoverall(iseason);
                        end

                    else
                        season_Ndatap(iseason) = 0;
                        season_mean(iseason) = NaN;
                        season_Ndatap_vs_LAT_2D(1:length(MLAT),iseason) = 0;
                        season_mean_vs_LAT_2D(1:length(MLAT),iseason) = NaN;

                    end



                case 'plotTimeHeightVap3'
                    am3_time_of_day_select=0;
                    Xbins=[0:25:1e3]; ichoose_Xbins=1;
                    ichoose_Xbins=0; nXpdf=17;
                    ioverride_pdf=1;
                    thresh_SZA=[45 55];
                    thresh_SZA=[45 55];
                    thresh_SZA=[0 90];
                    %thresh_SZA=[0 65];
                    %thresh_SZA=[50 55];
                    %thresh_SZA=[75 82];
                    %thresh_SZA=[46 56];
                    %thresh_SZA=[75 90];
                    thresh_CF=0.8;
                    thresh_CF=[0.8 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
                    %thresh_CF=[0.99 1.000001];  %cf screening is >thresh_CF(1) & <=thresh_CF(2)
                    %thresh_CF=[0.01 0.8];
                    %thresh_CF=[0.0 1.000001];
                    thresh_NP=10;
                    thresh_NP=50;
                    thresh_sensZA=45;
                    %thresh_sensZA=50;
                    thresh_sensZA=[0 41.4];
                    thresh_sensZA=[0 90];
                    %thresh_sensZA=[55 65];
                    %thresh_sensZA=[30 70];
                    %thresh_sensZA=58;

                    thresh_CTH = [-20 20];
                    thresh_CTH = [-20 2];

                    thresh_AZ = [50 130];

                    thresh_relAZ = [0 180];

                    thresh_CTT = [273-5 273+100];
                    %thresh_CTT = [273 273+100];
                    %thresh_CTT = [273 273+100];
                    thresh_CTT = [273-100 273+100];
                    thresh_CTP = 800; %Cloud top pressure (hPa)

                    %    thresh_sensZA=80;
                    thresh_maxSZA=200;

                    thresh_Nd_per_error = 100;
                    thresh_Reff_per_error = 50;
                    thresh_Reff_abs_error = 4;
                    thresh_Reff = 30;

                    thresh_stdW = [0 1e9];




                    %screen_type='NP + CF + MAX sensZA';
                    %                                    screen_type='NP + MAX sensZA';
                    screen_type='NP + CF + MEAN sensZA';
                    screen_type='NP + CF + MEAN sensZA + MEAN relAZ';
                    screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA';
                    screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT';
                    screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + mean_CTT + noice';
                    screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT + noice';
                    screen_type='NP + CF + MEAN sensZA + MEAN relAZ + stdLWP + MEAN solarZA + min_CTT';
                    screen_type ='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT';
                    %screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT + min_tau + mean_CTH';
                    %screen_type='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';




                    plotTimeHeightVap3
                    close(gcf);
                    if exist('ioverride_location_selection2') & ioverride_location_selection2==1
                        ioverride_location_selection=1;
                    end


                    %                    saveas_ps_fig_emf(gcf,[savename]);
                    seasonal_meanX(iseason)= X_mean_overall;
                    seasonal_meanY(iseason)= Y_mean_overall;
                    seasonal_NX(iseason) = sum(NX_vals);
                    seasonal_NY(iseason) = sum(NY_vals);

                    switch multi_case
                        case 'Regional'
                            pdflab =['LAT ' LAT_str ' LON ' LON_str];
                    end

                    time_highlight_path=[];

                    iytick_relabel=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
                    y_axis_type=''; %default
                    x_axis_type='';
                    i_set_dateticks=0;
                    iadd_nums_above=0;

                    xlims=0;
                    fsize=14;

                    idatetick=0; %flag to put times as HH:MM instead of decimal time

                    noplot=0;



                    graph=96;  %X_mean
                    man_choose_water_graph=1;
                    waterVapourMay2005
                    close(gcf);
                    %                    saveas_ps_fig_emf(gcf,[savename]);
                    season_vals_save(iseason).Yvals = ydat(1).y;
                    season_vals_save(iseason).Xmean = xdat(1).x;
                    season_vals_save(iseason).ylabelstr = ylabelstr;
                    season_vals_save(iseason).xlabelstr = xlabelstr;

                    graph=99; %
                    man_choose_water_graph=1;
                    waterVapourMay2005
                    close(gcf);
                    %                    saveas_ps_fig_emf(gcf,[savename]);
                    season_vals_save(iseason).Ndatap = xdat(1).x;

                    graph=977; %
                    man_choose_water_graph=1;
                    waterVapourMay2005
                    close(gcf);
                    %                    saveas_ps_fig_emf(gcf,[savename]);
                    season_vals_save(iseason).NdPDF = ydat(1).y;
                    season_vals_save(iseason).NdPDF_xbins = xdat(1).x;
                    season_vals_save(iseason).NdPDF_xlab = xlabelstr;
                    season_vals_save(iseason).pdflab = pdflab;
                    season_vals_save(iseason).pdf_ylab = ylab;

            end

        end


    case 'CTP'

        for ictp=1:length(CTPs)
            thresh_CTP = CTPs(ictp);
            eval(plot_script);
            saveas_ps_fig_emf(gcf,[savename]);
            ctp_mean(ictp)=Pmean;
        end

end




 clear ioverride_time_selection ioverride_years_time_screen ioverride_location_selection2 
catch time_error
    clear ioverride_time_selection ioverride_years_time_screen ioverride_location_selection2 
    rethrow(time_error);
end