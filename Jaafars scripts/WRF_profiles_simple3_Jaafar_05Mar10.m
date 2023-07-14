%This script is written to plot profile for several hour for WRF or
%RASS/SODAR (NOTE: YOU NEED TO LOAD BOTH FILES WRF AND RASS TO MAKE THE
%SCRIPT WORKING)
WRF_plot_simple_set_flags  %sets up flags and clears data (you always remeber to load the data first (running load data script)before running this script

%select case to run
graph=1; % (case 1)profiles at different locations (WRF)
%graph=2; % (case 2)our example x,y plot (RASS)

switch graph
    case 9999999
        % template to copy for new graph

        time=1; %time index

        ylab='Height (m)';
        xlab= 'WRF ice number concentration (L^{-1})';


        idat=0; %ignore this - it is a counter for the numebr of lines on the plot

        HGT=WRFUserARW(nca(1).nc,'Z',time,ilat(iloc),ilon(iloc)); %gets height data from WRF output
        %first profile
        idat=idat+1;
        xdat(idat).x=0.005*exp(0.304*(-tc)); %
        ydat(idat).y=HGT; %
        labs(idat).l='Profile 1';

        %second profile
        idat=idat+1;
        xdat(idat).x=0.005*exp(0.304*(-tc)); %x data
        ydat(idat).y=HGT; %y data
        labs(idat).l='Profile 2';

        %add more profiles here as required

        %change x-limits
        xlims=0;
        xlimits=1000*[0 0.025];

        zlims=1;
        zmin=100;
        zmax=3000;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane


        tstr=Times(time,:);
        iund=findstr('_',tstr);
        tstr(iund)=' ';
        titlenam = ['XXX for ' tstr];

        figname=titlenam;
        savename=figname;


        %%%%%%%%%%%%%%%%%%%% start of case 2 RASS Profile %%%%%%%%%%%%%%%%%%%%%
    case 2 %general plot of x against y (now is for RASS)

        %            time=12; %time index

        ylab='Height (m)';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% you chose one of these parameters to plot when case 2 is chosen %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %xlab='Wind Speed (m s^{-1})';
        %xlab='Wind Direction (degrees)';
        xlab='Temperature (^{o}C)';

        switch xlab
            case 'Wind Speed (m s^{-1})'
                idata_col=3;
                factor=1/100;
            case 'Temperature (^{o}C)'
                idata_col=12;
                factor=1/10; % Devide by 10 to make the decimal point from RASS data
            case 'Wind Direction (degrees)'
                idata_col=4;
                factor=1;
        end


        idat=0; %ignore this - it is a counter for the numebr of lines on the plot

        %first profile for case 2 only
        time=1; %set the time in hours required
        idat=idat+1;
        xdat(idat).x=data_structure(time+1).data(idata_col,:)*factor;
        ydat(idat).y=data_structure(time+1).data(1,:); %
        labs(idat).l='Hour 01';

        %second profile
        %           time=5; %set the time in hours required
        %          idat=idat+1;
        %         xdat(idat).x=data_structure(time+1).data(idata_col,:)*factor;
        %        ydat(idat).y=data_structure(time+1).data(1,:); %
        %       labs(idat).l='Hour 5';
        %Third profile
        %      time=10; %set the time in hours required
        %     idat=idat+1;
        %    xdat(idat).x=data_structure(time+1).data(idata_col,:)*factor;
        %   ydat(idat).y=data_structure(time+1).data(1,:); %
        %  labs(idat).l='Hour 10';

        %add more profiles here as required

        %change x-limits
        xlims=0; %1 means to restrict the x-axis to the values below, 0 means to leave alone
        xlimits=1e4*[1 3]; %set like [xmin xmax]

        zlims=0;  %1 means to restrict the x-axis to the values below, 0 means to leave alone
        zmin=0;
        zmax=10e3;

        nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
        lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane



        titlenam = ['RASS data for 13 March, 2009']; %here can change the title if you like
        %            titlenam = 'My plot';

        figname=titlenam;
        savename=figname;
        %%% end of case 2 ========================================================


    case 1  %profiles at different positions in WRF domain

        %% leave this bit for now %%

        dual=0;
        lor=-1;

        no_sort=0; %flag to stop sorting in height

        %num of markers - set to -1 for all data points to have a marker
        nmark=[-1 -1 -1 -1 -1 -1]; %if give an array then these apply for the different lines
        %    nmark=0; %set to zero for no markers

        incep=0;

        i_multi_wrf=1; %flag to say that want to plot from more than one WRF/analysis output file
        i_label_time=0; %flag to label lat/lon (=0) or not(=1)
        no_title=1; %switch off the title


        i_override_loc_lab=0;
        location_lab(1).l = ''; % Site name to be plotted in the graph as legend
        location_lab(2).l = 'B';
        location_lab(3).l = 'C';
        location_lab(4).l = 'D';
        %    location_lab(5).l = 'I';

        
%write the labels for the different WRF runs/analysis here        
        wrf_run_names(1).name = 'Run A';
        wrf_run_names(2).name = 'Run B';
        wrf_run_names(3).name = 'Run C';
        wrf_run_names(4).name = 'Run D';
        wrf_run_names(5).name = 'Run E';
        wrf_run_names(6).name = 'Run F';
        wrf_run_names(7).name = 'Run G';
        wrf_run_names(8).name = 'Run H';


        %             if i_multi_wrf==1
        %                 for ilab=1:length(rundir)
        %                     location_lab(ilab).l = rundir(ilab).dir;
        %                 end
        %             end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% start to set things after here for WRF compared with RASS  %%%%%%%
        %%%%%%%%%%%%% WRF Case 1 compared with RASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zlims=1;
        zmin=100;
        %    zmin=12000;
        zmax=3000;     %up to 3000m
        %    zmax=9000;



        %     time_array=[2 13]; %set time index here as single value for all
        %- OR, if want to plot different times for each profile
        %then make an array of the time indices, e.g. [time1 time2 time3] for each profile

        %here you you choose the location s of your profiles - can mix between lat,lon pairs or x,y pairs

        %x,y coordinates of locations to plot in km
        %    x_vals = [575 575]; %[x1 x2 x3]
        %   y_vals = [740 740]; %[y1 y2 y3]

        %x_vals=[];  %alternatively set to [] if don't want to use x,y but want to use LAT,LON instead
        %y_vals=[];  %(as will overwrite the old values if set again)

        %lat, lon coordinates of locations to plot
        %LAT_vals = [-65 -65]; %[LAT1 LAT2]
        %LON_vals = [-65 -60]; %[LON1 LON2]

        LAT_vals=[]; %set to these if are just using x,y co-ordinates
        LON_vals=[];


        %%%%%%%%%%%%%%%%%%%%% Chose which parameter to plot from WRF and RASS data %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %           var='pressure';
                 var='temperature';
        %var='equiv_potemp';
        %           var='potemp';         %potential temperature
        %    var='vapour';
        %    var='ice';
        %var='RH';
        %    var='wind speed';  %wind speed magnitude sqrt(U^2+V^2)
        %    var='wind speed component';
        %    var='cloud';
%        var='wind dir';
        %    var='Froude';
        %    var='density gradient';
        %    var='dpot/dz';
        %    var='westerly wind'


        ylab='Height (m)';

        %%%%%%%%%%%%%%% RASS plotting compared with WRF %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            iRASS_plot=0;
        iRASS_plot=1;  %turn on RASS plots


        idat=0; %ignore this - it is a counter for the numebr of lines on the plot


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   rest of the routine does the setting up and plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if iRASS_plot==1  %if you want to plot both profile (WRF and RASS) you chose 1
            switch var
                case 'wind speed'
                    idata_col=3;
                    factor=1/100;
                case 'temperature'
                    idata_col=12;
                    factor=1/10;
                case 'wind dir'
                    idata_col=4;
                    factor=1;
            end


            %%%%%%%%%% first profile for RASS only compared with WRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            %%%%%%%%%%%%%%%%%%%%%%% you need to set the time here for RASS only compared with WRF outputs %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            time=12 %set the time in hours required for first profile (RASS)
            idat=idat+1;
            xdat(idat).x=data_structure(time+1).data(idata_col,:)*factor;
            ydat(idat).y=data_structure(time+1).data(1,:); %
            labs(idat).l='OBS : 13 Feb 12:00LT'; % this is for first profile label and for the second profile see line 895

            %second profile
            %    time=12; %set the time in hours required
            %   idat=idat+1;
            %  xdat(idat).x=data_structure(time+1).data(idata_col,:)*factor;
            % ydat(idat).y=data_structure(time+1).data(1,:); %
            %labs(idat).l='Hour 12';
            %Third profile
            %           time=13; %set the time in hours required
            %          idat=idat+1;
            %         xdat(idat).x=data_structure(time+1).data(idata_col,:)*factor;
            %        ydat(idat).y=data_structure(time+1).data(1,:); %
            %       labs(idat).l='Hour 13';

            %add more profiles here as required

        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% end of RASS plotting and start of WRF plotting %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %% leave this bit for now %%

        i_previous=idat;

        dual=0;
        lor=-1;

        no_sort=0; %flag to stop sorting in height

        %num of markers - set to -1 for all data points to have a marker
        nmark=[-1 -1 -1 -1 -1 -1]; %if give an array then these apply for the different lines
        %    nmark=0; %set to zero for no markers

        incep=0;

        i_multi_wrf=1; %flag to say that want to plot from more than one WRF output file
        i_label_time=0; %flag to label lat/lon (=0) or not(=1)
        no_title=1; %switch off the title


        i_override_loc_lab=0;
        location_lab(1).l = ''; % Site name to be plotted in the graph as legend (A for site 1)
        location_lab(2).l = 'B';
        location_lab(3).l = 'C';
        location_lab(4).l = 'D';
        %    location_lab(5).l = 'I';


        if i_multi_wrf==1
            %                 for ilab=1:length(rundir)
            %                     location_lab(ilab).l = rundir(ilab).dir;
            %                 end
            n_multi_plots=length(nca);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% start to set TIMES and others after here for WRF only compared with RASS %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        time_array=[10];    % [2] means plot index= 2 (index 2 = 01 hour GMT and index=3 mean 02 hr) to find out
        %which time corresponding to which index number to put
        %it here [number] you type Times in MATLAB
        %command line, it will give you all WRF
        %times with it index (means first time say 00hr is index=1
        %and second time say 01hr is index=2 and so on upto 23hr)
        %and if you want to two times [2 13] means times 2 and 13
        % (set time index here as single value for all)
        %- OR, if want to plot different times for each profile
        %then make an array of the time indices, e.g. [time1 time2 time3] for each profile

        %here you you choose the location s of your profiles - can mix between lat,lon pairs or x,y pairs

        %x,y coordinates of locations to plot in km
        %            x_vals = [575 575]; %[x1 x2 x3]
        %            y_vals = [740 740]; %[y1 y2 y3]

        x_vals=[];  %alternatively set to [] if don't want to use x,y but want to use LAT,LON instead
        y_vals=[];  %(as will overwrite the old values if set again)

        %lat, lon coordinates of locations to plot
        LAT_vals = [23.9303]; % (this is site 1)if you want to plot two positions in WRF data [LAT1 LAT2][23.9303 23.9648]
        LON_vals = [38.3412]; % (this is site 1)if you want to plot two positions in WRF data [LON1 LON2][38.3412 38.3018]
        %LAT_vals = [23.969327]; %(this is site 3) if you want to plot two positions in WRF data [LAT1 LAT2][23.9303 23.9648]
        %LON_vals = [38.304498]; % (this is site 3)if you want to plot two positions in WRF data [LON1 LON2][38.3412 38.3018]
        %           LAT_vals=[]; %set to these if are just using x,y co-ordinates
        %           LON_vals=[];
        
        LAT_vals=[-65];
        LON_vals=[-65];


        nlat = size(lat2d.var,1);
        nlon = size(lat2d.var,2);
        dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
        dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

        i_grid = dx_grid * [0:nlon-1];
        j_grid = dy_grid * [0:nlat-1];

        LAT_vals2=[];
        LON_vals2=[];
        for iflt=1:length(x_vals)
            i_extra = findheight_nearest(i_grid,x_vals(iflt));
            j_extra = findheight_nearest(j_grid,y_vals(iflt));
            LAT_vals2(iflt) = lat2d.var(j_extra,i_extra);
            LON_vals2(iflt) = lon2d.var(j_extra,i_extra);
        end

        LAT = [LAT_vals LAT_vals2];
        LON = [LON_vals LON_vals2];


        %%%%%%% get the indices for the points (LOCATION) required but
        %%%%%%% you need to remove the earlier defined values in line
        %%%%%%% 336
        [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);   % Here you can chose another
        %location to plot by replacing LAT and LON by the values
        %such as LAT=23.9303 and LON=38.3412
        %            nca(1).nc=nc;
        is_met_em(1)=0;
        
        if isnan(ilat) | isnan(ilon)
            disp('Latitude or longitude or (x,y) are out of range!');
            break
        end






        for iloc=1:length(ilat)
            
            if length(time_array)==1
                time=time_array;
            else
                time=time_array(iloc);
            end
            
            
            for irun=1:n_multi_plots
                idat=idat+1;

                location_lab(idat).l = wrf_run_names(irun).name;

                Times=(nca(irun).nc{'Times'}(:));
                




                tstr=Times(time,:);
                iund=findstr('_',tstr);
                tstr(iund)=' ';
                
                switch file_type2{irun}
                    case 'wrfout'
                        [year,month,day,hour,mins,sec,month_text]=WRF_time_strings(Times,time);
                        is_met_em=0;
                        ydat(idat).y = WRFUserARW(nca(irun).nc,'Z',time,ilat(iloc),ilon(iloc));
                    case 'met_em'
                        [year,month,day,hour,mins,sec,month_text]=WRF_time_strings(Times',1);
                        is_met_em=1;
                        ydat(idat).y = get_height_wrf_analysis(nca(irun).nc,ilat,ilon,iloc);
                end

                

                %	ydat(idat).y = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));


                %%%%%%%  gets the data for the plots below for the different cases
                switch var
                    case 'pressure'
                        figname=['Pressure profile at ' tstr ' for ' filestr];

                        if is_met_em(1)
                            xdat(idat).x = nca(irun).nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100; %prob ignore this
                        else
                            xdat(idat).x = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));
                        end
                        % iydir = -1; %reverse the direction of the pressure axis so is right way around

                    case 'equiv_potemp'
                        figname=['Equivalent potential temperature profile at ' tstr ' for ' filestr];

                        if is_met_em(1)
                            %                        xdat(idat).x = nca(irun).nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                        else
                            potemp = nca(irun).nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                            P = nca(irun).nc{'P'}(time,:,ilat(iloc),ilon(iloc)) + nca(irun).nc{'PB'}(time,:,ilat(iloc),ilon(iloc));
                            T = potemp ./ ( (1e5./P).^0.286 );
                            qv = nca(irun).nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                            xdat(idat).x = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';

                        end

                        xlims=0;
                        xlimits=[275 290];


                    case 'potemp'
                        figname=['Potential temperature profile at ' tstr ' for ' filestr];
                        xlab='Potential temperature (K)';
                        if is_met_em(1)
                            %                        xdat(idat).x = nca(irun).nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                        else
                            xdat(idat).x = nca(irun).nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                        end

                        xlims=1;
                        xlimits=[270 290];
                        xlimits=[270 305];

                    case 'temperature'

                        figname=['Temperature profile at ' tstr ' for ' filestr];
                        xlab='Temperature (^{o}C)';

                        if is_met_em(1)
                            xdat(idat).x = nca(irun).nc{'TT'}(1,:,ilat(iloc),ilon(iloc)) - 273.15;
                            if incep==0
                                xdat(idat).x = xdat(idat).x(2:end);
                            end
                        else
                            xdat(idat).x = WRFUserARW(nca(irun).nc,'tc',time,ilat(iloc),ilon(iloc));
                            if size(xdat(idat).x,2)==1   % this script was added as recommeded by Dan to plot the temperature profile
                                xdat(idat).x=xdat(idat).x';
                            end
                            if size(ydat(idat).y,2)==1 % this script was added as recommeded by Dan to plot the temperature profile
                                ydat(idat).y=ydat(idat).y';
                            end

                            xdat(idat).x = [get_wrf_point_surface(nca(irun).nc,'T2',time,ilat(iloc),ilon(iloc))-273.15 xdat(idat).x];
                            terr_level = nca(irun).nc{'HGT'}(:,ilat(iloc),ilon(iloc));
                            terr_level=terr_level(time);
                            ydat(idat).y = [terr_level+2 ydat(idat).y]; %add air temp at 2 m
                        end

                        xlims=0;
                        xlimits=[-20 -9];

                    case 'vapour'
                        figname=['Vapour profile at ' tstr ' for ' filestr];
                        xlab='Vapour mixing ratio (g kg^{-1})';

                        if is_met_em(1)
                            rh = nca(irun).nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                            T = nca(irun).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                            P = nca(irun).nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                            qsat = satvappress(T,'goff','liq',P,1)/f;
                            xdat(idat).x = 1000 * rh/100 .* qsat;
                        else
                            xdat(idat).x = 1000*nca(irun).nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                        end

                    case 'ice'

                        xlab='Ice number concentration (L^{-1})';
                        %                    xlab='Snow number concentration (L^{-1})';
                        %                    xlab='Graupel number concentration (L^{-1})';
                        %                    xlab='Rain number concentration (L^{-1})';
                        %                    xlab='Cloud mixing ratio (g kg^{-1})';
                        %                    xlab='Ice mixing ratio (g kg^{-1})';
                        %                    xlab='Snow mixing ratio (g kg^{-1})';
                        %                    xlab='Graupel mixing ratio (g kg^{-1})';
                        %                    xlab='Total condensate mixing ratio (g kg^{-1})';
                        %                    xlab='Total number concentration (L^{-1})';
                        %                    xlab='Water supersaturation (%)';
                        %                    xlab='Ice supersaturation (%)';
                        figname=[xlab '  C ' tstr ' for ' filestr];

                        if is_met_em(1)
                            rh = nca(irun).nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                            T = nca(irun).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                            P = nca(irun).nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                            qsat = satvappress(T,'goff','liq',P,1)/f;
                            xdat(idat).x = 1000 * rh/100 .* qsat;
                        else

                            T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                            P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;
                            rho=density(P,T);

                            %numbers from WRF are in #/kg so multiply by the density to get #/m3 and then divide by 1000 to get #/L

                            switch xlab
                                case 'Ice number concentration (L^{-1})'
                                    xdat(idat).x = nca(irun).nc{'QNICE'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                case 'Snow number concentration (L^{-1})'
                                    xdat(idat).x = nca(irun).nc{'QNSNOW'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                case 'Graupel number concentration (L^{-1})'
                                    xdat(idat).x = nca(irun).nc{'QNGRAUPEL'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                case 'Rain number concentration (L^{-1})'
                                    xdat(idat).x = nca(irun).nc{'QNRAIN'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                case 'Ice mixing ratio (g kg^{-1})'
                                    xdat(idat).x = 1000*nca(irun).nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                                case 'Snow mixing ratio (g kg^{-1})'
                                    xdat(idat).x = 1000*nca(irun).nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                                case 'Graupel mixing ratio (g kg^{-1})'
                                    xdat(idat).x = 1000*nca(irun).nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                                case 'Cloud mixing ratio (g kg^{-1})'
                                    xdat(idat).x = 1000*nca(irun).nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                                case 'Total condensate mixing ratio (g kg^{-1})'
                                    xdat(idat).x = 1000*( nca(irun).nc{'QICE'}(time,:,ilat(iloc),ilon(iloc))+nca(irun).nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc))...
                                        +nca(irun).nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc))+nca(irun).nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc))...
                                        +nca(irun).nc{'QRAIN'}(time,:,ilat(iloc),ilon(iloc)) );
                                case 'Total number concentration (L^{-1})'
                                    xdat(idat).x = ( nca(irun).nc{'QNICE'}(time,:,ilat(iloc),ilon(iloc))+nca(irun).nc{'QNSNOW'}(time,:,ilat(iloc),ilon(iloc))...
                                        +nca(irun).nc{'QNGRAUPEL'}(time,:,ilat(iloc),ilon(iloc))+nca(irun).nc{'QNRAIN'}(time,:,ilat(iloc),ilon(iloc)) ).*rho/1000;
                                case 'Water supersaturation (%)'
                                    qv=nca(irun).nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                    qvs=satvappress(T,'goff','liq',P,1)/f;
                                    xdat(idat).x=100*(qv./qvs-1);
                                case 'Ice supersaturation (%)'
                                    qv=nca(irun).nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                    qvs=satvappress(T,'goff','ice',P,1)/f;
                                    xdat(idat).x=100*(qv./qvs-1);

                            end





                        end

                        xlims=0;
                        xlimits=[-20 1];

                    case 'RH'
                        figname=['RH profile at ' tstr ' for ' filestr];
                        xlab='Relative humidity (%)';

                        if is_met_em(1)
                            xdat(idat).x = nca(irun).nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                        else
                            qv = nca(irun).nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                            T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                            P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;
                            qsat = satvappress(T,'goff','liq',P,1)/f;
                            xdat(idat).x = qv./qsat *100;
                        end

                    case 'wind speed'
                        if is_met_em(1)
                            u=nca(idat-i_previous).nca(irun).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                            v=nca(idat-i_previous).nca(irun).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                        else
                            if i_multi_wrf==1
                                dire2 = [dire(iloc).dir rundir(iloc).dir];
                                cd(dire2);
                                u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                                v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                            else
                                u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                                v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                            end
                        end

                        xdat(idat).x= sqrt( u.^2 + v.^2 );

                        figname=['Wind speed profile at ' tstr ' for ' filestr];
                        if i_paper_labels==1
                            figname=['Wind speed profiles'];
                        end
                        xlab='Wind speed (m s^{-1})';

                        if incep==0 & is_met_em(1)==1
                            xdat(idat).x = xdat(idat).x(2:end);
                        end

                    case 'wind speed component'

                        xlims=0;
                        xlimits=[0 2.5];


                        if is_met_em(1)
                            u=nca(idat-i_previous).nca(irun).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                            v=nca(idat-i_previous).nca(irun).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                        else
                            if i_multi_wrf==1
                                dire2 = [dire(iloc).dir rundir(iloc).dir];
                                cd(dire2);
                                u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                                v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                            else
                                u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                                v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                            end
                        end

                        sp = sqrt( u.^2 + v.^2 );


                        %%%%%%%%%%%% get wind dir so can calculate easterly component
                        jnorth = ilat(iloc) + 10;
                        lons_north = lon2d.var(jnorth,:);
                        [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

                        %angle of the local north line relative to the grid
                        thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );

                        clear dir
                        for iuv=1:length(u)

                            theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                            if u(iuv)==0 & v(iuv)==0
                                dir(iuv) = 0;
                            elseif u(iuv)>=0 & v(iuv)>=0
                                dir(iuv) = theta2;
                            elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                                dir(iuv) = 180 + theta2;
                            elseif u(iuv)<=0 & v(iuv)<=0
                                dir(iuv) = 180 + theta2;
                            elseif u(iuv)<0 & v(iuv)>0
                                dir(iuv) = 360 + theta2; %theta2 is negative
                            end




                        end

                        dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north


                        %%%%% easterly component
                        xdat(idat).x = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)


                        figname=['Westerly wind speed profile at ' tstr ' for ' filestr];
                        if i_paper_labels==1
                            figname=['Westerly wind speed profiles'];
                        end
                        xlab='Westerly wind speed'; %dimensionless

                        if incep==0 & is_met_em(1)==1
                            xdat(idat).x = xdat(idat).x(2:end);
                        end


                    case 'wind dir'
                        xlab='Wind direction (degrees)';
                        figname=['Wind dir profile at ' tstr ' for ' filestr];

                        if is_met_em(1)
                            u=nca(irun).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                            v=nca(irun).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                        else
                            if i_multi_wrf==1
%                                dire2 = [dire(iloc).dir rundir(iloc).dir];
%                                cd(dire2);
                            end
                            u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                            v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                        end

                        jnorth = ilat(iloc) + 10;
                        lons_north = lon2d.var(jnorth,:);
                        [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

                        %angle of the local north line relative to the grid
                        thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );




                        for iuv=1:length(u)

                            theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                            if u(iuv)==0 & v(iuv)==0
                                xdat(idat).x(iuv) = 0;
                            elseif u(iuv)>=0 & v(iuv)>=0
                                xdat(idat).x(iuv) = theta2;
                            elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                                xdat(idat).x(iuv) = 180 + theta2;
                            elseif u(iuv)<=0 & v(iuv)<=0
                                xdat(idat).x(iuv) = 180 + theta2;
                            elseif u(iuv)<0 & v(iuv)>0
                                xdat(idat).x(iuv) = 360 + theta2; %theta2 is negative
                            end




                        end



                        xdat(idat).x = xdat(idat).x + 180 - thetaN*180/pi; %add 180 to make it the direction wind is coming from
                        % take away thetaN to give direction relative to north
                        i360 = find(xdat(idat).x>=360);
                        xdat(idat).x(i360) = xdat(idat).x(i360) - 360;

                        if incep==0 & is_met_em(1)==1
                            xdat(idat).x = xdat(idat).x(2:end);
                        end

                    case 'cloud'
                        cloud=nca(irun).nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                        cloud=cloud+nca(irun).nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                        cloud=cloud+nca(irun).nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                        cloud=cloud+nca(irun).nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                        xdat(idat).x=1000*cloud;


                        figname=['Cloud mixing ratio profile at ' tstr ' for ' filestr];
                        xlab='Total condensed water mixing ratio (g kg^{-1})';


                    case 'density gradient'

                        xlims=0;
                        xlimits=[0 2.5];

                        %get pressure in Pa
                        if is_met_em(1)
                            P = nca(irun).nc{'PRES'}(time,:,ilat(iloc),ilon(iloc));
                        else
                            P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc))*100;
                        end
                        %get temperature in K
                        if is_met_em(1)
                            T = nca(idat-i_previous).nca(irun).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                        else
                            T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc))+273.15;
                        end

                        %density gradient
                        xdat(idat).x = diff(density(P,T))./diff(ydat(idat).y);
                        ydat(idat).y = ydat(idat).y(2:end);

                        %density
                        %                xdat(idat).x = density(P,T);




                        figname=['Density gradient profile at ' tstr ' for ' filestr];
                        if i_paper_labels==1
                            figname=['Density gradient profiles'];
                        end
                        xlab='Density gradient (kg m^{-4})';


                        if incep==0 & is_met_em(1)==1
                            xdat(idat).x = xdat(idat).x(2:end);
                        end

                    case 'dpot/dz'

                        xlims=0;
                        xlimits=1000*[0 0.025];

                        %%%% potemp
                        if is_met_em(1)
                            %                        pot = nca(irun).nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                        else
                            pot = nca(irun).nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                        end


                        xdat(idat).x = 1000*(diff(pot)./diff(ydat(idat).y));
                        ydat(idat).y = ydat(idat).y(2:end);




                        figname=['Potential temperature gradient profile at ' tstr ' for ' filestr];
                        if i_paper_labels==1
                            figname=['Potential temperature gradient profiles'];
                        end
                        xlab='Potential temperature gradient (K km^{-1})'; %




                        if incep==0 & is_met_em(1)==1
                            xdat(idat).x = xdat(idat).x(2:end);
                        end

                    case 'westerly wind'

                        if is_met_em(1)
                            u=nca(idat-i_previous).nca(irun).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                            v=nca(idat-i_previous).nca(irun).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                        else
                            if i_multi_wrf==1
                                dire2 = [dire(iloc).dir rundir(iloc).dir];
                                cd(dire2);
                                u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                                v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                            else
                                u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                                v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));

                            end
                        end

                        sp = sqrt( u.^2 + v.^2 );


                        %%%%%%%%%%%% get wind dir so can calculate easterly component
                        jnorth = ilat(iloc) + 10;
                        lons_north = lon2d.var(jnorth,:);
                        [temp inorth] = min( abs(lons_north - lon2d.var(ilat(iloc),ilon(iloc)) ) );

                        %angle of the local north line relative to the grid
                        thetaN = atan ( (inorth - ilon(iloc)) / (jnorth - ilat(iloc)) );

                        clear dir
                        for iuv=1:length(u)

                            theta2 = 180/pi * atan ( u(iuv) ./ v(iuv) );

                            if u(iuv)==0 & v(iuv)==0
                                dir(iuv) = 0;
                            elseif u(iuv)>=0 & v(iuv)>=0
                                dir(iuv) = theta2;
                            elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                                dir(iuv) = 180 + theta2;
                            elseif u(iuv)<=0 & v(iuv)<=0
                                dir(iuv) = 180 + theta2;
                            elseif u(iuv)<0 & v(iuv)>0
                                dir(iuv) = 360 + theta2; %theta2 is negative
                            end




                        end

                        dir = dir*pi/180 - thetaN; %subtract thetaN to make it the bearing from north


                        %%%%% easterly component
                        xdat(idat).x = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)





                        figname=['Westerly wind component profile at ' tstr ' for ' filestr];
                        if i_paper_labels==1
                            figname=['Westerly wind component profiles'];
                        end
                        xlab='Westerly wind component (m s^{-1})'; %




                        if incep==0 & is_met_em(1)==1
                            xdat(idat).x = xdat(idat).x(2:end);
                        end










                end

                if no_sort==0;
                    [ydat(idat).y I]=sort(ydat(idat).y);
                    xdat(idat).x = xdat(idat).x(I);
                end

                if ydat(idat).y(1)<0 & is_met_em(1)==1 & incep==1
                    ydat(idat).y=ydat(idat).y(2:end);
                    xdat(idat).x=xdat(idat).x(2:end);
                end



                abc=['ABCDEFGHIJKLM'];

                labs(idat).l=['WRF :', location_lab(idat-i_previous).l ' ' day ' ' month_text ' ' hour ':' mins];  %simple label more suitable for papers
                % Previously this command was written ['WRF:', location_lab(idat-i_previous).l ' ' day ' ' month_text ' ' hour ':' mins]

            end %for irun=1....

        end
        %        labs(idat).l=abc(idat);


end






savename=figname;
titlenam = figname;


if no_title==1
    titlenam='';
end





if zlims==0
    zmin='';
    zmax='';
end

if noplot==0
    if subplotting==0
        hf=figure('name',figname,'Position',posit);
        fsize=18;
        fsize=14;
        %        fsize=26;
        %        fsize=30;
        ixlab=1;
    else
        if nsub==1
            posit(4)=length(idirs)*posit(4)/2.3;
            hf=figure('name',figname,'Position',posit);
        end
        subplot(xsub,ysub,nsub);
        fsize=12;
        if nsub==length(idirs)
            ixlab=1;
        else
            ixlab=0;
        end
    end

    %now call the plotting function
    [h,ax,ax2]=plotXY6(xdat,ydat,labs,nmark,lwidth,lor,logflag,xlab,ylab,[zmin zmax],1,dual,secyA,secyB,lab2,fsize,ixlab,ixdir,iydir,xloc);
    %function [H1,ax2]=plotXY3(xdat,ydat,labs,nmark,lwidth,leglor,logflag,xlab,ylab,ylims,zline) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers
    %put nmark as -1 for markers for all points
    % LEGEND(...,Pos) places the legend in the specified
    %     location:
    %         0 = Automatic "best" placement (least conflict with data)
    %         1 = Upper right-hand corner (default)
    %         2 = Upper left-hand corner
    %         3 = Lower left-hand corner
    %         4 = Lower right-hand corner
    %        -1 = To the right of the plot
    if gridon==1
        grid on;
    end

    if xlims==1
        set(gca,'xlim',xlimits);
    end

    if (ixtime==1)
        xx=get(gca,'xticklabels');
        xx=str2num(xx);
        xx=num2str(mod(xx,24));
        set(gca,'xticklabels',xx);
    end

    if (ititle==1)
        %    title(titlenam,'fontsize',18,'verticalalignment','middle');

        htit=title(titlenam,'fontsize',fsize,'verticalalignment','baseline');

        if dual==2
            ptit=get(htit,'position');
            set(htit,'position',[ptit(1) ptit(2)*1.015]);
        end
        %        title(titlenam,'fontsize',18,'verticalalignment','top');

        %    title(titlenam,'fontsize',18);
    end

    if idirlabel==1
        xlims=get(gca,'xlim');
        ylims=get(gca,'ylim');
        text(xlims(1),ylims(1)-(ylims(2)-ylims(1))/13,direcDan(idir).dir);
    end

    if add_points==1
        add_points_to_plot(xpos,ypos,point_labs,8,11)
    end

    if iaxis_square==1;
        axis square
    end

end


