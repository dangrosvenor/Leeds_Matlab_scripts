WRF_plot_simple_set_flags  %sets up flags and clears data

%select case to run
graph=1; %profiles at different locations
%graph=2; %our example x,y plot 

    switch graph
        case 9999999
            % template to copy for new graph
            
            time=3; %time index                                    

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
            zmin=1500;
            zmax=3000;
            
            nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

            
            tstr=Times(time,:);
            iund=findstr('_',tstr);
            tstr(iund)=' ';
            titlenam = ['XXX for ' tstr];

            figname=titlenam;
            savename=figname;

    case 2 %general plot of x against y
            
            time=3; %time index                                    

            ylab='Height (m)';
            xlab= 'Random plot data (m)';

            
            idat=0; %ignore this - it is a counter for the numebr of lines on the plot

            HGT=WRFUserARW(nca(1).nc,'Z',time,ilat(iloc),ilon(iloc)); %gets height data from WRF output
%first profile
            idat=idat+1;
            xdat(idat).x=2*HGT; %
            ydat(idat).y=HGT; %
            labs(idat).l='Jaafar';

%second profile            
            idat=idat+1;
            xdat(idat).x=0.5*HGT; %x data
            ydat(idat).y=HGT; %y data
            labs(idat).l='Rudra';

%add more profiles here as required

%change x-limits
            xlims=1; %1 means to restrict the x-axis to the values below, 0 means to leave alone
            xlimits=1e4*[1 3]; %set like [xmin xmax]

            zlims=1;  %1 means to restrict the x-axis to the values below, 0 means to leave alone
            zmin=0;
            zmax=10e3;
            
            nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.
            lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

            
            tstr=Times(time,:);
            iund=findstr('_',tstr);
            tstr(iund)=' ';
            
            
            titlenam = ['XXX data for ' tstr]; %here can change the title if you like
%            titlenam = 'My plot';

            figname=titlenam;
            savename=figname;
            

        case 1  %profiles at different positions in WRF domain
            
            %% leave this bit for now %%

            dual=0;
            lor=-1;

            no_sort=0; %flag to stop sorting in height

            %num of markers - set to -1 for all data points to have a marker
            nmark=[-1 -1 -1 -1 -1 -1]; %if give an array then these apply for the different lines
            %    nmark=0; %set to zero for no markers

            incep=0;

            i_multi_wrf=0; %flag to say that want to plot from more than one WRF output file
            i_label_time=0; %flag to label lat/lon (=0) or not(=1)
            no_title=1; %switch off the title


            i_override_loc_lab=0;
            location_lab(1).l = 'A';
            location_lab(2).l = 'B';
            location_lab(3).l = 'C';
            location_lab(4).l = 'D';
            %    location_lab(5).l = 'I';


            if i_multi_wrf==1
                for ilab=1:length(rundir)
                    location_lab(ilab).l = rundir(ilab).dir;
                end
            end

%%%%%%%% start to set things after here  %%%%%%%

            zlims=1;
            zmin=0;
            %    zmin=12000;
            zmax=5000;
            %    zmax=9000;



            time_array=[3 3 3]; %set time index here as single value for all
                                %- OR, if want to plot different times for each profile
                                %then make an array of the time indices, e.g. [time1 time2 time3] for each profile
            
%here you you choose the location s of your profiles - can mix between lat,lon pairs or x,y pairs

%x,y coordinates of locations to plot in km
            x_vals = [500 500 500]; %[x1 x2]
            y_vals = [200 300 400]; %[y1 y2]  
            
            %x_vals=[];  %alternatively set to [] if don't want to use x,y but want to use LAT,LON instead
            %y_vals=[];  %(as will overwrite the old values if set again)
            
            %lat, lon coordinates of locations to plot
            %LAT_vals = [-65 -65]; %[LAT1 LAT2]
            %LON_vals = [-65 -60]; %[LON1 LON2]
                        
           LAT_vals=[]; %set to these if are just using x,y co-ordinates
           LON_vals=[];
      

%            var='pressure';
 %           var='temperature';
            %var='equiv_potemp';
            var='potemp';         %potential temperature
            %    var='vapour';
            %    var='ice';
            %    var='RH';
  %              var='wind speed';  %wind speed magnitude sqrt(U^2+V^2)
            %    var='wind speed component';
            %    var='cloud';
            %    var='wind dir';
            %    var='Froude';
            %    var='density gradient';
            %    var='dpot/dz';
            %    var='westerly wind'


            ylab='Height (m)';
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   rest of the routine does the setting up and plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







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


            %%%%%%% get the indices for the points required
            [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);

            nca(1).nc=nc;
            is_met_em(1)=0;

            i=0;
            i_previous=i;




                for iloc=1:length(ilat)

                    if length(time_array)==1
                        time=time_array;
                    else
                        time=time_array(iloc);
                    end



                    i=i+1;

                    Times=(nca(1).nc{'Times'}(:));



                    tstr=Times(time,:);
                    iund=findstr('_',tstr);
                    tstr(iund)=' ';

                    [year,month,day,hour,mins,sec,month_text]=WRF_time_strings(Times,time);

                    ydat(i).y = WRFUserARW(nca(1).nc,'Z',time,ilat(iloc),ilon(iloc));

                    %	ydat(i).y = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));
                    
                    
%%%%%%%  gets the data for the plots below for the different cases
                    switch var
                        case 'pressure'
                            figname=['Pressure profile at ' tstr ' for ' filestr];

                            if is_met_em(1)
                                xdat(i).x = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100; %prob ignore this
                            else
                                xdat(i).x = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));
                            end
                            % iydir = -1; %reverse the direction of the pressure axis so is right way around

                        case 'equiv_potemp'
                            figname=['Equivalent potential temperature profile at ' tstr ' for ' filestr];

                            if is_met_em(1)
                                %                        xdat(i).x = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                            else
                                potemp = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                                P = nc{'P'}(time,:,ilat(iloc),ilon(iloc)) + nc{'PB'}(time,:,ilat(iloc),ilon(iloc));
                                T = potemp ./ ( (1e5./P).^0.286 );
                                qv = nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                xdat(i).x = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';

                            end

                            xlims=0;
                            xlimits=[275 290];


                        case 'potemp'
                            figname=['Potential temperature profile at ' tstr ' for ' filestr];
                            xlab='Potential temperature (K)';
                            if is_met_em(1)
                                %                        xdat(i).x = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                            else
                                xdat(i).x = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                            end

                            xlims=1;
                            xlimits=[270 290];
                            xlimits=[270 305];

                        case 'temperature'

                            figname=['Temperature profile at ' tstr ' for ' filestr];
                            xlab='Temperature (^{o}C)';

                            if is_met_em(1)
                                xdat(i).x = nca(i-i_previous).nc{'TT'}(1,:,ilat(iloc),ilon(iloc)) - 273.15;
                                if incep==0
                                    xdat(i).x = xdat(i).x(2:end);
                                end
                            else
                                xdat(i).x = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc));
                                xdat(i).x = [get_wrf_point_surface(nc,'T2',time,ilat(iloc),ilon(iloc))-273.15 xdat(i).x];
                                terr_level = nc{'HGT'}(:,ilat(iloc),ilon(iloc));
                                terr_level=terr_level(time);
                                ydat(i).y = [terr_level+2 ydat(i).y]; %add air temp at 2 m
                            end

                            xlims=0;
                            xlimits=[-20 -9];

                        case 'vapour'
                            figname=['Vapour profile at ' tstr ' for ' filestr];
                            xlab='Vapour mixing ratio (g kg^{-1})';

                            if is_met_em(1)
                                rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                                T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                                P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                                qsat = satvappress(T,'goff','liq',P,1)/f;
                                xdat(i).x = 1000 * rh/100 .* qsat;
                            else
                                xdat(i).x = 1000*nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
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
                                rh = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                                T = nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                                P = nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                                qsat = satvappress(T,'goff','liq',P,1)/f;
                                xdat(i).x = 1000 * rh/100 .* qsat;
                            else

                                T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                                P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;
                                rho=density(P,T);

                                %numbers from WRF are in #/kg so multiply by the density to get #/m3 and then divide by 1000 to get #/L

                                switch xlab
                                    case 'Ice number concentration (L^{-1})'
                                        xdat(i).x = nc{'QNICE'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                    case 'Snow number concentration (L^{-1})'
                                        xdat(i).x = nc{'QNSNOW'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                    case 'Graupel number concentration (L^{-1})'
                                        xdat(i).x = nc{'QNGRAUPEL'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                    case 'Rain number concentration (L^{-1})'
                                        xdat(i).x = nc{'QNRAIN'}(time,:,ilat(iloc),ilon(iloc)).*rho/1000;
                                    case 'Ice mixing ratio (g kg^{-1})'
                                        xdat(i).x = 1000*nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                                    case 'Snow mixing ratio (g kg^{-1})'
                                        xdat(i).x = 1000*nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                                    case 'Graupel mixing ratio (g kg^{-1})'
                                        xdat(i).x = 1000*nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                                    case 'Cloud mixing ratio (g kg^{-1})'
                                        xdat(i).x = 1000*nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                                    case 'Total condensate mixing ratio (g kg^{-1})'
                                        xdat(i).x = 1000*( nc{'QICE'}(time,:,ilat(iloc),ilon(iloc))+nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc))...
                                            +nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc))+nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc))...
                                            +nc{'QRAIN'}(time,:,ilat(iloc),ilon(iloc)) );
                                    case 'Total number concentration (L^{-1})'
                                        xdat(i).x = ( nc{'QNICE'}(time,:,ilat(iloc),ilon(iloc))+nc{'QNSNOW'}(time,:,ilat(iloc),ilon(iloc))...
                                            +nc{'QNGRAUPEL'}(time,:,ilat(iloc),ilon(iloc))+nc{'QNRAIN'}(time,:,ilat(iloc),ilon(iloc)) ).*rho/1000;
                                    case 'Water supersaturation (%)'
                                        qv=nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                        qvs=satvappress(T,'goff','liq',P,1)/f;
                                        xdat(i).x=100*(qv./qvs-1);
                                    case 'Ice supersaturation (%)'
                                        qv=nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                        qvs=satvappress(T,'goff','ice',P,1)/f;
                                        xdat(i).x=100*(qv./qvs-1);

                                end





                            end

                            xlims=0;
                            xlimits=[-20 1];

                        case 'RH'
                            figname=['RH profile at ' tstr ' for ' filestr];
                            xlab='Relative humidity (%)';

                            if is_met_em(1)
                                xdat(i).x = nc{'RH'}(1,:,ilat(iloc),ilon(iloc));
                            else
                                qv = nc{'QVAPOR'}(time,:,ilat(iloc),ilon(iloc));
                                T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc)) + 273.15;
                                P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc)) *100;
                                qsat = satvappress(T,'goff','liq',P,1)/f;
                                xdat(i).x = qv./qsat *100;
                            end

                        case 'wind speed'
                            if is_met_em(1)
                                u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                                v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
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

                            xdat(i).x= sqrt( u.^2 + v.^2 );

                            figname=['Wind speed profile at ' tstr ' for ' filestr];
                            if i_paper_labels==1
                                figname=['Wind speed profiles'];
                            end
                            xlab='Wind speed (m s^{-1})';

                            if incep==0 & is_met_em(1)==1
                                xdat(i).x = xdat(i).x(2:end);
                            end

                        case 'wind speed component'

                            xlims=0;
                            xlimits=[0 2.5];


                            if is_met_em(1)
                                u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                                v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
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
                            xdat(i).x = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)


                            figname=['Westerly wind speed profile at ' tstr ' for ' filestr];
                            if i_paper_labels==1
                                figname=['Westerly wind speed profiles'];
                            end
                            xlab='Westerly wind speed'; %dimensionless

                            if incep==0 & is_met_em(1)==1
                                xdat(i).x = xdat(i).x(2:end);
                            end


                        case 'wind dir'
                            xlab='Wind direction (degrees)';
                            figname=['Wind dir profile at ' tstr ' for ' filestr];

                            if is_met_em(1)
                                u=nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                                v=nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                            else
                                if i_multi_wrf==1
                                    dire2 = [dire(iloc).dir rundir(iloc).dir];
                                    cd(dire2);
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
                                    xdat(i).x(iuv) = 0;
                                elseif u(iuv)>=0 & v(iuv)>=0
                                    xdat(i).x(iuv) = theta2;
                                elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
                                    xdat(i).x(iuv) = 180 + theta2;
                                elseif u(iuv)<=0 & v(iuv)<=0
                                    xdat(i).x(iuv) = 180 + theta2;
                                elseif u(iuv)<0 & v(iuv)>0
                                    xdat(i).x(iuv) = 360 + theta2; %theta2 is negative
                                end




                            end



                            xdat(i).x = xdat(i).x + 180 - thetaN*180/pi; %add 180 to make it the direction wind is coming from
                            % take away thetaN to give direction relative to north
                            i360 = find(xdat(i).x>=360);
                            xdat(i).x(i360) = xdat(i).x(i360) - 360;

                            if incep==0 & is_met_em(1)==1
                                xdat(i).x = xdat(i).x(2:end);
                            end

                        case 'cloud'
                            cloud=nc{'QCLOUD'}(time,:,ilat(iloc),ilon(iloc));
                            cloud=cloud+nc{'QICE'}(time,:,ilat(iloc),ilon(iloc));
                            cloud=cloud+nc{'QSNOW'}(time,:,ilat(iloc),ilon(iloc));
                            cloud=cloud+nc{'QGRAUP'}(time,:,ilat(iloc),ilon(iloc));
                            xdat(i).x=1000*cloud;


                            figname=['Cloud mixing ratio profile at ' tstr ' for ' filestr];
                            xlab='Total condensed water mixing ratio (g kg^{-1})';


                        case 'density gradient'

                            xlims=0;
                            xlimits=[0 2.5];

                            %get pressure in Pa
                            if is_met_em(1)
                                P = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc));
                            else
                                P = WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc))*100;
                            end
                            %get temperature in K
                            if is_met_em(1)
                                T = nca(i-i_previous).nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                            else
                                T = WRFUserARW(nc,'tc',time,ilat(iloc),ilon(iloc))+273.15;
                            end

                            %density gradient
                            xdat(i).x = diff(density(P,T))./diff(ydat(i).y);
                            ydat(i).y = ydat(i).y(2:end);

                            %density
                            %                xdat(i).x = density(P,T);




                            figname=['Density gradient profile at ' tstr ' for ' filestr];
                            if i_paper_labels==1
                                figname=['Density gradient profiles'];
                            end
                            xlab='Density gradient (kg m^{-4})';


                            if incep==0 & is_met_em(1)==1
                                xdat(i).x = xdat(i).x(2:end);
                            end

                        case 'dpot/dz'

                            xlims=0;
                            xlimits=1000*[0 0.025];

                            %%%% potemp
                            if is_met_em(1)
                                %                        pot = nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                            else
                                pot = nc{'T'}(time,:,ilat(iloc),ilon(iloc)) + 300;
                            end


                            xdat(i).x = 1000*(diff(pot)./diff(ydat(i).y));
                            ydat(i).y = ydat(i).y(2:end);




                            figname=['Potential temperature gradient profile at ' tstr ' for ' filestr];
                            if i_paper_labels==1
                                figname=['Potential temperature gradient profiles'];
                            end
                            xlab='Potential temperature gradient (K km^{-1})'; %




                            if incep==0 & is_met_em(1)==1
                                xdat(i).x = xdat(i).x(2:end);
                            end

                        case 'westerly wind'

                            if is_met_em(1)
                                u=nca(i-i_previous).nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                                v=nca(i-i_previous).nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
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
                            xdat(i).x = sp.*sin(dir); %get the wind direcion in the east component (as is approx perpendicular to peninsula)





                            figname=['Westerly wind component profile at ' tstr ' for ' filestr];
                            if i_paper_labels==1
                                figname=['Westerly wind component profiles'];
                            end
                            xlab='Westerly wind component (m s^{-1})'; %




                            if incep==0 & is_met_em(1)==1
                                xdat(i).x = xdat(i).x(2:end);
                            end

                       





                    end

                    if no_sort==0;
                        [ydat(i).y I]=sort(ydat(i).y);
                        xdat(i).x = xdat(i).x(I);
                    end

                    if ydat(i).y(1)<0 & is_met_em(1)==1 & incep==1
                        ydat(i).y=ydat(i).y(2:end);
                        xdat(i).x=xdat(i).x(2:end);
                    end



                    abc=['ABCDEFGHIJKLM'];

                    labs(i).l=[location_lab(i-i_previous).l ' ' day ' ' month_text ' ' hour ':' mins];  %simple label more suitable for papers




                end
                %        labs(i).l=abc(i);


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

