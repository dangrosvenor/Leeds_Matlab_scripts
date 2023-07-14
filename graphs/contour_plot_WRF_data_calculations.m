


                    %********************************************************************
%%%                case 'wrf_plot'


                    xinds = 1:size(lat2d(1).var,2); %lat2d is arranged in LAT,LON order so x is second index (LON)
                    yinds = 1:size(lat2d(1).var,1);


                    zz(1).z = yinds;
                    timesTH(1).t = xinds;

                    dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
                    dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

                    zz(1).z = (zz(1).z - 1)*dy_grid;
                    timesTH(1).t = (timesTH(1).t - 1)*dx_grid;





                    %             savemem=0;
                    %             switch savemem
                    %             case 0
                    %              pdat(1).p = squeeze(WRFUserARW(nc,'Z',time,ih_wrf));
                    %		pdat(1).p = WRFUserARW(nc,'p',time,ih_wrf);
                    %		pdat(1).p = WRFUserARW(nc,'tc',time,ih_wrf);   %temperature
                    %        pdat(1).p = nc{'TT'}(time,ih_wrf,:,:)-273.15;   %temperature
                    %        pdat(1).p = nc{'RH'}(time,ih_wrf,:,:);   %RH
                    %        pdat(1).p = nc{'T'}(time,ih_wrf,:) + 300;   %potential temperature
                    %        pdat(1).p = ( nc{'P'}(time,ih_wrf,:) + nc{'PB'}(time,ih_wrf,:) ) / 100;
                    %         pdat(1).p = ( nc{'PRES'}(time,ih_wrf,:) );
                    %        pdat(1).p = squeeze(Zlev) - terrain; %height above terrain

                    %    pdat(1).p = sp_nudging_z8 - sp_no_nudging_z8;

                    %        pdat(1).p=squeeze(Zlev);
                    %		pdat(1).p = nc{'W'}(time,ih_wrf,:,:);
                    %                pdat(1).p = -nc{'LH'}(time,:,:); %latent heat (W/m2)  -minus so that positive gives energy to surface
                    %                 pdat(1).p = -nc{'HFX'}(time,:,:); %sensible heat (W/m2)

                    %                pdat(1).p = -nc{'LH'}(time,:,:) - nc{'HFX'}(time,:,:); %latent heat (W/m2)  minus so that positive gives energy to surface
                    %                pdat(1).p = nc{'TSLB'}(time,1,:,:)-273.15; %soil temperature (degC) - level one is closest to the surface, 4 is deep below

                    %                pdat(1).p = nc{'Q2'}(time,:,:)*1000; %vapour MR air temperature (degC)


                    %%%%%%%%%%%%%%%%%%%%%%%%%% difference plots %%%%%%%%%%%%%%%%%%
                    %                 pdat(1).p=temp_ecmwf-temp_ncep
                    %                 pdat(1).p=(p270_ecmwf_d02-p270_ncep_d02);
                    %                 pdat(1).p=(psfc_ecmwf_d01-psfc_ncep_d01)/100;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %                pdat(1).p = 0.01* nc{'PSFC'}(time,:,:); %surface pressure




                    %                pdat(1).p = nc{'GRDFLX'}(time,:,:);

                    %                pdat(1).p = nc{'QFX'}(time,:,:);

                    %                pdat(1).p = nc{'SWDOWN'}(time,:,:);
                    %                pdat(1).p = nc{'ALBEDO'}(time,:,:);

                    %                 pdat(1).p = nc{'SEAICE'}(time,:,:);
                    %                  pdat(1).p = nc{'LU_INDEX'}(time,:,:);
                    %                 pdat(1).p = nc{'SNOWC'}(time,:,:);

                    %                 pdat(1).p = nc{'LANDMASK'}(time,:,:); %includes the ice shelves - same in ncep and ecmwf
                    %                 pdat(1).p = nc{'LANDSEA'}(time,:,:); %large diffs - ncep looks more blocky and doesn't cover the ice shelf like ecmwf does
                    %ecmwf looks better - but perhaps LANDMASK takes precedent - this looks good and the same in both
                    %                 pdat(1).p = nc{'SEAICE'}(time,:,:);   %bit different but perhaps not enough to explain diffs - same over the ice shelf
                    %                 pdat(1).p = nc{'SLOPECAT'}(time,:,:);   %large difference especially over the ice shelf
                    %                 pdat(1).p = nc{'SNOALB'}(time,:,:);     %large changes similar to above
                    %                 pdat(1).p = nc{'GREENFRAC'}(time,1,:,:); %large changes similar to above

                    %                 pdat(1).p = nc{'ALBEDO12M'}(time,12,:,:); %%large changes similar to above
                    %                 pdat(1).p = nc{'SOILCBOT'}(time,16,:,:); % %same
                    %                 pdat(1).p = nc{'SOILCTOP'}(time,14,:,:); % %same
                    %                  pdat(1).p = nc{'SPECHUMD'}(time,:,:); % %not present for ECMWF
                    %                  pdat(1).p = nc{'SLPY'}(time,:,:); %SLPX and SLPY very similar for both ncep and ecmwf
                    %                  pdat(1).p = nc{'HGT_M'}(time,:,:); %same
                    %                  pdat(1).p = nc{'SNOW'}(time,:,:); %very different - units wrong for ecmwf - 10 kg/m2 instead of 10m!
                    %                  pdat(1).p = nc{'SM'}(time,1,:,:); %think is soil moisture - only two levels for ncep, 4 for ecmwf
                    %                  pdat(1).p = nc{'ST'}(time,1,:,:)-273.15; %soil temps - since have 4 depth layers in ecmwf, 2 in ncep

                    %                pdat(1).p = nc{'ST010200'}(time,:,:);
                    %                pdat(1).p = nc{'ST000010'}(time,:,:);
                    %                pdat(1).p = nc{'SM000010'}(time,:,:);
                    %                pdat(1).p = nc{'SM010200'}(time,:,:);

                    %                pdat(1).p = nc{'ST010200'}(time,:,:);
                    %                pdat(1).p = nc{'ST000010'}(time,:,:);
                    %                pdat(1).p = nc{'SM000010'}(time,:,:);
                    %                pdat(1).p = nc{'ST100255'}(time,:,:);
                    %                 pdat(1).p = nc{'SOILTEMP'}(time,:,:)-273.15;
                    %                 pdat(1).p = nc{'LANDUSEF'}(time,16,:,:);
                    %                 pdat(1).p = nc{'SOILHGT'}(time,:,:); %large diffs between ncep and ecmwf


                    %                pdat(1).p = nc{'SKINTEMP'}(time,:,:) - 273.15; %skin temperature for met_em files
                    %                pdat(1).p = nc{'SST'}(time,:,:) - 273.15; %skin temperature for met_em files - only for ecmwf files - values are very strange
                    %range from 1.6e29 to -1.2e30



                    %                pdat(1).p = 0.01* nc{'PMSL'}(time,:,:); %pressure

                    %                pdat(1).p = nc{'TSLB'}(time,1,:,:) - 273.15; %skin temperature for wrfout files
                    %                 pdat(1).p = nc{'XLAND'}(time,:,:);

                    %                pdat(1).p = nc{'ALBEDO'}(time,:,:); %
                    %                pdat(1).p = nc{'SOILTB'}(time,:,:); %
                    %                pdat(1).p = nc{'VEGFRA'}(time,:,:);
                    %                pdat(1).p = nc{'ISLTYP'}(time,:,:);
                    %                pdat(1).p = f* (nc{'QICE'}(itime,ih_wrf,:,:)+nc{'QSNOW'}(itime,ih_wrf,:,:)+nc{'QGRAUP'}(itime,ih_wrf,:,:)+nc{'QVAPOR'}(itime,ih_wrf,:,:) );



                    P = 100*WRFUserARW(nc,'p',time,ih_wrf);    %Pa
                    T = WRFUserARW(nc,'tc',time,ih_wrf)+273.15;   %temperature K
                    rho=density(P,T);
                    %        pdat(1).p = (nc{'QNRAIN'}(time,ih_wrf,:,:)+nc{'QNICE'}(time,ih_wrf,:,:)+nc{'QNSNOW'}(time,ih_wrf,:,:)+nc{'QNGRAUPEL'}(time,ih_wrf,:,:)).*rho/1000; %total number
                    %        pdat(1).p = (nc{'QNICE'}(time,ih_wrf,:,:)).*rho/1000; %total number
                    %        pdat(1).p = (nc{'QNSNOW'}(time,ih_wrf,:,:)).*rho/1000; %total number
                    %        pdat(1).p = (nc{'QNRAIN'}(time,ih_wrf,:,:)).*rho/1000; %total number
                    %        pdat(1).p = (nc{'QNGRAUPEL'}(time,ih_wrf,:,:)).*rho/1000; %total number

                    %        pdat(1).p = 1000*(nc{'QCLOUD'}(time,ih_wrf,:,:));
                    %        pdat(1).p = 1000*(nc{'QICE'}(time,ih_wrf,:,:));
                    %        pdat(1).p = 1000*(nc{'QSNOW'}(time,ih_wrf,:,:));
                    %        pdat(1).p = 1000*(nc{'QVAPOR'}(time,ih_wrf,:,:));
                    %        pdat(1).p = 1000*(nc{'QGRAUP'}(time,ih_wrf,:,:));
                    %        pdat(1).p = 1000*(nc{'QRAIN'}(time,ih_wrf,:,:));

                    %        pdat(1).p = (nc{'SNOWNC'}(time,:,:));
                    %        pdat(1).p = (nc{'RAINNC'}(time,:,:));

                    %        pdat(1).p = (nc{'RAINNC'}(time,:,:)) - (nc{'RAINNC'}(time-1,:,:));   %precip tendency (since the last output time)

                    %        pdat(1).p = 1000*(nc{'QCLOUD'}(time,ih_wrf,:,:)+nc{'QICE'}(time,ih_wrf,:,:)+nc{'QSNOW'}(time,ih_wrf,:,:)+nc{'QGRAUP'}(time,ih_wrf,:,:)); %total condensate
                    %        pdat(1).p = 1000*(nc{'QVAPOR'}(time,ih_wrf,:,:)+nc{'QCLOUD'}(time,ih_wrf,:,:)+nc{'QICE'}(time,ih_wrf,:,:)+nc{'QSNOW'}(time,ih_wrf,:,:)+nc{'QGRAUP'}(time,ih_wrf,:,:)); %total condensate


                    %N.B. when setting short_plot_name earlier, take off the file_str at the end but keep on here (below)
                    switch short_plot_name
                        case  ['Water supersaturation (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
                            qv=nc{'QVAPOR'}(time,ih_wrf,:,:);
                            qvs=satvappress(T,'goff','liq',P,1)/f;
                            pdat(1).p =100*(qv./qvs-1);
                        case ['Ice supersaturation (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
                            qv=nc{'QVAPOR'}(time,ih_wrf,:,:);
                            qvs=satvappress(T,'goff','ice',P,1)/f;
                            pdat(1).p=100*(qv./qvs-1);
                        case ['Temperature (^{o}C) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
                            pdat(1).p=WRFUserARW(nc,'tc',time,ih_wrf);
                        case ['Hallet Mossop flag at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr]
                            qv=nc{'QVAPOR'}(time,ih_wrf,:,:);
                            qvs=satvappress(T,'goff','ice',P,1)/f;
                            si=100*(qv./qvs-1);

                            qvsw=satvappress(T,'goff','liq',P,1)/f;
                            sw =100*(qv./qvsw-1);

                            ql=1000*nc{'QCLOUD'}(time,ih_wrf,:,:); %g/kg
                            qi=1000*nc{'QICE'}(time,ih_wrf,:,:);

                            tc=WRFUserARW(nc,'tc',time,ih_wrf);
                            ihm=find(tc<=-3 & tc>=-8 & ql>0.01 & qi>2e-4);
                            pdat(1).p=zeros(size(tc));
                            pdat(1).p(ihm)=1;

                        case ['Sea Ice flag for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                            pdat(1).p=nc{'SEAICE'}(time,:,:);
                        case ['Skin temperature for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (^{o}C) for ' filestr];
                            pdat(1).p = nc{'TSK'}(time,:,:) - 273.15; %skin temperature for wrfout files
                        case ['Wind speed (m s^{-1}) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                            u=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
                            v=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
                            pdat(1).p = sqrt(u.^2+v.^2);
                        case ['Surface pressure (hPa) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                            pdat(1).p = 0.01* nc{'PSFC'}(time,:,:); %surface pressure
                        case ['Pressure at level ' num2str(ih_wrf) ' for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC (hPa) for ' filestr]
                            pdat(1).p = P/100;
                        case ['Relative Humidity (%) at level ' num2str(ih_wrf) ' (~' medZ 'm above terrain) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                            qv = nc{'QVAPOR'}(time,ih_wrf,:,:); %get vapour
                            qsat = satvappress(T,'goff','liq',P,1)/f;
                            pdat(1).p = qv./qsat *100;
                        case ['2m air temperature (^{o}C) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                            pdat(1).p = nc{'T2'}(time,:,:)-273.15; %2m air temperature (degC)
                        case ['10m wind speed (m s^{-1}) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                             u10 = nc{'U10'}(time,:,:); %10 m winds
                             v10 = nc{'U10'}(time,:,:); %
                             pdat(1).p = sqrt(u10.^2+v10.^2);
                         case ['2m relative humidity (%) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr];
                            qv = nc{'Q2'}(time,:,:); %get vapour
                            T2=nc{'T2'}(time,:,:); %2m air temperature (K)
                            TH2=nc{'TH2'}(time,:,:); %potential temperature
                            P2= (T2./TH2 * 1000e2^0.286 ).^(1/0.286);
                            qsat = satvappress(T2,'goff','liq',P2,1)/f;
                            pdat(1).p = qv./qsat *100;   
                        case ['Terrain height ' Times(time,9:10) ' Jan ' Times(time,12:16) ' (m) for ' filestr];
                            pdat(1).p = nc{'HGT'}(time,:,:); %terrain height
                    end

'here';
                    %pdat(1).p = nc{'GLW'}(time,:,:); %Downward LW flux
                    %pdat(1).p = lwdown_tot; %Net LW flux - mean over whole 3 days
                    %pdat(1).p = temp_tot; %Net LW flux
                    %pdat(1).p = nc{'SWDOWN'}(time,:,:); %Downward SW flux
                    %pdat(1).p = nc{'OLR'}(time,:,:); %TOA outgoing LW

                    %pdat(1).p = melt_tot_total; %total melt in mm
                    %pdat(1).p = sw_tot; %
                    %pdat(1).p = lw_tot; %
                    %pdat(1).p = 1000*cond_tot; %average total condensate (column total g/m2)
                    % pdat(1).p =  temp_tot-273.15; %
                    %pdat(1).p = sh_tot; %average melt in mm for 12,15,18 and 21 UTC
                    %pdat(1).p = rh_iz_tot; %average relative humidity


                    % iadd_overlay=0;
                    % i2_save2=i2_save(:);
                    % i2_save2(i2_save2==0)=[];
                    % i2_save2=unique(i2_save2);
                    % [ix_overlay, iy_overlay] = ind2sub(size(melt_tot),i2_save2);
                    % x_overlay=timesTH(1).t(iy_overlay);
                    % y_overlay=zz(1).z(ix_overlay);

                    % [ix_overlay, iy_overlay] = find(n_tot2_restrict>=1);
                    % x_overlay=timesTH(1).t(ix_overlay);
                    % y_overlay=zz(1).z(iy_overlay);

                    iwind_field=0;
                    if iwind_field==1
                        u=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
                        v=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
                        pdat(1).p = sqrt(u.^2+v.^2);
                    end

                    i_cloud_layer=0;
                    if i_cloud_layer==1;
                        pdat(1).p = 1000*nc{'QCLOUD'}(time,ih_wrf,:,:) + 1000*nc{'QRAIN'}(time,ih_wrf,:,:)...
                            + 1000*nc{'QICE'}(time,ih_wrf,:,:) + 1000*nc{'QSNOW'}(time,ih_wrf,:,:)...
                            + 1000*nc{'QGRAUP'}(time,ih_wrf,:,:);
                    end

                    i_condensate=0;
                    recalc=1;
                    if i_condensate==1
                        if recalc==1
                            ih_inds=1:10;

                            M = 28.97*1.67E-27;
                            k = 1.38E-23;
                            P = WRFUserARW(nc,'p',time); P.var=P.var*100;
                            T = WRFUserARW(nc,'tc',time); T.var=T.var+273.15;
                            rho=P.var(ih_inds,:,:).*M./k./T.var(ih_inds,:,:);
                            Z = squeeze(WRFUserARW(nc,'Z',time));
                            dz_grid = Z.var(ih_inds+1,:,:)-Z.var(ih_inds,:,:);



                            cond_dat = squeeze (sum( (nc{'QCLOUD'}(time,ih_inds,:,:) + nc{'QRAIN'}(time,ih_inds,:,:)...
                                + nc{'QICE'}(time,ih_inds,:,:) + nc{'QSNOW'}(time,ih_inds,:,:)...
                                + nc{'QGRAUP'}(time,ih_inds,:,:)) .* rho*1000*dx_grid*1000*dy_grid.*dz_grid ,1) );
                        else
                            fprintf(1,'\n****** WARNING - not recalculating melt energy values *******');
                        end


                        pdat(1).p = cond_dat;

                    end

                    i_max_cloud=0;
                    recalc=1;
                    if i_max_cloud==1
                        if recalc==1

                            %         M = 28.97*1.67E-27;
                            %         k = 1.38E-23;
                            %
                            %         P = nc{'P'}(time,1:ih_inds,:) + nc{'PB'}(time,1:ih_inds,:);
                            %         potemp = nc{'T'}(time,1:ih_inds,:) + 300;
                            %         T = potemp ./ ( (1e5./P).^0.286 );
                            %
                            %         rho=P.*M./k./T;
                            %         z_temp=(nc{'PH'}(time,1:ih_inds+1,:) + nc{'PHB'}(time,1:ih_inds+1,:) )./9.81;
                            %         Z = 0.5.*(z_temp(1:end-1,:)+z_temp(2:end,:));
                            %
                            %         dz_grid = Z(2:end,:,:)-Z.var(1:end-1,:,:);
                            %         cond_dat = squeeze ( sum( nc{'QCLOUD'}(time,1:ih_inds,:,:) .* rho*1000*dx_grid*1000*dy_grid.*dz_grid , 2) );

                            max_var='cloud';
                            max_var='ice no';
                            %        max_var='snow no';

                            %whether to calculate the density or not
                            switch max_var
                                case {'ice no','snow no'}
                                    potemp = nc{'T'}(time,1:ih_wrf,:,:) + 300;
                                    P = nc{'P'}(time,1:ih_wrf,:,:) + nc{'PB'}(time,1:ih_wrf,:,:);
                                    T = potemp ./ ( (1e5./P).^0.286 );
                                    rho=density(P,T);

                            end

                            %no calculate the data
                            switch max_var
                                case 'cloud'
                                    var_dat = 1000*squeeze(nc{'QCLOUD'}(time,1:ih_wrf,:,:));
                                case 'ice no'
                                    var_dat = squeeze(nc{'QNICE'}(time,1:ih_wrf,:,:)).*rho/1000;
                                case 'snow no'
                                    var_dat = squeeze(nc{'QNSNOW'}(time,1:ih_wrf,:,:)).*rho/1000;
                            end

                            %calculate the max of the data
                            if ih_wrf>1
                                max_var_dat = squeeze(max(var_dat));
                            else
                                fprintf(1,'ih_wrf=1 SO MAX OVER HEIGHT DOESN''T MAKE SENSE!');
                                break
                            end


                        else
                            fprintf(1,'\n****** WARNING - NOT recalculating MAX values *******');
                        end


                        pdat(1).p = max_var_dat;

                    end

                    iequiv=0;
                    if iequiv==1 %equivalent potential temperature
                        potemp = nc{'T'}(time,ih_wrf,:) + 300;
                        P = nc{'P'}(time,ih_wrf,:) + nc{'PB'}(time,ih_wrf,:);
                        T = potemp ./ ( (1e5./P).^0.286 );
                        qv = nc{'QVAPOR'}(time,ih_wrf,:);
                        pdat(1).p = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                    end


                    isnow_depth=0;
                    if isnow_depth==1;
                        %    snowH=nc{'SNOWH'}(:); %snow height is wrong for polar WRF Noah set up (density applied to whole snow layer)
                        %    pdat(1).p = squeeze(snowH(end,:,:)-snowH(1,:,:));

                        snow=nc{'SNOW'}(:,:);
                        pdat(1).p = squeeze(snow(end,:,:)-snow(1,:,:));
                    end

                   



                    %%%%%%%%%%%%%%%%%%%%%%%%
                    constant_height=0; %for plotting a constant height (rather than model level) by in interpolation
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    recalc=1;
                    reread=1;
                    con_height_winds=0; %set to one to also interpolate the winds - otherwise it uses the closest model level to h_wrf

                    %choose the surface to interpolate onto
                    surface='height';
                    surface='pressure';
                    surface='potemp';

                    %choose the data to interpolate
                    con_height_case = 'wind';
                    %    con_height_case = 'pressure';

                    if constant_height==1

                        if con_height_winds==0
                            switch file_type
                                case 'wrfout'
                                    Zlev=WRFUserARW(nc,'Z',time,1,1);
                                    terrain = nc{'HGT'}(time,1,:);
                                    terrain = terrain(1);
                                case 'met_em'
                                    Zlev=nc{'GHT'}(time,:,1,1);
                                    terrain = nc{'HGT_M'}(time,1,:);
                                    terrain = terrain(1);
                            end

                            ih_wrf_quiver = findheight_nearest(Zlev-terrain,h_wrf*1000);
                        end

                        if recalc==1
                            clear wrf_zint wrf_zint_u wrf_zint_v
                            if reread==1
                                %	wrf_dat = 1000*nc{'QVAPOR'}(time,:,:,:);
                                %                wrf_dat = WRFUserARW(nc,'tc',time); %load in all data for interpolation (set h_wrf in first set of settings)
                                %				wrf_dat =  squeeze(WRFUserARW(nc,'p',time));
                                %                wrf_dat=eta_get_p(nc)/100; %for wrfinput files

                                switch con_height_case
                                    case 'wind'
                                        wrf_dat = sqrt( (0.5 * ( nc{'U'}(time,:,:,1:end-1) + nc{'U'}(time,:,:,2:end) )).^2 ...
                                            + (0.5 * ( nc{'V'}(time,:,1:end-1,:) + nc{'V'}(time,:,2:end,:) ) ).^2 );

                                    case 'pressure'
                                        wrf_dat =  squeeze(WRFUserARW(nc,'p',time));
                                end

                                if con_height_winds==1
                                    u_wind_dat=0.5 * ( nc{'U'}(time,:,:,1:end-1) + nc{'U'}(time,:,:,2:end) );
                                    v_wind_dat=0.5 * ( nc{'V'}(time,:,1:end-1,:) + nc{'V'}(time,:,2:end,:) );
                                end


                                if prod(size(wrf_dat))==1  %if WRFUserARW returns the structure as wrf_dat.var
                                    wrf_dat = wrf_dat.var;
                                end

                                switch surface
                                    case 'pressure'
                                        Z_wrf = squeeze(WRFUserARW(nc,'p',time)); %if want plot at a constant pressure
                                        h_wrf = p_wrf; %
                                    case 'height'
                                        Z_wrf = squeeze(WRFUserARW(nc,'Z',time));
                                        Z_wrf.var = Z_wrf.var/1000;
                                    case 'potemp'
                                        Z_wrf.var = nc{'T'}(time,:) + 300;   %potential temperature;
                                        h_wrf = potemp_wrf;
                                end

                            end


                            %h_wrf set in first settings bit

                            igriddata=0;  %don't use gridddata as it takes too long and/or crashes
                            if igriddata==1  %if want to use the griddata3 method of interpolation instead

                                x_grid = ([1:size(lat2d(1).var,2)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
                                y_grid = ([1:size(lat2d(1).var,1)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

                                [X2D,Y2D] = MESHGRID(x_grid,y_grid);
                                X=permute(repmat(X2D,[1 1 size(Z_wrf.var,1)]) , [3 1 2]);
                                Y=permute(repmat(Y2D,[1 1 size(Z_wrf.var,1)]) , [3 1 2]);


                                %    I=[1:size(Z_wrf.var,1)];

                                %                         wrf_zint(ilat,ilon) = interp1(Z_wrf.var(:,ilat,ilon),wrf_dat(:,ilat,ilon),h_wrf*1000,[],'extrap');
                                wrf_zint = griddata3(X,Y,Z_wrf.var,wrf_dat,X2D,Y2D,h_wrf*ones(size(X2D)));
                                if con_height_winds==1
                                    wrf_zint_u(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),u_wind_dat(1:length(I),ilat,ilon),h_wrf);
                                    wrf_zint_v(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),v_wind_dat(1:length(I),ilat,ilon),h_wrf);
                                end


                            else

                                for ilat=1:size(wrf_dat,2);
                                    for ilon=1:size(wrf_dat,3);

                                        switch surface
                                            case 'potemp'
                                                [temp,I]=unique(Z_wrf.var(:,ilat,ilon));
                                                Z_wrf.var(1:length(I),ilat,ilon)=temp;
                                                wrf_dat(1:length(I),ilat,ilon)=wrf_dat(I,ilat,ilon);
                                                if con_height_winds==1
                                                    u_wind_dat(1:length(I),ilat,ilon)=u_wind_dat(I,ilat,ilon);
                                                    v_wind_dat(1:length(I),ilat,ilon)=v_wind_dat(I,ilat,ilon);
                                                end
                                            otherwise
                                                I=[1:size(Z_wrf.var,1)];
                                        end

                                        %                         wrf_zint(ilat,ilon) = interp1(Z_wrf.var(:,ilat,ilon),wrf_dat(:,ilat,ilon),h_wrf*1000,[],'extrap');
                                        wrf_zint(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),wrf_dat(1:length(I),ilat,ilon),h_wrf);
                                        if con_height_winds==1
                                            wrf_zint_u(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),u_wind_dat(1:length(I),ilat,ilon),h_wrf);
                                            wrf_zint_v(ilat,ilon) = interp1(Z_wrf.var(1:length(I),ilat,ilon),v_wind_dat(1:length(I),ilat,ilon),h_wrf);
                                        end
                                    end
                                end

                            end

                        end

                        pdat(1).p = wrf_zint;
                        if con_height_winds==1
                            u_quiver = wrf_zint_u;
                            v_quiver = wrf_zint_v;
                        end
                    end
                    %    pdat(1).p=squeeze(WRFUserARW(nc,'Z',time,0,0,ih_wrf));  %horiz_slice at first level

                    if constant_height==0
                        ih_wrf_quiver = ih_wrf;
                    end


                    if constant_height==0 | con_height_winds==0

                        %ih_wrf_quiver = 4;

                        switch file_type
                            case 'wrfout'
                                if strcmp(short_plot_name,['10m wind speed (m s^{-1}) for ' Times(time,9:10) ' Jan ' Times(time,12:16) ' UTC for ' filestr])==1
                                    u_quiver=nc{'U10'}(time,:,:);
                                    v_quiver=nc{'V10'}(time,:,:);
                                else
                                    u_quiver=0.5 * ( nc{'U'}(time,ih_wrf_quiver,:,1:end-1) + nc{'U'}(time,ih_wrf_quiver,:,2:end) );
                                    v_quiver=0.5 * ( nc{'V'}(time,ih_wrf_quiver,1:end-1,:) + nc{'V'}(time,ih_wrf_quiver,2:end,:) );
                                end
                            case 'met_em'
                                UU = nc{'UU'}(1,ih_wrf_quiver,:,:);
                                VV = nc{'VV'}(1,ih_wrf_quiver,:,:);
                                u_quiver=0.5 * ( UU(:,1:end-1) + UU(:,2:end) );
                                v_quiver=0.5 * ( VV(1:end-1,:) + VV(2:end,:) );
                        end

                    end



                    x_quiver=timesTH(1).t;
                    y_quiver=zz(1).z;


                    if(ixlim==1 & iylim==1)
                        [xinds(1),xinds(2)] = findheight(timesTH(1).t,xlims(1),xlims(2));
                        [yinds(1),yinds(2)] = findheight(timesTH(1).t,ylims(1),ylims(2));
                        xinds = [xinds(1):xinds(2)];
                        yinds = [yinds(1):yinds(2)];
                        timesTH(1).t = timesTH(1).t(xinds);
                        zz(1).z = zz(1).z(yinds);
                        ixlim=0;
                        iylim=0;
                        pdat(1).p=pdat(1).p(yinds,xinds);
                    end

