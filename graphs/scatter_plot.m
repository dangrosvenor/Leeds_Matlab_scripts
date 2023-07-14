scatter_case = 'LWC CAS';
%scatter_case = 'LWC PADS';
scatter_case = 'MODIS';
%scatter_case ='template';
%scatter_case ='Nd vs LWC difference from adiabatic';


if ~exist('savedir')
    savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';
end


%set default flags
fsize=14;
ixlims=0;
iylims=0;
ione_to_one_line=0; %whether to draw a one-to-one line or not
iadd_line=0; %flag to add normal line defined by xline and yline
icolour = 0; %flag to add a variable for the colour scaling of the scatter points
iset_cmap=0;
iclims=0;
iy_reverse=0;
ilimit_data=0;
scatter_properties_string = ['''k''']; 
iadd_nums_above=0; %flag to say whether to print numbers inside data points
ioverride_DY=1; %override the displacement of the text
DY=0; %(set to zero)
ijoin_points=0; %flag to join the points together
iplot_error_bar=0; %flag to say whether to draw horizontal error bars
i_override_vertical_coord=0;
iaxis_square=1;


f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

switch scatter_case
    case {'LWC CAS','LWC PADS'}
        eval(['time_flt=time_flt' flight_no ';']);

        if size(dat_flt,2)==15
            %for flt_19
            col_temp=6;
            col_alt=11;
            col_lat=2;
            col_lon=3;
            col_press=4;
            col_wind=9;
            col_winddir=10;
        else
            %for Feb2010 flights
            col_temp=5;
            col_alt=12;
            col_lat=2;
            col_lon=3;
            col_press=6;
            col_wind=9;
            col_winddir=10;
            col_frostpoint_hygro=7;
            col_frostpoint_humi=8;
            col_airspeed=4;
        end


        temperature_data = dat_flt(:,col_temp);        

end



switch scatter_case
    case 'template'
        
        xdat(1).x = lwc_save(2).x - lwc_save(1).x;
        ydat(1).y = Nd_save(1).x;
        
        xlabel_name = 'LWC_{ad} - LWC (g m^{-3})';
        ylabel_name = 'Nd (cm^{-3})';

        ixlims=0;
        iylims=0;

        xlims = [0 1.2];
        ylims = [0 8.5];

        icolour = 0;
        %                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
        %colour data

        colour_data_var='';

        title_name = [scatter_case];

        ione_to_one_line=0; %draw a one-to-one line
        
     case 'Nd vs LWC difference from adiabatic'
        
        xdat(1).x = lwc_save(2).x - lwc_save(1).x;
        ydat(1).y = Nd_save(1).x;
        
        xlabel_name = 'LWC_{ad} - LWC (g m^{-3})';
        ylabel_name = 'Nd (cm^{-3})';

        ixlims=0;
        iylims=0;

        xlims = [0 1.2];
        ylims = [0 8.5];

        icolour = 0;
        %                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
        %colour data

        colour_data_var='';

        title_name = [scatter_case ' ' profile_set];

        ione_to_one_line=0; %draw a one-to-one line    
                
                
    case 'MODIS'
        plot_type='CF35 vs CF06'; %Cloud Fraction from MOD35 vs. that from MOD06
        plot_type='Nd with W calculated vs with MODIS W';
        plot_type='W calculated vs from MOD06';
        plot_type='W histogram vs from MOD06';
%        plot_type='Nd from mean tau,reff vs. mean Nd from joint tau,reff histogram';
%         plot_type='W vs. no. pixels';        
         plot_type='Tau vs. no. pixels'
%        plot_type='Tau vs. reff';  %for different W
%        plot_type='W vs. tau';
         plot_type='N vs.'; 
%         plot_type='std(N) vs.';
         plot_type='Nd CDP vs. Nd MODIS'; 

        switch plot_type
            case 'Nd CDP vs. Nd MODIS'
                scatter_properties_string = ['100,''b'',''filled'''];
                
             %y-axis   
                Nd_case = 'Mean tau mean re';
                Nd_case = 'Individual pixels'; highlight_case = 'MODIS';
%                Nd_case = 'log Individual pixels'; highlight_case = 'MODIS';                
%                Nd_case = 'CDP';   highlight_case = 'CDP';

highlight_case='';
                
             %x-axis
                Puijo_instrument = 'CDP';
%                Puijo_instrument = 'CDP NaN match';
%                Puijo_instrument = 'DMPS';   %run read_Nacc_data
%                Puijo_instrument = 'log DMPS';   %run read_Nacc_data                
%                Puijo_instrument = 'DMPS Nd';
%                Puijo_instrument = 'DMPS Nd NaN match';

           DMPS_type = 'Irshad';
           DMPS_type = 'Irshad new';           
%           DMPS_type = 'marine';
%           DMPS_type = 'continental';
%           DMPS_type = 'all';
%           DMPS_type = 'Dan';


                iplot_error_bar=1;
              
                
                re_wavelength = '1.6';
                re_wavelength = '2.1';
                re_wavelength = '3.7';
%                re_wavelength = 're 2.1 *0.8';

                std_rel_abs = 'abs';
                
                ifilter_std=1;
                std_thresh = 50; %relative std (%)                
                std_thresh = 150; %absolute std (cm3)
%                std_thresh = 140; %absolute std (cm3)                
                
                switch Nd_case
                    case 'Mean tau mean re'
                        switch re_wavelength
                            case '2.1'
                                y_axis_vals = 'Nd from grid vals timeseries3';
                            case '1.6'
                                y_axis_vals = 'Nd_{1.6} from grid vals timeseries3';
                            case '3.7'
                                y_axis_vals = 'Nd_{3.7} from grid vals timeseries3';
                            case 're 2.1 *0.8'
                                y_axis_vals = 'Nd from grid vals re*0.8 timeseries3';                                                                                     
                        end
                        
                    case 'Individual pixels'
                        switch re_wavelength
                            case '2.1'
                                y_axis_vals = 'Nd from individual pixels timeseries3';
                            case '1.6'
                                y_axis_vals = 'Nd from individual pixels 1.6\mum timeseries3';
                            case '3.7'
                                y_axis_vals = 'Nd from individual pixels 3.7\mum timeseries3';
                        end
                        
                    case 'log Individual pixels'
                        switch re_wavelength
                            case '2.1'
                                y_axis_vals = 'log Nd from individual pixels timeseries3';
                            case '1.6'
                                y_axis_vals = 'Nd from individual pixels 1.6\mum timeseries3';
                            case '3.7'
                                y_axis_vals = 'Nd from individual pixels 3.7\mum timeseries3';
                        end    
                        
   
                    case 'CDP'
                         y_axis_vals = 'CDP Nd';
                        
                end
                
                ioverride_pdf_varchoose=1;
                
                switch Puijo_instrument
                    case 'CDP'
                        x_axis_vals = 'CDP Nd';
                    case 'CDP NaN match'
                        x_axis_vals = 'CDP Nd NaN match';
                    case 'DMPS'
                        x_axis_vals = 'DMPS Nacc';
                    case 'log DMPS'
                        x_axis_vals = 'log DMPS Nacc';    
                    case 'DMPS Nd'
                        x_axis_vals = 'DMPS Nd estimate';
                    case 'DMPS Nd NaN match'
                        x_axis_vals = 'DMPS Nd estimate CDP NaN match';
                end

                pdf2D_plot_commands
                
                yerror = sqrt(std_dev.^2 + instrument_err.^2);
                Combined_Nd_stddev_and_instrument_error.timeseries3 = yerror;
                

                if ifilter_std==1  
                    switch std_rel_abs
                        case 'abs'
                            istd=find(squeeze(yerror(ihtot)) <std_thresh);
                        case 'rel'
                            istd=find(100*squeeze(yerror(ihtot))./Y <std_thresh);
                    end
                else
                    istd=1:length(X);
                end
                
                
                xdat(1).x = X(istd);
                ydat(1).y = Y(istd);
                
                error_barH(1).dat=stdNd_CDP.timeseries3(ihtot(istd));
                error_barV(1).dat=yerror(ihtot(istd));

                
                
                
                ixlims=1;
                iylims=1;

                maxXY = max(xdat(1).x);
                maxXY = max(cat(1,ydat(1).y,maxXY));
                
                xlims = [0 maxXY+100];
                ylims = [0 maxXY+100];
                
                xlims = [0 1000];
                ylims = [0 1000];
                
%                xlims = [0 1600];
%                ylims = [0 500];
                
%                xlims = [0 450];
%                ylims = [0 450];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of
%                code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                if strcmp(thresh_str,'xxx')~=1
                    thresh_str2=[' for ' thresh_str];
                else
                    thresh_str2='';
                end
                title_name = [plot_type colour_data_var thresh_str2];
                
                ione_to_one_line=1; %draw a one-to-one line
                i_override_vertical_coord=0;
                
                xlabel_name = xlabelstr;
                ylabel_name = ylabelstr;
                
                ihighlightX = 0;
                ihighlightY = 0;  


                
                switch highlight_case
                    case 'MODIS'
                        clear highlight_pointsX highlight_pointsY
                        highlight_pointsY(1) = 176.6; %MOIDIS y-values
                        highlight_pointsY(2) = 139.1;
                        ihighlightY = 1;

                    case 'CDP'
                        highlight_pointsY(1) = 370;
                        highlight_pointsY(2) = 399.4; %for CDP
                        ihighlightY = 1;
                end

                
                
                
                
             case 'std(N) vs.'  
  
                
                thresh=150;                
                thresh=0.3;
                thresh=0;
                
                %this script does the calculations
                N_H_calc_std_histo
                
                
                cf = Cloud_Fraction_Liquid.data;
                WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                tau = Cloud_Optical_Thickness_Liquid_Mean.data;
%                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres
                sangle = Scattering_Angle_Mean.data;
                

                
                
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
%                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
                 ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data

                xdat(1).x = Nstd_norm(ihtot);               
                xlabel_name = 'Std. dev of N_d divided by mean N_d)';
                
Nvs_case='WMOD';
Nvs_case='CF';
%Nvs_case='Scattering angle';

plot_type=[plot_type ' ' Nvs_case];

                switch Nvs_case
                    case 'WMOD'
                        ylabel_name = 'W times five sixths (kg m^{-2})';
                        ydat(1).y = WMOD(ihtot); %MOD06
                    case 'CF'
                        ylabel_name = 'Cloud Fraction';
                        ydat(1).y = cf(ihtot); %MOD06 CF
                    case 'Scattering angle'
                        ylabel_name = 'Scattering angle';
                        ydat(1).y = sangle(ihtot); 
                end
                
                ixlims=0;
                iylims=0;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                if strcmp(thresh_str,'xxx')~=1
                    thresh_str2=[' for ' thresh_str '.GTE.' num2str(thresh)];
                else
                    thresh_str2='';
                end
                title_name = [plot_type colour_data_var thresh_str2 ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=0; %draw a one-to-one line
                i_override_vertical_coord=0;
                
            case 'N vs.'  
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                MODIS_N_H_calc %runs the routine to calculate Nd    
                sangle = Scattering_Angle_Mean.data;
                
                timeseries_case=''; %for timeseries
                timeseries_case='3'; %for timeseries3 
                
                LAT=eval(['Cloud_Optical_Thickness_Liquid_Mean.timeseries' timeseries_case '_LAT']);
                LON=eval(['Cloud_Optical_Thickness_Liquid_Mean.timeseries' timeseries_case '_LON']);


                LAT_val = [-64.5000];
                LAT_val = [-60.5000];
                LAT_val = [-59.5000];                
%                LAT_val = [-89.5:1:89.5];                
                %LAT_val = [-12.5000];

                ilat=findheight_nearest(LAT,LAT_val);
                
                LON_val = [-179.5:179.5];
                
                ilon=findheight_nearest(LON,LON_val);                

                
                
                x_axis_vals = 'Nd from grid vals timeseries';
                x_axis_vals = 'Nd from grid vals timeseries3';                
%                x_axis_vals = 'Nd from histogram vals timeseries';
%                x_axis_vals = 'Nd from histogram vals timeseries3';                
%                x_axis_vals = 'Std. dev of Nd from histogram timeseries';
%                x_axis_vals = 'Std. dev of Nd from histogram timeseries3';                
                
                
                Nvs_case='WMOD';
                Nvs_case='CF';
                %Nvs_case='Scattering angle';
                Nvs_case='SZA std dev timeseries';
                Nvs_case='Mean SZA timeseries';
                Nvs_case='Mean SZA timeseries3';                
                %Nvs_case='SZA std dev timeseries2';
%                Nvs_case='SZA std dev timeseries3';                
%                Nvs_case = 'Nd from histogram vals timeseries';
%                Nvs_case = 'Nd from histogram vals timeseries3';

                
                thresh=150;                
                thresh=0.3;
                thresh=0;
                
                cf = Cloud_Fraction_Liquid.data;
                WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                tau = Cloud_Optical_Thickness_Liquid_Mean.data;
%                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres



                cf_time=Cloud_Fraction_Liquid.timeseries(ilat,:);
                NP_time=Cloud_Fraction_Liquid_Pixel_Counts.timeseries(ilat,:)./cf_time;
                
                tau_time = Cloud_Optical_Thickness_Liquid_Mean.timeseries(ilat,:);
                reff_time = Cloud_Effective_Radius_Liquid_Mean.timeseries(ilat,:)*1e-6; %convert to metres
                Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)

                [N_time,H,W,k,Q,cw]=MODIS_N_H_func(tau_time,reff_time,Wflag,0);
                
                
                %timeseries3
                cf_time3=Cloud_Fraction_Liquid.timeseries3(ilat,:,:);
                NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3(ilat,:,:)./cf_time3;
                
                tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,:,:);
                reff_time3 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,:,:)*1e-6; %convert to metres
                Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)

                [N_time3,H,W,k,Q,cw]=MODIS_N_H_func(tau_time3,reff_time3,Wflag,0);
                
                
               


                switch x_axis_vals
                    case 'Nd from grid vals'
                        xlabel_name = 'N_d (cm^{-3})';
                        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
                        xdat(1).x = N;

                    case 'Nd from grid vals timeseries'
                        xlabel_name = 'N_d (cm^{-3}) from grid values timeseries';
                        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data
                        xdat(1).x = N_time;

                    case 'Nd from grid vals timeseries3'
                        xlabel_name = 'N_d (cm^{-3}) from grid values timeseries3';
                        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(1,:)))]; thresh_str='xxx'; %all the data
                        xdat(1).x = N_time3;

                        
                    case 'Nd from histogram vals timeseries'
                        xlabel_name = 'N_d (cm^{-3}) from histogram timeseries';
                        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data
                        
                        [histo_output]=...
                            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,ihtot) );
                        xdat(1).x=histo_output.N_histo_mean;
                        %                        xdat(1).x=histo_output.N_histo_std;
                        
                    case 'Nd from histogram vals timeseries3'
                        xlabel_name = 'N_d (cm^{-3}) from histogram timeseries3';
                        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(1,:)))]; thresh_str='xxx'; %all the data
                        
%                        [histo_output]=...
%                            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );
%                        xdat(1).x=histo_output.N_histo_mean;
                        
xdat(1).x = Nd_timeseries.mean(ilat,ilon,:);
                        
                        
                    case 'Std. dev of Nd from histogram timeseries'
                        
                        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data
                        
                        [histo_output]=...
                            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,ihtot) );

                        xlabel_name = 'Std dev N_d (cm^{-3}) from histogram timeseries';
                        xdat(1).x=histo_output.N_histo_std;

                        xlabel_name = 'Normalised std dev N_d from histogram timeseries';
                        xdat(1).x=histo_output.N_std_norm;
                        
                     case 'Std. dev of Nd from histogram timeseries3'
                        
                        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:)))]; thresh_str='xxx'; %all the data
                        
                        [histo_output]=...
                            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );

                        xlabel_name = 'Std dev N_d (cm^{-3}) from histogram timeseries';
                        xdat(1).x=histo_output.N_histo_std;

                        xlabel_name = 'Normalised std dev N_d from histogram timeseries';
                        xdat(1).x=histo_output.N_std_norm;    


                end


                                
%%%%%  choosing only certain points (lat-lon grid)                 
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
%                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
%                 ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                
                
%%%%%  choosing only certain points (timeseries)                
%                ihtot = find(NP_time>50 & cf_time>0.8); thresh_str='NP.GT.50 AND CF.GT.0.8';
                
%%%%%  choosing only certain points (timeseries3)                
                ihtot = find(NP_time3>50 & cf_time3>0.8); thresh_str='NP.GT.50 AND CF.GT.0.8';                
        
        
        xdat(1).x=xdat(1).x(ihtot);

                                         



plot_type=[x_axis_vals ' vs ' Nvs_case];

extra_title_info = [' for day ' modis_day_str ', Y' modis_year_str];

if length(ilat)==1
    LAT_str=num2str(LAT(ilat));
else
    LAT_str=[num2str(LAT(ilat(1))) ' to ' num2str(LAT(ilat(end)))];
end

if length(ilon)==1
    LON_str=num2str(LON(ilon));
else
    LON_str=[num2str(LON(ilon(1))) ' to ' num2str(LON(ilon(end)))];
end


                switch Nvs_case
                    case 'WMOD'
                        ylabel_name = 'W times five sixths (kg m^{-2})';
                        ydat(1).y = WMOD(:); %MOD06
                        
                    case 'CF'
                        ylabel_name = 'Cloud Fraction';
                        ydat(1).y = cf(:); %MOD06 CF

                    case 'Scattering angle'
                        ylabel_name = 'Scattering angle (degrees)';
                        ydat(1).y = sangle(:);

                    case 'SZA std dev timeseries'
                        ylabel_name = 'Std dev of SZA';
                        ydat(1).y = Solar_Zenith_Standard_Deviation.timeseries(ilat,:);  
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
                    case 'Mean SZA timeseries'
                        ylabel_name = 'Mean SZA';
                        ydat(1).y = Solar_Zenith_Mean.timeseries(ilat,:);   
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
                    case 'Mean SZA timeseries3'
                        ylabel_name = 'Mean SZA';
                        ydat(1).y = Solar_Zenith_Mean.timeseries3(ilat,ilon,:,:);   
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
                        
                    case 'SZA std dev timeseries2'
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];                                               
                        ylabel_name = 'Std dev of SZA';
                        ydat(1).y = Solar_Zenith_Standard_Deviation.timeseries(ilat,:);  
%                        xdat(1).x = N_time(ihtot);

                    case 'SZA std dev timeseries3'
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];                                               
                        ylabel_name = 'Std dev of SZA';
                        ydat(1).y = Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:,:);  
%                        xdat(1).x = N_time(ihtot);

                    case 'Nd from histogram vals timeseries'
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];                                               
                        ylabel_name = 'N_d from histo timeseries (cm^{-3})';
                        
                        [histo_output]=...
                            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,:) );
                        ydat(1).y=histo_output.N_histo_mean;                       
                        
                     case 'Nd from histogram vals timeseries3'
                        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];                                               
                        ylabel_name = 'N_d from histo timeseries3 (cm^{-3})';
                        
                        [histo_output]=...
                            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );
                        ydat(1).y=histo_output.N_histo_mean;                       
                        
                        
                end
                
                %%% only select the required points
                ydat(1).y = ydat(1).y(ihtot);
                
                ixlims=0;
                iylims=0;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                if strcmp(thresh_str,'xxx')~=1
                    thresh_str2=[' for ' thresh_str];
                else
                    thresh_str2='';
                end
                title_name = [plot_type colour_data_var thresh_str2 extra_title_info];                
                
                ione_to_one_line=0; %draw a one-to-one line
                i_override_vertical_coord=0;
                
                
                
             case 'W vs. tau'  
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                
                thresh=150;                
                thresh=0.3;
                thresh=0;
                
%                cf = Cloud_Fraction_Liquid.data;
                WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                tau = Cloud_Optical_Thickness_Liquid_Mean.data;
%                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres
                
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
%                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
                 ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data

                xdat(1).x = tau(ihtot);               
                xlabel_name = 'Optical Depth';
                
%                xdat(1).x = reff(ihtot);               
%                xlabel_name = 'Reff (\mum)';
                
                
                
                ylabel_name = 'W times five sixths (kg m^{-2})';
                ydat(1).y = WMOD(ihtot); %MOD06 
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                if strcmp(thresh_str,'xxx')~=1
                    thresh_str2=[' for ' thresh_str '.GTE.' num2str(thresh)];
                else
                    thresh_str2='';
                end
                title_name = [plot_type colour_data_var thresh_str2 ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=0; %draw a one-to-one line
                i_override_vertical_coord=0;
                
                
                
             case 'Tau vs. reff'  
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                
                thresh=150;                
                thresh=0.5;
                
%                cf = Cloud_Fraction_Liquid.data;
                WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                tau = Cloud_Optical_Thickness_Liquid_Mean.data;
                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres
                
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
%                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
                 ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths';
                
                xdat(1).x = tau(ihtot);               
                xlabel_name = 'Optical Depth';
                                                   
                ylabel_name = 'Reff (\mum)';
                ydat(1).y = reff(ihtot); %MOD06 
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                title_name = [plot_type colour_data_var ' for ' thresh_str '.GTE.' num2str(thresh) ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=0; %draw a one-to-one line
                i_override_vertical_coord=0;
                
            case 'Tau vs. no. pixels'
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                MODIS_N_H_calc
                
                thresh=150;                
                thresh=0.5;
                                
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels GT ';
%                %only plot for data with more a threshold no. pixels in total
%                ihtot = find(cf>=thresh);  thresh_str='CF GTE '; %only plot for data with more a threshold value (of CF in this case)
                 ihtot = find(WMOD>=thresh);  thresh_str='W times five-sixths GTE ';
                 ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
                
                xdat(1).x = tau(ihtot);               
                xlabel_name = 'Optical Depth';
                                                   
                ylabel_name = 'No. pixels';
                ydat(1).y = totN(ihtot); %MOD06 
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                title_name = [plot_type colour_data_var ' for ' thresh_str num2str(thresh) ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=0; %draw a one-to-one line
                i_override_vertical_coord=0;
                
                
                
            case 'W vs. no. pixels'    
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                
                thresh=150;                
                thresh=0.0;
                
                cf = Cloud_Fraction_Liquid.data;
                
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
                ihtot = [1:prod(size(N_histo_mean))]; %all the data
                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
                
                xdat(1).x = totN(ihtot); 
              
                xlabel_name = 'No. pixels';
                                                   
                WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                ylabel_name = 'W from MOD06 times five-sixths (kg m^{-2})';
                fact= (5/9) / (2/3);  % fact=5/6  convert the MODIS 2/3 formulua to 5/9
%                fact=1;
                ydat(1).y = fact * WMOD(ihtot); %MOD06 
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                title_name = [plot_type colour_data_var ' for ' thresh_str '.GTE.' num2str(thresh) ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=0; %draw a one-to-one line
                i_override_vertical_coord=0;
                
            case 'W histogram vs from MOD06'    
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                MODIS_N_H_calc
                
                thresh=150;                
                thresh=0.0;
                thresh=50;
                
                cf = Cloud_Fraction_Liquid.data;
                
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
                ihtot = find(cf>=thresh);  thresh_str='CF GTE '; %only plot for data with more a threshold value (of CF in this case)
                ihtot = find(tau<=thresh);  thresh_str='tau LTE '; %only plot for data with more a threshold value (of CF in this case)
                
                xdat(1).x = W_histo_mean(ihtot); 
              
                xlabel_name = 'W from histo five ninths formula (kg m^{-2})';
                                                   
                WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                ylabel_name = 'W from MOD06 times five-sixths (kg m^{-2})';
%                ylabel_name = 'W from MOD06 using mean reff and tau (kg m^{-2})';
                
                fact= (5/9) / (2/3);  % fact=5/6  convert the MODIS 2/3 formulua to 5/9
%                fact=1;
                ydat(1).y = fact * WMOD(ihtot); %MOD06 
%                 ydat(1).y = W(ihtot);
                
                
                ixlims=1;
                iylims=1;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                if icolour==0                
                    colour_data_var='';
                else                    
                    colour_data_var=' with tot pixels';
                    colour_data=totN(ihtot);
                    
                end
                
%                title_name = [plot_type ' for flight ' flight_no ' with '
%                colour_data_var];
                title_name = [plot_type colour_data_var ' for ' thresh_str num2str(thresh) ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=1; %draw a one-to-one line
                i_override_vertical_coord=0;
                
            case 'W calculated vs from MOD06'    
                cf = Cloud_Fraction_Liquid.data;
                tau = Cloud_Optical_Thickness_Liquid_Mean.data;
                reff = Cloud_Effective_Radius_Liquid_Mean.data*1e-6; %convert to metres
                rhoL=1000; %water density kg/m3
                WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
                Wcalc=2/3*rhoL.*tau.*reff;
%                Wcalc=5/9*rhoL.*tau.*reff;

                xdat(1).x = Wcalc(:); %calculate W                
                xlabel_name = 'W from two thirds formula w. no cf scaling (kg m^{-2})';
                                                                              
                ylabel_name = 'W from MOD06 (kg m^{-2})';
                ydat(1).y = WMOD(:); %MOD06 
                
                
                ixlims=1;
                iylims=1;

                xlims = [0 2];
                ylims = [0 2];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                colour_data_var='';  

%                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];
                title_name = [plot_type ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=1; %draw a one-to-one line
                i_override_vertical_coord=0;
                        
            case 'Nd with W calculated vs with MODIS W'
                set_MODIS_NH_flags=1;
                Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
                MODIS_N_H_calc %runs the routine to calculate Nd                
                
                xdat(1).x = N(:); %MOD35                               
                xlabel_name = 'N_d, W calculated';
                                                              
                set_MODIS_NH_flags=1; %gets reset every time
                Wflag='MODIS'; %use the MODIS LWP
                MODIS_N_H_calc %runs the routine to calculate Nd
                
                ylabel_name = 'N_d, MODIS W';
                ydat(1).y = N(:); %MOD06 
                
                
                 ixlims=0;
                iylims=0;

                xlims = [0 1.2];
                ylims = [0 8.5];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                colour_data_var='';  

%                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];
                title_name = [plot_type ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=1; %draw a one-to-one line
                i_override_vertical_coord=0;
                
            case 'Nd from mean tau,reff vs. mean Nd from joint tau,reff histogram'
                
                set_MODIS_NH_flags=1; %gets reset every time
                Wflag='calc'; %use the calculated LWP (default)
                N_H_calc_histo %runs the routine to calculate Nd from the joint histogram
                
                ihtot = find(totN>20);  %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
                
                xdat(1).x = N_histo_mean(ihtot); %mean value from binned tau, reff for each lat,lon                              
                xlabel_name = 'N_d, histogram (cm^{-3})';
                                                              
                set_MODIS_NH_flags=1; %gets reset every time
                Wflag='calc'; %use the calculated LWP
                MODIS_N_H_calc %runs the routine to calculate Nd (using mean tau & reff)
                
                ylabel_name = 'N_d, average tau,reff (cm^{-3})';
                ydat(1).y = N(ihtot); %MOD06 
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 1.2];
                ylims = [0 8.5];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                colour_data_var='';  

%                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];
                title_name = [plot_type ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=1; %draw a one-to-one line
                i_override_vertical_coord=0;
                
    
            case 'CF35 vs CF06'; %Cloud Fraction from MOD35 vs. that from MOD06
        
                i_override_vertical_coord=0;
                
                xdat(1).x = Cloud_Fraction_Day_Mean.data(:); %MOD35
                ydat(1).y = Cloud_Fraction_Combined.data(:); %MOD06                               
                
                xlabel_name = 'Cloud Fraction MOD35';
                ylabel_name = 'Cloud Fraction MOD06';                                                            

                ixlims=0;
                iylims=0;

                xlims = [0 1.2];
                ylims = [0 8.5];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                colour_data_var='';  

%                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];
                title_name = [plot_type ' for day ' modis_day_str ', Y' modis_year_str];
                
                ione_to_one_line=1; %draw a one-to-one line
        end
                

    case 'LWC CAS'

        %        eval(['dat_flt = dat_flt' flight_no ';']);  %put the data for the required flight here
        %        eval(['time_flt = time_flt' flight_no ';']);

        plot_type='CAS LWC vs hotwire LWC';

        %         plot_type='CAS LWC vs conc.';
 %                plot_type='CAS LWC vs vertical windspeed';        
%                 plot_type='Hotwire LWC vs vertical windspeed';                         
        %         plot_type='CAS LWC vs conc. cutoff';
        %         plot_type='Hotwire LWC vs conc.';
        %         plot_type='Hotwire LWC vs mode';
        %         plot_type='Hotwire LWC vs mean';
        %        plot_type='CAS LWC vs mode';
        %        plot_type='CAS LWC vs mean';
        %        plot_type='Mode vs number';
 
        plot_type='CAS LWC';            
        plot_type='Hotwire LWC';       

%                plot_type='CAS number vs height';
%        plot_type='Aircraft altitude';   
%        plot_type='Ice mean diameter (\mum)';
        plot_type='Ice number (L^{-1})';
        plot_type='Ice number binned (L^{-1})';        
%        plot_type='Ice mass (mg m^{-3})';
%        plot_type='Wind direction (degrees) vs height';
%        plot_type='Vapour mixing ratio (kg kg^{-1}) (Hygrometer)';
%        plot_type='Vapour mixing ratio (kg kg^{-1}) (Humicap)';    
%         plot_type='Equivalent potential temperature (K) (Humicap)';    
%         plot_type='Equivalent potential temperature (K) (Saturation when LWC)';
%         plot_type='Estimated displacement (m) vs height';
%         plot_type='Vertical windspeed (m s^{-1})';
%         plot_type='CAS hotwire matches';
%          plot_type='Number of droplets vs. mode diameter';


if ~exist('ichoose_scatter_vert_coord') | ichoose_scatter_vert_coord==0
%%%%%%%%%   *** choose the vertical coordinate ***  %%%%%%%%%%%


i_override_vertical_coord=1; %flag to say to ignore what is written in individual bits
%and use the common bit of code for setting the vertical coordinate

        vertical_coord='pressure';
        vertical_coord='temperature';   
%        vertical_coord='height'; 
%        vertical_coord='potemp';  
%        vertical_coord='equiv potemp fp';  
%        vertical_coord='equiv potemp humi';  
%        vertical_coord='equiv potemp sat';  %changed to use qsat          
%         vertical_coord='vertical velocity';      
%         vertical_coord='lwc';      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

else
    clear ichoose_scatter_vert_coord
end

        if exist('dat_flt')



            times=[14.65 16.097];
            times=[14.4 14.95];
            times=[14.625 14.71];
            times=[20.3791 20.4107]; %flight 102 - adiabatic wave crest A1 - 70-80km
            times=[20.5122 20.5344]; %flight 102 - adiabatic wave crest A2 - 108-115 km
            times=[20+38/60 20+38.75/60];
            times=[20.6312 20.6514];         
            times=[20.8143 20.8646];   %flight 102 192-208 km];   
            
            times=[19.6861 20.3949]; %flight 102 - ascent
%            times=[21.6772 22.0281]; %flight 102 - descent
            times=[21 21+20/60]; %flight 102 - period of high ice concs
            times=[13+13/60 14+0/60]; %flight 100 - initial ascent
            times=[14+0/60 14+20/60]; %flight 100 - mid-level leg 1         
            times=[14+20/60 14+37/60]; %flight 100 - mid-level leg 2
            times=[14+37/60 14+53/60]; %flight 100 -             
            
            times=[13+17/60 15+51/60];
%            times=[13+40/60 13+50/60];            
            times=[20+20/60 20+40/60]; %flight 102 20:20 - 20:40
            times=[18+55/60 19+25/60]; %flight 104
            times=[19+8.5/60 19+22/60]; %flight 104     
            times=[19+8.5/60 19+15/60]; %flight 104    
            times=[19+3/60 19+17/60]; %flight 104                
%            times=[];
            
            i_times_cut_out=0; %flag to say that times are the ones we want to cut out rather than
            %keep
            
            %Note on CIP_time_Jonny and CIP_time_Jonny2 - the former (1) is
            %86401 elements long, whereas the latter (2) is 86400 elements
            %long. (1) is corresponds to the PSDs, (2) is for ice_no_Jonny,
            %etc. This is some kind of bug in Jonny's code
            
            if length(times)==0
                inds=1:length(dat_flt(:,1));
                indsCAS=1:length(CAS_time_all);
                if exist('CIP_time_Jonny2')
                    indsJonny2 = find(CIP_time_Jonny2/3600>=time_flt(1) & CIP_time_Jonny2/3600<=time_flt(end));                
                end
                if exist('CIP_time_Jonny')
                    indsJonny = find(CIP_time_Jonny/3600>=time_flt(1) & CIP_time_Jonny/3600<=time_flt(end));                
                end
            else
                [time_0,time_1]=findheight_nearest(dat_flt(:,1)/1000/3600,times(1),times(2));
                inds=time_0:time_1;
                indsCAS = find(CAS_time_all>=dat_flt(time_0,1)/1e3 & CAS_time_all<=dat_flt(time_1,1)/1e3);
                if exist('CIP_time_Jonny2')
                    indsJonny2 = find(CIP_time_Jonny2>=dat_flt(time_0,1)/1e3 & CIP_time_Jonny2<=dat_flt(time_1,1)/1e3);
                end
                if exist('CIP_time_Jonny')
                    indsJonny = find(CIP_time_Jonny>=dat_flt(time_0,1)/1e3 & CIP_time_Jonny<=dat_flt(time_1,1)/1e3);
                end
                
                if i_times_cut_out==1
                    indsCAS2=1:length(CAS_time_all);
                    indsCAS2(indsCAS)=[];
                    indsCAS=indsCAS2;
                    
                    if exist('CIP_time_Jonny2')
                        indsJonny = 1:length(CIP_time_Jonny2);
                        indsJonny(indsJonny2)=[];
                        indsJonny2=indsJonny;
                    end
                    
                end
            end
            

                
            
            time_base_scatter=CAS_time_all(indsCAS);



            %%% these are set in plotTimeHeightVap3.m - might want to override here
            % or to make sure they are the same if requried
            if ~exist('cut_off_size')
                cut_off_size=1; %size (microns) below which to ignore counts for particle concentration
            end
            if ~exist('air_speed_type')
                %air_speed_type = 'aircraft';
                %air_speed_type = 'constant 60m/s';
                air_speed_type = 'CIP probe';
            end

            cut_off_size=1; %cut off size for CAS_total_number in microns

            air_speed_type = 'constant';
            air_speed_type = 'CIP probe';
            instrument = 'CAS'; %to determine the returned airspeed
            airspeed_constant=10; %not used unless air_speed_type set to 'constant'
            %%%    ------------------------------------   %%%

            %get the sample volume and total concentrations, plus air speed if required.
%             [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
%                 ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
%                 ,CAS_mean_diameter,LWC_dist_cas_cutoff]...
%                 =cas_sample_volume_and_stats(dat_flt,CAS_time_all...
%                 ,CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type...
%                 ,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes);
            

            CAS_LWC_cut_off_sizes=[10 30]; %upper cut off for CAS LWC (microns)
            CAS_LWC_cut_off_sizes=[10 20]; %upper cut off for CAS LWC (microns)
            CAS_LWC_cut_off_sizes=[2 50]; %upper cut off for CAS LWC (microns)            
%            CAS_LWC_cut_off_sizes=[0 50]; %upper cut off for CAS LWC (microns)

       
            if exist('iapply_ratios') & iapply_ratios==1
                clear iapply_ratios
                disp('*** WARNING - applying ratios ***');
                ratio_array = repmat([1; mean_ratio],[1 size(CAS_counts_all,1)])';
                ratio_array = repmat([1; mean_ratio_1_3],[1 size(CAS_counts_all,1)])';                
%                ratio_array = repmat([1; max_ratio],[1 size(CAS_counts_all,1)])';               
%                ratio_array = repmat([1; test_ratio],[1 size(CAS_counts_all,1)])';  
%                ratio_array = repmat([1; test_ratio4],[1 size(CAS_counts_all,1)])';  
%                ratio_array = repmat([1; test_ratio_10_15],[1 size(CAS_counts_all,1)])';                  
               ratio_array = repmat([1; test_ratio_10_15_2],[1 size(CAS_counts_all,1)])';                  
               
                CAS_counts_scatter=CAS_counts_all./ratio_array;
                iapplied_ratios=1;
            else
                CAS_counts_scatter=CAS_counts_all;
                iapplied_ratios=0;                
            end
            
            if exist('ishift_bins') & ishift_bins==1
                clear ishift_bins
                disp('*** WARNING - shifting bins ***');
                
%                shift=-1; %how much to shift by
                clear CAS_bins_scatter                
%                CAS_bins_scatter = CAS_bins/1.1;
                                
                CAS_bins_scatter(1)=CAS_bins(1)-0.1; 

                CAS_bins_scatter(2:length(CAS_bins))=CAS_bins(1:end-1);
                CAS_bins_scatter=CAS_bins_scatter';
                ishifted_bins=1;
            else
                CAS_bins_scatter=CAS_bins;
                ishifted_bins=0;                
            end
            
                     
            
             time_timeseries = CAS_time_all;
            
            [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
                ,CAS_total_number_cutoff ...                                    
                ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
                        ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
                =cas_sample_volume_and_stats2...
                (dat_flt,time_timeseries,...
                CAS_bins_scatter,CAS_counts_scatter,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,airspeed_constant); 

            
     alt=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                 
                 

        end

        switch plot_type
            case 'Number of droplets vs. mode diameter'
                i_override_vertical_coord=0;
                
                CAS_LWC = LWC_dist_cas_cutoff';

%                 %hotwire LWC
%                 hotwire_LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
%                 hotwire_LWC(hotwire_LWC<-0.05)=NaN;
%                 hotwire_LWC(hotwire_LWC>2.5)=NaN; %is an LWC point at 4 g/m3 that must be wrong
%                 
%                 %only plot and calculate for certain conditions
                 max_num=maxALL(CAS_total_number_cutoff);
                 inds_ratio = find(CAS_total_number_cutoff>=max_num*0.0);
%                 ratio = CAS_LWC(inds_ratio)./hotwire_LWC(inds_ratio);

%                 
%                 cas_lwc_xdat = CAS_mode_diameter(inds_ratio);
                
                xdat(1).x = CAS_mode_diameter(inds_ratio);
                ydat(1).y = CAS_total_number_cutoff(inds_ratio);
                
                %example of binning the data
                [meanvals_ratio,med_vals,maxvals_ratio,max_inds_ratio,mid_points,nvals]=bin_data(xdat(1).x,ydat(1).y,CAS_bins);
                
                
                inds_ratio = find(CAS_total_number_cutoff>=max_num*0.25);                
                xdat(1).x = CAS_mode_diameter(inds_ratio);
                ydat(1).y = CAS_total_number_cutoff(inds_ratio);               
                [meanvals_ratio2,med_vals2,maxvals_ratio2,max_inds_ratio2,mid_points,nvals2]=bin_data(xdat(1).x,ydat(1).y,CAS_bins);

                inds_ratio = find(CAS_total_number_cutoff>=max_num*0.5);                
                xdat(1).x = CAS_mode_diameter(inds_ratio);
                ydat(1).y = CAS_total_number_cutoff(inds_ratio);               
                [meanvals_ratio3,med_vals3,maxvals_ratio3,max_inds_ratio3,mid_points,nvals3]=bin_data(xdat(1).x,ydat(1).y,CAS_bins);
                                                
                inds_ratio = find(CAS_total_number_cutoff>=max_num*0.65);                
                xdat(1).x = CAS_mode_diameter(inds_ratio);
                ydat(1).y = CAS_total_number_cutoff(inds_ratio);               
                [meanvals_ratio4,med_vals4,maxvals_ratio4,max_inds_ratio4,mid_points,nvals4]=bin_data(xdat(1).x,ydat(1).y,CAS_bins);

                inds_ratio = find(CAS_total_number_cutoff>=max_num*0.75);                
                xdat(1).x = CAS_mode_diameter(inds_ratio);
                ydat(1).y = CAS_total_number_cutoff(inds_ratio);               
                [meanvals_ratio5,med_vals5,maxvals_ratio5,max_inds_ratio5,mid_points,nvals5]=bin_data(xdat(1).x,ydat(1).y,CAS_bins);

                inds_ratio = find(CAS_total_number_cutoff>=max_num*0.85);                
                xdat(1).x = CAS_mode_diameter(inds_ratio);
                ydat(1).y = CAS_total_number_cutoff(inds_ratio);               
                [meanvals_ratio6,med_vals6,maxvals_ratio6,max_inds_ratio6,mid_points,nvals6]=bin_data(xdat(1).x,ydat(1).y,CAS_bins);

                
                
                
                
                xlabel_name = 'CAS dN/dlogD mode size (\mum)';
                ylabel_name = ['CAS total no. ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum'];
                                


                              
                colour_data_var='';
                
                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];

                ixlims=0;
                iylims=0;

                xlims = [0 1.2];
                ylims = [0 8.5];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
                
            case 'CAS hotwire matches'
                i_override_vertical_coord=0;
                
                CAS_LWC = LWC_dist_cas_cutoff';

                %hotwire LWC
                hotwire_LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                hotwire_LWC(hotwire_LWC<-0.05)=NaN;
                hotwire_LWC(hotwire_LWC>2.5)=NaN; %is an LWC point at 4 g/m3 that must be wrong
                
                %only calculate ratio when we have a reasonable amount of LWC
                inds_ratio = find(CAS_LWC>0.02 & hotwire_LWC>0.05);
                ratio = CAS_LWC(inds_ratio)./hotwire_LWC(inds_ratio);
                %ratio(ratio<0)=0;
                
                cas_lwc_xdat = CAS_mode_diameter(inds_ratio);
                
                xdat(1).x = cas_lwc_xdat;
                ydat(1).y = ratio;
                
                %example of binning the data
                [meanvals_ratio,med_vals,maxvals_ratio,max_inds_ratio,mid_points,nvals]=bin_data(cas_lwc_xdat,ratio,CAS_bins);
                
                %trying different thresholds
                inds_ratio2 = find(CAS_LWC>0.02 & hotwire_LWC>0.1);
                ratio2 = CAS_LWC(inds_ratio2)./hotwire_LWC(inds_ratio2);
                cas_lwc_xdat2 = CAS_mode_diameter(inds_ratio2);                
                [meanvals_ratio2,med_vals2,maxvals_ratio2,max_inds_ratio2,mid_points,nvals2]=bin_data(cas_lwc_xdat2,ratio2,CAS_bins);
                
                inds_ratio3 = find(CAS_LWC>0.1 & hotwire_LWC>0.1);
                ratio3 = CAS_LWC(inds_ratio3)./hotwire_LWC(inds_ratio3);
                cas_lwc_xdat3 = CAS_mode_diameter(inds_ratio3);                
                [meanvals_ratio3,med_vals3,maxvals_ratio3,max_inds_ratio3,mid_points,nvals3]=bin_data(cas_lwc_xdat3,ratio3,CAS_bins);
                
                
                
                xlabel_name = 'CAS mode size (/mum)';
                ylabel_name = 'Ratio CAS to hotwire';
                                


                              
                colour_data_var='';
                
                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];

                ixlims=0;
                iylims=1;

                xlims = [0 1.2];
                ylims = [0 8.5];
                
                icolour = 0;
%                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                

                                
                                
                
            case 'Vertical windspeed (m s^{-1})';
                xdat(1).x = interp1(time_flt*3600,w2_turb,CAS_time_all);
                
                xlabel_name = plot_type;
                
                colour_data_var='LWC';
                
                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];

                ixlims=0;
                iylims=1;

                xlims = [0 1.2];
                ylims = [2.1 3.5];
                
                LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                LWC(LWC<-0.05)=NaN;
                LWC(LWC>2.5)=NaN; %is an LWC point at 4 g/m3 that must be wrong


                icolour = 1;
                colour_data = LWC;  %should perhaps do a common bit of code for choosing different
                %colour data
                
            case 'Aircraft altitude'
                eval(['X_flt = X_flt' flight_no ';']);
                eval(['Y_flt = Y_flt' flight_no ';']);
                D_flt = sqrt( (X_flt-X_pos).^2 + (Y_flt-Y_pos).^2 );

                colour_data_var='Wind direction (degrees)';
%                colour_data_var='Wind speed (m s^{-1})';
%                colour_data_var='Mean ice size (microns)';
                colour_data_var='Total ice number (L^{-1})';



                di=20;
                di=1;
                inds2=[di:di:length(inds)];
                inds=inds(inds2);
                
                x_axis_type='Distance (km)';
                x_axis_type='Time (UTC)';                
                
                switch x_axis_type
                    case 'Distance (km)'
                        xdat(1).x = D_flt(inds); %distance
                    case 'Time (UTC)'
                        xdat(1).x = time_flt(inds); %distance
                end

                ydat(1).y = dat_flt(inds,col_alt);

                title_name = [plot_type ' for flight ' flight_no ' with ' colour_data_var];
                xlabel_name = x_axis_type;
                ylabel_name = 'Altitude (m)';

                ixlims=0;
                iylims=0;

                xlims = [0 1.2];
                ylims = [0 1.2];

                iclims=0;

                icolour = 1;
                iset_cmap=0;
                lb_map=lbmap(256,'redblue'); %nice colormap for colorblind people
                cmap=flipud(lb_map);

                switch colour_data_var
                    case 'Wind direction (degrees)'
                        colour_data = dat_flt(inds,col_winddir)+180;
                        iclims=1;
                        clims=[0 360];
                    case 'Wind speed (m s^{-1})'
                        colour_data = dat_flt(inds,col_wind);
                        iclims=1;
                        clims=[0 20];
                    case 'Mean ice size (microns)'
                        %                        [i1,i2]=findheight_nearest(CIP_time_Jonny/3600,min(time_flt100),max(time_flt100));
                        %                        inds=[i1:i2];


                        colour_data = interp1(CIP_time_Jonny/3600,mean_ice_size,time_flt);

                    case 'Total ice number (L^{-1})'

                        colour_data = interp1(CIP_time_Jonny2/3600,1000*ice_no_Jonny,time_flt(inds));
                        
                        indsJonny2 = find( CIP_time_Jonny/3600>=time_flt((1)) & CIP_time_Jonny/3600<=time_flt((end)) );
                        colour_data = 1000*ice_no_Jonny(indsJonny2);
                        
                        xdat(1).x = interp1(time_flt,xdat(1).x,CIP_time_Jonny(indsJonny2)/3600);
                        ydat(1).y = interp1(time_flt,ydat(1).y,CIP_time_Jonny(indsJonny2)/3600);                        
                        

                end




            case 'CAS LWC vs hotwire LWC'



                %CAS LWC
                %                xdat(1).x = LWC_dist_cas;
                xdat(1).x = (1/1)*LWC_dist_cas_cutoff;

                %hotwire LWC
                ydat(1).y = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                ydat(1).y(ydat(1).y<-0.05)=NaN;

                

                title_name = ['Flight ' flight_no];
                xlabel_name = ['LWC CAS size dist ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum (g m^{-3})'];
                ylabel_name = 'LWC hotwire (g m^{-3})';

                ixlims=1;
                iylims=1;

                xlims = [0 1.2];
                ylims = [0 1.2];
                                
%                xlims = [0 2.5];
%                ylims = [0 2.5];
                
                 icolour = 1;
                colour_data = CAS_mean_diameter;
%                colour_data = CAS_total_number{icas_count};


                
                idat=1;
                ismooth_x(idat)=0;
                ismooth_y(idat)=0;                
                
                if ismooth_x(idat)==1 | ismooth_y(idat)==1
                    nfilter=40; bfilter=ones([1 nfilter])*1/nfilter;
                    if ismooth_x==1
                        xdat(idat).x=filter(bfilter,1,xdat(idat).x);
                        title_name = [title_name ' SMOOTHED X'];
                    end
                    if ismooth_y==1
                        ydat(idat).y=filter(bfilter,1,ydat(idat).y);
                        title_name = [title_name ' SMOOTHED Y'];                        
                    end

                    ydat(idat).y(end-nfilter+1:end)=[]; xdat(idat).x(end-nfilter+1:end)=[];
                    colour_data(end-nfilter+1:end)=[];
                    ydat(idat).y(1:nfilter)=[]; xdat(idat).x(1:nfilter)=[];
                    colour_data(1:nfilter)=[];
                end
                
                

                ione_to_one_line=1; %draw a one-to-one line

               

                iadd_line=1;

                nanx = isnan(xdat(1).x);
                nany = isnan(ydat(1).y');
                not_nan = find(nanx==0 & nany==0);

                P = polyfit(xdat(1).x(not_nan),ydat(1).y(not_nan)',1);

                xline = [0:xlims(2)/50:xlims(2)];
                yline = P(1)*xline + P(2);

                if iadd_line==1
                    title_name = [title_name ' ' num2str(P(1)) 'x+' num2str(P(2))];
                end
                
                
                
                if iapplied_ratios==1
                    title_name = [title_name ' CORR. RATIOS APPLIED'];
                end
                
                if ishifted_bins==1
                    title_name = [title_name ' BINS SHIFTED ONE LOWER'];
                end
                
        case 'CAS LWC vs vertical windspeed'



                %CAS LWC

                xdat(1).x = interp1(time_flt*3600,w2_turb,CAS_time_all);
                
                ydat(1).y = (1/1)*LWC_dist_cas_cutoff;

                

                title_name = ['Flight ' flight_no];
                xlabel_name = ['Vertical windspeed (m s^{-1})'];
                ylabel_name = 'LWC CAS (g m^{-3})';

                ixlims=0;
                iylims=1;

                xlims = [0 1.2];
                ylims = [0 0.5];
                                
%                xlims = [0 2.5];
%                ylims = [0 2.5];
                
                 icolour = 0;
                colour_data = CAS_mean_diameter;
%                colour_data = CAS_total_number{icas_count};
    
 case 'Hotwire LWC vs vertical windspeed'



                %CAS LWC
                %                xdat(1).x = LWC_dist_cas;
                xdat(1).x = interp1(time_flt*3600,w2_turb,CAS_time_all);

                %hotwire LWC
                ydat(1).y = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                ydat(1).y(ydat(1).y<-0.05)=NaN;

                

                title_name = ['Flight ' flight_no];
                xlabel_name = ['Vertical windspeed (m s^{-1})'];
                ylabel_name = 'LWC hotwire (g m^{-3})';

                ixlims=0;
                iylims=1;

                xlims = [0 1.2];
                ylims = [0 0.5];
                                
%                xlims = [0 2.5];
%                ylims = [0 2.5];
                
                 icolour = 0;
                colour_data = CAS_mean_diameter;
%                colour_data = CAS_total_number{icas_count};

            case 'PADS LWC vs hotwire LWC'



                %CAS LWC
                %                xdat(1).x = LWC_dist_cas;
                xdat(1).x = (1/1)*LWC_dist_PACS_cutoff;

                %hotwire LWC
                ydat(1).y = interp1(PACS_TAS_time,LWC',PACS_time);
                ydat(1).y(ydat(1).y<-0.05)=NaN;

                title_name = ['PADS flight 433'];
                xlabel_name = ['LWC CAS size dist ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum (g m^{-3})'];
                ylabel_name = 'LWC hotwire (g m^{-3})';

                ixlims=1;
                iylims=1;

                xlims = [0 2.2];
                ylims = [0 2.2];

                ione_to_one_line=1; %draw a one-to-one line

                icolour = 0;
                colour_data = CAS_mean_diameter;

                iadd_line=1;

                nanx = isnan(xdat(1).x);
                nany = isnan(ydat(1).y');
                not_nan = find(nanx==0 & nany==0);

                P = polyfit(xdat(1).x(not_nan),ydat(1).y(not_nan)',1);

                xline = [0:xlims(2)/50:xlims(2)];
                yline = P(1)*xline + P(2);

                if iadd_line==1
                    title_name = [title_name ' ' num2str(P(1)) 'x+' num2str(P(2))];
                end

            case 'Wind direction (degrees) vs height'     
                
                xdat(1).x = dat_flt(:,col_winddir)+180;
                ydat(1).y = dat_flt(:,col_alt);

                title_name = ['Wind direciton for flight ' flight_no ' vs height'];                    
                ylabel_name = 'Altitude (km)';
                xlabel_name = plot_type;
                
                ixlims=1;
                iylims=0;

                xlims = [0 360];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line




                icolour = 0;
                colour_data = 1000*ice_mass_Jonny(indsJonny2);   
                
            case 'Equivalent potential temperature (K) (Humicap)'
                
                T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                qv=interp1(dat_flt(:,1)/1e3,qv_flt_humi,CAS_time_all(indsCAS))';
%                qv = SatVapPress(T,'goff','liq',P,1)/f;
                
                xdat(1).x = equivalent_potemp(T,P,qv);

                title_name = [plot_type flight_no ' vs height'];                    
                xlabel_name = plot_type;
                
                switch vertical_coord
                    case 'height'    
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                        ylabel_name = 'Altitude (km)'; 
                        iylims=1;                        
                        ylims = [2.5 3.5];
                    case 'temperature'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        ylabel_name = 'Temperature (^oC)';   
                        iy_reverse=1;                        
                    case 'pressure'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ylabel_name = 'Pressure (hPa)';                        
                        iy_reverse=1;    
                    case 'potemp'
                        ylabel_name = 'Potential Temperature (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ydat(1).y = T.*(1000./P).^0.286;                         
                        iy_reverse=0;     
                end
                
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 360];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line

                LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                LWC(LWC<-0.05)=NaN;
                LWC(LWC>2.5)=NaN; %is an LWC point at 4 g/m3 that must be wrong


                icolour = 1;
                colour_data = LWC;
                
           case 'Equivalent potential temperature (K) (Saturation when LWC)'
                
                T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                
                LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                i_ignore = find(LWC<-0.05 & LWC>2.5); 
                LWC(i_ignore)=NaN;
                iLWC=find(LWC>=0.08);
                
                qv=interp1(dat_flt(:,1)/1e3,qv_flt_humi,CAS_time_all(indsCAS))';
                %only use saturation qv if LWC greater than a threshold
                qv(iLWC) = SatVapPress(T(iLWC),'goff','liq',P(iLWC),1)/f;
                
                xdat(1).x = equivalent_potemp(T,P,qv);

                title_name = [plot_type flight_no ' vs height'];                    
                xlabel_name = plot_type;
                
                switch vertical_coord
                    case 'height'    
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                        ylabel_name = 'Altitude (km)'; 
                        iylims=1;                        
                        ylims = [2.5 3.5];
                    case 'temperature'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        ylabel_name = 'Temperature (^oC)';   
                        iy_reverse=1;                        
                    case 'pressure'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ylabel_name = 'Pressure (hPa)';                        
                        iy_reverse=1;    
                    case 'potemp'
                        ylabel_name = 'Potential Temperature (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ydat(1).y = T.*(1000./P).^0.286;                         
                        iy_reverse=0;   
                    case 'vertical velocity'
                        ylabel_name = 'Updraught (m s^{-1})';
                        ydat(1).y = interp1(time_flt*3600,w2_turb,CAS_time_all);
                        iy_reverse=0;       
                end
                
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 360];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line

                LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                LWC(LWC<-0.05)=NaN;
                LWC(LWC>2.5)=NaN; %is an LWC point at 4 g/m3 that must be wrong
%                LWC(LWC>0.1)=NaN; %is an LWC point at 4 g/m3 that must be wrong
                
                w_obs = interp1(time_flt*3600,w2_turb,CAS_time_all);
                w_obs(LWC<0.1)=NaN;
                
                alt=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                alt(LWC<0.1)=NaN;
                alt(LWC>2.5)=NaN;
                alt(alt<2.5)=NaN;
                
                

                icolour = 1;
                colour_data = LWC; 
%                colour_data = w_obs;   
%                colour_data = alt;   
                
          case 'Estimated displacement (m) vs height'
                  %is an LWC point at 4 g/m3 that must be wrong                                
                LWC = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                i_ignore = find(LWC<-0.05 & LWC>2.5); 
                LWC(i_ignore)=NaN;
                iLWC=find(LWC>=0.08);
                
                T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
% for qv will use satuartion qv when have LWC and vapour measurement when do no
% since the qv measurement is likely to be unreliable in cloud
                clear qv
                qv=interp1(dat_flt(:,1)/1e3,qv_flt_humi,CAS_time_all(indsCAS))';
                %replace value for points with LWC
                qv(iLWC) = SatVapPress(T(iLWC),'goff','liq',P(iLWC),1)/f;                
                equiv = equivalent_potemp(T,P,qv);
                heights=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))';
                
              

                
                w_obs = interp1(time_flt*3600,w2_turb,CAS_time_all);
                w_obs(i_ignore)=NaN;
                
                time_obs=CAS_time_all;
                time_obs(i_ignore)=NaN;
                
            %idealised soundings from wave_temperature... Matlab model
                eq_use = eq_h_fp;
                eq_use = eq_h_humi; %the above two are very similar
                h_use = h_wave_model;
            %actual sounding
%                eq_use = equiv_sound_humi;
%                h_use = Zsound;
                
                min_LWC=0;
                i_origin_height=0;
                
                for i=1:length(equiv)
                    if equiv(i)>=min(eq_use) & equiv(i)<=max(eq_use) & LWC(i)>min_LWC
%                        isource = findheight(eq_h_fp,equiv(i));

                        if i_origin_height==0 %plot displacement
                            xdat(1).x(i) = heights(i) - interp1(eq_use,h_use,equiv(i)); 
                            title_name = [plot_type flight_no ' vs height'];  
                            xlabel_name = plot_type;
                        else %plot height of origin
                            xdat(1).x(i) = interp1(eq_use,h_use,equiv(i));     
                            title_name = ['Origin height vs height'];
                            xlabel_name = 'Origin height (m)';
                        end
                    else
                        xdat(1).x(i)=NaN;  %might want to guess using an idealised profile
                        %but mark on them with crosses
                    end
                end

                                 


                
                switch vertical_coord
                    case 'height'    
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                        ylabel_name = 'Altitude (km)'; 
                        iylims=1;
                        ylims = [2 4];
                    case 'temperature'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        ylabel_name = 'Temperature (^oC)';   
                        iy_reverse=1;                        
                    case 'pressure'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ylabel_name = 'Pressure (hPa)';                        
                        iy_reverse=1;    
                    case 'potemp'
                        ylabel_name = 'Potential Temperature (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ydat(1).y = T.*(1000./P).^0.286;                         
                        iy_reverse=0;    
                    case 'lwc'    
                        ydat(1).y = LWC;
                        ylabel_name = 'LWC (g m^{-3})'; 
                        iylims=0;                        
                        ylims = [2.5 3.5];    
                        
                end
                
                
                
                ixlims=0;
                xlims = [0 360];                
                
                

                ione_to_one_line=0; %draw a one-to-one line

%                LWC(LWC<0.1)=NaN; %is an LWC point at 4 g/m3 that must be wrong

                icolour = 1;
                colour_data = LWC;    
%                colour_data = alt;     
%                colour_data = w_obs;                                 
%                colour_data = time_obs/3600; 

%only plot the data within a certain range
    %colour_data(abs(colour_data)>1)=NaN;
%    colour_data(abs(colour_data)<0.1)=NaN;
                iclims=1;
                clims=[2 3.5];
                clims=[0 0.3];                
                
                
                
                
                
                
            case 'Vapour mixing ratio (kg kg^{-1}) (Hygrometer)'
                xdat(1).x = interp1(dat_flt(:,1)/1e3,qv_flt_fp,CAS_time_all(indsCAS))';

                title_name = ['Vapour mixing ratio (hygrometer) for flight ' flight_no ' vs height'];                    
                xlabel_name = plot_type;
                
                switch vertical_coord
                    case 'height'    
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                        ylabel_name = 'Altitude (km)'; 
                        iylims=1;                        
                        ylims = [2.5 3.5];
                    case 'temperature'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        ylabel_name = 'Temperature (^oC)';   
                        iy_reverse=1;                        
                    case 'pressure'
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ylabel_name = 'Pressure (hPa)';                        
                        iy_reverse=1;    
                    case 'potemp'
                        ylabel_name = 'Potential Temperature (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ydat(1).y = T.*(1000./P).^0.286;                         
                        iy_reverse=0;     
                end
                
                
                
                ixlims=0;
                iylims=0;

                xlims = [0 360];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line




                icolour = 1;
                colour_data = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))';

            case 'Ice mean diameter (\mum)'     
                
                time_base_scatter = CIP_time_Jonny2(indsJonny2); %nneds to be in secs
                xdat(1).x = mean_ice_size(indsJonny2);
                ydat(1).y = interp1(time_flt,dat_flt(:,col_alt),CIP_time_Jonny2(indsJonny2)/3600)'/1e3;
                ydat(1).y = interp1(time_flt,dat_flt(:,col_temp),CIP_time_Jonny2(indsJonny2)/3600)';
                ylabel_name = 'Altitude (km)';
                ylabel_name = 'Temperature (^{o}C)';                
                xlabel_name = plot_type;
                

                
                iy_reverse=1;
                
                ixlims=1;
                iylims=0;

                xlims = [0 1500];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line

            

                
                icolour = 1;
                colour_data = 1000*ice_no_Jonny(indsJonny2);
                title_name = [datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,15) ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,15) ...
                    ', Ice mean diam& no.(colours; \mum), flight ' flight_no];
                iclims=1;
                clims=[0 10];
                
                inds_limit = find(colour_data>0.5);
                ilimit_data=1; %limit the date to inds_limit
                
%                colour_data = interp1(time_flt,dat_flt(:,col_alt),CIP_time_Jonny2(indsJonny2)/3600)/1000';
%                title_name = ['Ice number & altitude (colours; km) for flight ' flight_no];
%                iclims=1;
%                clims=[2 3.3];


             case 'Ice number binned (L^{-1})'     %note this uses the same time bins as calculated for the 
                 %plot_flight_segments_for_different_flights &
                 %timetimeseries_dan plot. So run that first to get the
                 %correct envelope peiod of time and the time bins used -
                 %then can narrow down the range for the scatter plot,  but
                 %keep the same time bins - this means that the plotted
                 %values are the same in the scatter plot as in the
                 %timeseries (for the same binning period)
                 
                         ice_tot = 'ice';
%                        ice_tot = 'tot';

                        ismooth(idat)=0;  
                        
                        %want to use the same bins as in the timeseries for comparison                        
                        %%%%%%%%%%
                        NT=60; %no. secs
                        %%%%%%%%%%
                        %run timeseries_dan first to get inds2
                        
                        ijoin_points=0; %put a line through the points
                        
                                              
%                        [i1,i2]=findheight_nearest(CIP_time_Jonny2/3600,min(time_flt),max(time_flt));
%                        inds2=[i1:i2];

                        %inds2 are taken from the plot_flight_segments_for_different_flights plot
                            
                        
%stuff for the size range of ice required      
                        size_range_Jonny = [87.5 1600];
%                        size_range_Jonny = [0 1600];                        
                        i1 = find(CIP_size_bins_Jonny>=size_range_Jonny(1));
                        i2 = find(CIP_size_bins_Jonny<=size_range_Jonny(2))-1;                        
                        
                        inds_size_Jonny=[i1(1):i2(end)]; %for the PSD array
                        actual_sizes=num2str(CIP_size_bins_Jonny([i1(1) i2(end)]));
                        size_label=[actual_sizes(1,:) ' to ' actual_sizes(2,:) ' \mum'];
                        
                        
%get the ice data and the number of crystals for the (shorter) scatter_plot range
                        %below uses indsJonnny, which is the shorter range chosen
                        %in scatter_plot
                        switch ice_tot
                            case 'ice'
                                ylab=['Binned ice number (L^{-1})'];
                                ice_dat = 1000*sum(ice_PSD(inds_size_Jonny,indsJonny),1);
                            case 'tot'
                                ylab=['Total number (L^{-1})'];
                                ice_dat = 1000*sum(tot_PSD(inds_size_Jonny,indsJonny),1);
                        end                        
                        ice_accepts = sum(ice_no_particles_PSD(inds_size_Jonny,indsJonny),1); %no. crystals
                        tas_dat = TAS_Jonny(indsJonny);

                        CIP_time_plot = CIP_time_Jonny(indsJonny)/3600;
                                                
                   %is this used 
                        [i1,i2]=findheight_nearest(CIP_time_Jonny/3600,min(time_flt),max(time_flt));
                        inds2_timebins=[i1:i2];                                                
                    %??
                    
                    %now calculate the time bins that would have been used
                    %in plot_flight_segments_for_different_flights
                    time_bins_old = [CIP_time_Jonny(inds2(1))/3600:NT/3600:CIP_time_Jonny(inds2(end))/3600];
                        
                    %find the indices in the old time bins that contain the
                    %new shorter scatter plot range (indicated by
                    %times(1:2)
                    ibins_old = find( time_bins_old >= times(1) & time_bins_old <= times(2) );                                                                
                    time_bins = time_bins_old(ibins_old); %
                    
%or use this below to just calculate its own time bins
                    %time_bins = [CIP_time_Jonny2(indsJonny2(1))/3600:NT/3600:CIP_time_Jonny2(indsJonny2(end))/3600];


    %now do some binning
    %first the mean airspeed in the time bin
                [meanvals_tas,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny(indsJonny)/3600,tas_dat,time_bins);
    %now the number of crystals
                [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny(indsJonny)/3600,ice_accepts,time_bins);
                num_dat(1).dat = sum_vals; %number of each statistic

%                         %absolute errors (estimated when N=0)
%                         i0=find(ice_accepts==0);
                 Aerr_all = ice_dat.*sqrt(1./ice_accepts +0.04); %combined absolute error due to Poisson and 20% in TAS
                 Aerr_all = ice_dat.*sqrt(1./ice_accepts); %combined absolute error due to Poisson
                 %                         SV_mean = tas_dat*160/1000; %mean sample volume in litres assuming 160 mm2 for area
%                         Aerr_all(i0)= 1.8./SV_mean(i0); %set absolute error to 1.8 crystals when have no crystals
%                         %then work out an error in the concentration from the mean sample volume     
%

%now the binned ice concentrations (as already calculated in plot_flight...
                [meanvals_Nice,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(CIP_time_Jonny(indsJonny)/3600,ice_dat,time_bins);

                
%%%                
%estimate the error based on Poisson stats (sqrt of the number of crystals)                
                        Aerr = meanvals_Nice.*sqrt(1./num_dat(1).dat); %combined absolute error due to Poisson and 20% in TAS
%                        Aerr = meanvals_Nice.*sqrt(1./num_dat(1).dat +0.04); %combined absolute error due to Poisson and 20% in TAS                        

%estimate an error bar for zero counts
                        i0=find(num_dat(1).dat==0);
                        SV_mean = meanvals_tas*160/1000; %mean sample volume in litres assuming 160 mm2 for area
                        Aerr(i0)= 1.8./SV_mean(i0)/NT; %set absolute error to 1.8 crystals (from Poisson) when have no crystals
                        %then work out an error in the concentration from the mean sample volume    
                                                
%                [meanvals_Nice,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals,std_dev,mean_error]=bin_data2(CIP_time_Jonny(indsJonny)/3600,ice_dat,time_bins,Aerr_all);
                        error_bar(1).dat=Aerr;
%                        error_bar(1).dat=std_dev;
                        iplot_error_bar=1;
                       
 %%%                       
                
                        
                        
%                       ioverride_DY=1;
%                       DY=0; %
                        iadd_nums_above=1;

%bin the temperature data - use the averages of the temperature in each NT second time bin                      
                [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(time_flt,dat_flt(:,col_temp),time_bins);
                temperature_data = meanvals; %mean temperature of each bin
                                                
                        
                ixlims=0;
                iylims=0;

                xlims = [0 18]; %flight 100
                xlims = [0 0.8]; %flight 99                
%                xlims = [0 20];                
                ylims = [0 0.9];    
                        
                        
                
                
%additional option to bin all of the NT second time averages into
%temperature bins
                itemp_means=0;                        
                if itemp_means==1
                    %time_base_scatter = CIP_time_Jonny2(indsJonny2); %nneds to be in secs
%                    time_base_scatter = mid_points; %time bins                   
                    
                    temp_bins=[-22.25:0.5:-16];                                    
                    
%calculate the combined error of all the different NT second samples
                    [meanvals_tempbins,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals,mean_ig_zeros,std_dev,mean_error]=bin_data2(temperature_data,meanvals_Nice,temp_bins,Aerr);                                        

%now calculate the total number of crystals in each temperature bin                    
                    [meanvals,med_vals,maxvals,max_inds,mid_points,nvals,sum_vals]=bin_data(temperature_data,num_dat(1).dat,temp_bins);
                    temp_bins_mids=mid_points;
                    num_dat(1).dat=sum_vals;
                    error_test = meanvals_tempbins.*sqrt(1./num_dat(1).dat);

%                    error_bar(1).dat=mean_error;
%                    error_bar(1).dat=std_dev;
                    error_bar(1).dat=error_test;
                    iplot_error_bar=1;
                    
                    ydat(1).y = temp_bins_mids;
                    xdat(1).x = meanvals_tempbins;
%                    xdat(1).x = mean_ig_zeros; %mean if ignore the zeros
                    
                    iylims=1;              
                    ylims = [temp_bins(1) temp_bins(end)];    
                    
                    ijoin_points=1;
                
                
                else
 %                   time_base_scatter = mid_points; %time bins
                    ydat(1).y = temperature_data;
                    xdat(1).x = meanvals_Nice;

                end
%                ydat(1).y = interp1(time_flt,dat_flt(:,col_temp),CIP_time_Jonny2(indsJonny2)/3600)';
%                ylabel_name = 'Altitude (km)';
                ylabel_name = 'Temperature (^{o}C)';                
                xlabel_name = [plot_type ', ' num2str(NT) ' sec bins'];
                
                ylab = [num2str(NT) ' sec ' ylab];
                

                
                iy_reverse=1;
                


                ione_to_one_line=0; %draw a one-to-one line

                scatter_properties_string = ['100,''g'',''filled'''];  

                
                icolour = 0;
                colour_data = 1000*ice_mass_Jonny(indsJonny2);
                title_name = ['Ice mass (colours; mg m^{-3}) and number for flight ' flight_no ' for ' datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,13)...
                    ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,13)];


%                colour_data = mean_ice_size(indsJonny2);
                title_name = [datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,15) ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,15) ...
                    ', Ice no.&mean diam (colours; \mum), flight ' flight_no];
                iclims=0;
                clims=[0 500];
                clims=[0 900];
                
%                colour_data = interp1(time_flt,dat_flt(:,col_alt),CIP_time_Jonny2(indsJonny2)/3600)/1000';
%                title_name = ['Ice number & altitude (colours; km) for flight ' flight_no];
%                iclims=1;
%                clims=[2 3.3];

                 title_name = ['Ice no.(' size_label ') & no. samples for flight ' flight_no ' for ' datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,13)...
                    ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,13)];   

                 i_override_vertical_coord=0;



                
            case 'Ice number (L^{-1})'     
                
                time_base_scatter = CIP_time_Jonny2(indsJonny2); %nneds to be in secs
                xdat(1).x = 1000*ice_no_Jonny(indsJonny2);
                ydat(1).y = interp1(time_flt,dat_flt(:,col_alt),CIP_time_Jonny2(indsJonny2)/3600)'/1e3;
                ydat(1).y = interp1(time_flt,dat_flt(:,col_temp),CIP_time_Jonny2(indsJonny2)/3600)';
                ylabel_name = 'Altitude (km)';
                ylabel_name = 'Temperature (^{o}C)';                
                xlabel_name = plot_type;
                

                
                iy_reverse=1;
                
                ixlims=1;
                iylims=0;

                xlims = [0 18]; %flight 100
                xlims = [0 0.8]; %flight 99                
%                xlims = [0 20];                
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line

            

                
                icolour = 1;
                colour_data = 1000*ice_mass_Jonny(indsJonny2);
                title_name = ['Ice mass (colours; mg m^{-3}) and number for flight ' flight_no ' for ' datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,13)...
                    ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,13)];


                colour_data = mean_ice_size(indsJonny2);
                title_name = [datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,15) ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,15) ...
                    ', Ice no.&mean diam (colours; \mum), flight ' flight_no];
                iclims=1;
                clims=[0 500];
                clims=[0 900];
                
%                colour_data = interp1(time_flt,dat_flt(:,col_alt),CIP_time_Jonny2(indsJonny2)/3600)/1000';
%                title_name = ['Ice number & altitude (colours; km) for flight ' flight_no];
%                iclims=1;
%                clims=[2 3.3];



                
                
            case 'Ice mass (mg m^{-3})'     
                
                xdat(1).x = 1000*ice_mass_Jonny(indsJonny2);
                ydat(1).y = interp1(time_flt,dat_flt(:,col_alt),CIP_time_Jonny2(indsJonny2)/3600)'/1e3;

                title_name = ['Ice number (colours; L^{-1}) and mass for flight ' flight_no ' for ' datestr(CIP_time_Jonny2(indsJonny2(1))/3600/24,13)...
                    ' to ' datestr(CIP_time_Jonny2(indsJonny2(end))/3600/24,13)];
                ylabel_name = 'Altitude (km)';
                xlabel_name = plot_type;
                
                ixlims=0;
                iylims=0;

                xlims = [0 0.9];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line




                icolour = 1;
                colour_data = 1000*ice_no_Jonny(indsJonny2);   
                
                inds_limit = find(colour_data>0.1);
                ilimit_data=1; %limit the date to inds_limit



            case 'CAS LWC vs conc.'


                %CAS LWC
                xdat(1).x = LWC_dist_cas;

                %total concentration above 1 micron diameter
                ydat(1).y = CAS_total_number{icas_count};


                title_name = ['Flight ' flight_no];
                xlabel_name = 'LWC CAS size dist (g m^{-3})';
                ylabel_name = 'Total number (cm^{-3})';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

                icolour = 1;
                colour_data = CAS_mean_diameter;


            case 'CAS LWC vs conc. cutoff'


                %CAS LWC
                xdat(1).x = LWC_dist_cas_cutoff

                %total concentration above 1 micron diameter
                ydat(1).y = CAS_total_number{icas_count};


                title_name = ['Flight ' flight_no];
                xlabel_name = ['LWC CAS ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum cut off (g m^{-3})'];
                ylabel_name = 'Total number (cm^{-3})';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

                icolour = 1;
                colour_data = CAS_mean_diameter;



            case 'Hotwire LWC vs conc.'


                %CAS LWC
                xdat(1).x = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                xdat(1).x(xdat(1).x<-0.05)=NaN;

                %total concentration above 1 micron diameter
                ydat(1).y = CAS_total_number{icas_count};


                title_name = ['Flight ' flight_no];
                xlabel_name = 'Hotwire LWC (g m^{-3})';
                ylabel_name = 'Total number (cm^{-3})';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

                icolour = 0;
                colour_data = CAS_mean_diameter;

            case 'Hotwire LWC vs mode'

                %CAS LWC
                xdat(1).x = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                xdat(1).x(xdat(1).x<-0.05)=NaN;

                %total concentration above 1 micron diameter
                ydat(1).y = CAS_mode_diameter;


                title_name = ['Flight ' flight_no];
                xlabel_name = 'Hotwire LWC (g m^{-3})';
                ylabel_name = 'Mode diameter (\mum)';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

            case 'Hotwire LWC vs mean'

                %CAS LWC
                xdat(1).x = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all);
                xdat(1).x(xdat(1).x<-0.05)=NaN;

                %total concentration above 1 micron diameter
                ydat(1).y = CAS_mean_diameter;


                title_name = ['Flight ' flight_no];
                xlabel_name = 'Hotwire LWC (g m^{-3})';
                ylabel_name = 'Mean diameter (\mum)';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

                icolour = 1;
                colour_data = CAS_mean_diameter;

            case 'CAS LWC vs mode'

                %CAS LWC
                xdat(1).x = LWC_dist_cas;
                ydat(1).y = CAS_mode_diameter;


                title_name = ['Flight ' flight_no];
                xlabel_name = 'CAS LWC (g m^{-3})';
                ylabel_name = 'Mode diameter (\mum)';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

            case 'CAS LWC vs mean'

                %CAS LWC
                xdat(1).x = LWC_dist_cas;
                ydat(1).y = CAS_mean_diameter;


                title_name = ['Flight ' flight_no];
                xlabel_name = 'CAS LWC (g m^{-3})';
                ylabel_name = 'Mean diameter (\mum)';

                ixlims=1;
                iylims=0;

                xlims=[0 0.9];

                icolour = 1;
                colour_data = CAS_mean_diameter;

            case 'Mode vs number'


                %CAS LWC
                xdat(1).x = CAS_mode_diameter;
                ydat(1).y = CAS_total_number{icas_count};


                title_name = ['Flight ' flight_no];
                xlabel_name = 'Mode diameter (\mum)';
                ylabel_name = 'Total number (cm^{-3})';

                ixlims=0;
                iylims=0;


            case 'CAS LWC vs height'

                icut_off_lwc=1;

                %CAS LWC
                switch icut_off_lwc
                    case 1
                        xdat(1).x = (1/1)*LWC_dist_cas_cutoff(indsCAS);
                        xlabel_name = ['LWC CAS size dist ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum (g m^{-3})'];
                    otherwise
                        xdat(1).x = LWC_dist_cas{1}(indsCAS);
                        xlabel_name = ['LWC CAS size dist (g m^{-3})'];
                end

                %height
                ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;

                title_name = ['Flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                ylabel_name = 'Altitude (km)';

                ixlims=1;
                iylims=1;

                xlims = [0 0.9];
                ylims = [3 3.5];

                ione_to_one_line=0; %draw a one-to-one line




                icolour = 1;
                colour_data = CAS_mean_diameter(indsCAS);
                colour_data = CAS_total_number{icas_count}(indsCAS);
                
             case 'CAS LWC'

                icut_off_lwc=1;

                %CAS LWC
                switch icut_off_lwc
                    case 1
                        xdat(1).x = (1/1)*LWC_dist_cas_cutoff(indsCAS);
                        xlabel_name = ['LWC CAS size dist ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum (g m^{-3})'];
                    otherwise
                        xdat(1).x = LWC_dist_cas{1}(indsCAS);
                        xlabel_name = ['LWC CAS size dist (g m^{-3})'];
                end

               
                title_name = ['Flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                
                ixlims=1;
                iylims=0;
                xlims = [0 0.9];
                
                
                switch vertical_coord
                    case 'temperature'
                        ylabel_name = 'Temperature (^oC)';
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        iy_reverse=1;
                    case 'height'
                        ylabel_name = 'Height (km)';
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))/1000';
                        iylims=1;
                        ylims = [2.5 3.5];
                    case 'pressure'
                        ylabel_name = 'Pressure (hPa)';
                        ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        iy_reverse=1;
                    case 'potemp'
                        ylabel_name = 'Potential Temperature (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        ydat(1).y = T.*(1000./P).^0.286;                         
                        iy_reverse=0;    
                     case 'equiv potemp humi'
                        ylabel_name = 'Equivalent Potential Temperature (humi) (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        qv=interp1(dat_flt(:,1)/1e3,qv_flt_humi,CAS_time_all(indsCAS))';
                                                                                     
                        ydat(1).y = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                        
                        
                        iy_reverse=0; 
                        
                        iylims=1;
                        ylims = [295 305];
                        
                     case 'equiv potemp fp'
                        ylabel_name = 'Equivalent Potential Temperature (fp) (K)';
                        T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                        P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                        qv=interp1(dat_flt(:,1)/1e3,qv_flt_fp,CAS_time_all(indsCAS))';
                                                                                     
                        ydat(1).y = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
                        
                        
                        iy_reverse=0;    
                        
                        iylims=1;
                        ylims = [295 305];
                    
                end




                ione_to_one_line=0; %draw a one-to-one line

                icolour = 1;
                colour_data = CAS_mean_diameter(indsCAS);
                colour_data = CAS_total_number{icas_count}(indsCAS);
                

                
                
   
                

            case 'Hotwire LWC'

                icut_off_lwc=1;

                xdat(1).x = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all(indsCAS));
                xdat(1).x(xdat(1).x<-0.05)=NaN;
                xlabel_name = 'Hotwire LWC (g m^{-3})';
                
                ixlims=1;
                xlims = [0 0.55];
                iylims=0;                
                
                              
                title_base_str = ['LWC vs ' vertical_coord];
                
%                title_name = ['Hotwire LWC&CAS no. conc. for D>' num2str(cut_off_size) '\mum (colours;cm^{-3}) for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
%                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                
%                title_name = ['LWC&temperature (colours;^{o}C) for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
%                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];

                title_name = ['LWC&altitude (colours;km) for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                
  title_name = ['LWC&altitude (colours;km) for flight ' flight_no];

                
        switch vertical_coord
            case {'equiv potemp sat','equiv potemp fp','equiv potemp humi'}
                title_name = [title_base_str ' &altitude (colours;km) for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];

                colour_data = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                iclims=1;
                clims=[2.5 3.5];
                
            case 'height'
                title_name = [title_base_str ' &equiv potemp sat (colours;K) for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                
                T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';
                %assume that are at dew point when have LWC (as hygrometer measurements
                %likely to be unreliable                        
                    qv = SatVapPress(T,'goff','liq',P,1)/f;   
                
                colour_data = equivalent_potemp(T,P,qv);
                iclims=1;
                clims=[295 305];
        end
                             


                ione_to_one_line=0; %draw a one-to-one line
                                                                       

                icolour = 1;
%                colour_data = CAS_mean_diameter(indsCAS);
%                colour_data = CAS_total_number{icas_count}(indsCAS);
                %                colour_data = CAS_time_all(indsCAS)/3600;
%                colour_data = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';
                colour_data = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))'/1e3;
                    iclims=1;
                    clims=[2000 3300]/1e3; %flight 102
                    clims=[0 1200]/1e3; %filght 100      
                    clims=[0 4500]/1e3; %filght 101                          
%                colour_data = equivalent_potemp(T,P,qv);
                
    case 'Hotwire LWC vs temperature'

                icut_off_lwc=1;

                xdat(1).x = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all(indsCAS));
                xdat(1).x(xdat(1).x<-0.05)=NaN;
                xlabel_name = 'Hotwire LWC (g m^{-3})';

                %height
                ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all(indsCAS))';

                
                title_name = ['Hotwire LWC vs temperature & droplet number (colours; m) for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                
                ylabel_name = 'Temperature (^{o}C)';

                ixlims=1;
                iylims=0;

                xlims = [0 0.9];
                ylims = [0 0.9];
                
                iy_reverse=1;

                ione_to_one_line=0; %draw a one-to-one line



                icolour = 1;
                colour_data = CAS_mean_diameter(indsCAS);
                colour_data = CAS_total_number{icas_count}(indsCAS);
                %                colour_data = CAS_time_all(indsCAS)/3600;
                colour_data = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))';

            case 'Hotwire LWC vs pressure'

                icut_off_lwc=1;

                xdat(1).x = interp1(CIP_time_all,LWC_CAS_all',CAS_time_all(indsCAS));
                xdat(1).x(xdat(1).x<-0.05)=NaN;
                xlabel_name = 'Hotwire LWC (g m^{-3})';

                %height
                ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all(indsCAS))';

                
                title_name = ['Hotwire LWC&tot num (colours; m) vs Temp. for flight ' flight_no ' for ' datestr(CAS_time_all(indsCAS(1))/3600/24,13)...
                    ' to ' datestr(CAS_time_all(indsCAS(end))/3600/24,13)];
                
                ylabel_name = 'Pressure (hPa)';

                ixlims=1;
                iylims=0;

                xlims = [0 0.9];
                ylims = [0 0.9];
                
                iy_reverse=1;

                ione_to_one_line=0; %draw a one-to-one line



                icolour = 1;
                colour_data = CAS_mean_diameter(indsCAS);
                colour_data = CAS_total_number{icas_count}(indsCAS);
                %                colour_data = CAS_time_all(indsCAS)/3600;
                colour_data = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all(indsCAS))';
                
            case 'CAS number vs height'

                xdat(1).x = CAS_total_number{icas_count};
                xlabel_name = ['Total number (cm^{-3})'];

                %height
                ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all)'/1e3;

                title_name = ['Flight ' flight_no];
                ylabel_name = 'Altitude (km)';

                ixlims=0;
                iylims=0;

                xlims = [0 0.9];
                ylims = [0 0.9];

                ione_to_one_line=0; %draw a one-to-one line




                icolour = 1;
                colour_data = CAS_mean_diameter;


        end
        
%     switch vertical_coord
%         case 'height'
%             ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),time_base_scatter)'/1e3;
%             ylabel_name = 'Altitude (km)';
%             iylims=1;
%             ylims = [2.5 3.5];
%         case 'temperature'
%             ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),time_base_scatter)';
%             ylabel_name = 'Temperature (^oC)';
%             iy_reverse=1;
%             iylims=1;
%             ylims = [-16 -7];
%         case 'pressure'
%             ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
%             ylabel_name = 'Pressure (hPa)';
%             iy_reverse=1;
%         case 'potemp'
%             ylabel_name = 'Potential Temperature (K)';
%             T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),time_base_scatter)';
%             P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
%             ydat(1).y = T.*(1000./P).^0.286;
%             iy_reverse=0;
% 
%         case 'equiv potemp sat'
%             ylabel_name = 'Equivalent Potential Temperature at Saturation (K)';
%             T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),time_base_scatter)';
%             P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
%             %                        qv=interp1(dat_flt(:,1)/1e3,qv_flt_humi,time_base_scatter)';
%             %or assume that are at dew point when have LWC (as hygrometer measurements
%             %likely to be unreliable
% 
%             qv = SatVapPress(T,'goff','liq',P,1)/f;
% 
%             %                        ydat(1).y = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
%             ydat(1).y = equivalent_potemp(T,P,qv);
% 
%             iylims=1;
%             ylims = [295 305];
%             iy_reverse=0;
% 
%         case 'equiv potemp humi'
%             ylabel_name = 'Equivalent Potential Temperature (K)';
%             T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),time_base_scatter)';
%             P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
%             qv=interp1(dat_flt(:,1)/1e3,qv_flt_humi,time_base_scatter)';
%             %or assume that are at dew point when have LWC (as hygrometer measurements
%             %likely to be unreliable
% 
%             %                        qv = SatVapPress(T,'goff','liq',P,1)/f;
% 
%             ydat(1).y = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
% 
%             iylims=1;
%             ylims = [295 305];
%             iy_reverse=0;
% 
% 
%             %bit of extra code to calculate the cloud base by working backwards from a give cloudy point
%             %with known LWC - then can get an idea of the consistency of the background profile or the LWC
%             %measurement (e.g. CAS). Or an idea of how much entrainment there might have been if are confident
%             %of the background profile
%             %so want to solve for the path of the parcel (P,T,qv) working towards the surface
%             ii=find(xdat(1).x>0.27 & xdat(1).x<0.32); %%where xdat contains LWC - finding a particular LWC point here
%             LWCeq=(xdat(1).x(ii));
%             Teq=T(ii);
%             Peq=P(ii);
%             Qeq=qv(ii);
%             Pseq=[Peq(1):0.1e2:1000e2];
%             [LWCad_eq,Tad_eq]=adLWC_PaulLawson_simple(Peq(1),Teq(1),Pseq);
% 
%             %%% need to convert LWC into kg/kg
%             rho_eq=density(Peq(1),Teq(1));
%             LWCeq_kg=LWCeq(1)/1000 /rho_eq;
%             rho_ad=density(Pseq,Tad_eq+273.15);
%             LWCad_eq_kg = LWCad_eq/1000 ./ rho_ad;
% 
%             Qad_eq = Qeq(1)-LWCad_eq_kg; %vapour MR during the descent
% 
% 
%             [minLWC,ii2]=min(abs(LWCad_eq_kg+LWCeq_kg(1))); %looking for when the combination is zero (i.e. all LWC evaporated)
%             %so from index ii2 onwards we have are following the dry adiabat
%             pot0 = (Tad_eq(ii2)+273.15) * (1000e2/Pseq(ii2))^0.286;
%             Tad_eq(ii2+1:end) = pot0./ (1000e2./Pseq(ii2+1:end)).^0.286 - 273.15;
%             %        rho_cb = density(Pseq(ii2),Tad_eq(ii2));
%             Qad_eq(ii2+1:end) = Qad_eq(ii2);
% 
%         case 'equiv potemp fp'
%             ylabel_name = 'Equivalent Potential Temperature (fp) (K)';
%             T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),time_base_scatter)';
%             P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
%             %                        qv=interp1(dat_flt(:,1)/1e3,qv_flt_fp,time_base_scatter)';
% 
%             qv = SatVapPress(T,'goff','liq',P,1)/f;
% 
%             ydat(1).y = ( (T + 2.453e6*qv/1004).*(1e5./P).^0.286 )';
% 
% 
%             iy_reverse=0;
% 
%             iylims=1;
%             ylims = [295 305];
%     end

                
          

        switch plot_type
            case {'Hotwire LWC','CAS LWC'}
                if strcmp(vertical_coord,'height')==1
                iadd_line=1; %flag to add the moist adiabat line 
                %defined by xline and yline
                                
                
                method='choose cb get temp from flight data';
                %                 method='assume temperature gradient with height';
                method = 'choose H get cb and temp from flight data';
                
%                method=''; %set to switch off the adiabat line
                
                switch method
                    case 'choose H get cb and temp from flight data'
                        H_cb = 350;
                        H_cb = 100; %flight 100 tests
                        %H_cb = 250;                        
                        H_cb = 550; %flight 113
                        H_cb = 530; %flight 113    
                        H_cb = 100; %flight 113 
                        H_cb = 3140; %flight 102
                        H_cb = 3055;
                        H_cb = 3155;                        
                        H_cb = 3000;
                        H_cb = 3250;     
                        H_cb = 2850;  
                        H_cb = 2900;  
                        H_cb = 2750;                          
                        H_cb = 2600;                          

use_all_data_for_mean_tb=2;

switch use_all_data_for_mean_tb
    case 0
        Z=dat_flt(inds,col_alt);

        tol=10;
        iz = find(Z<(H_cb+tol) & Z>(H_cb-tol));
        tb = mean(dat_flt(inds(iz),col_temp))+273.15;
        cb = mean(dat_flt(inds(iz),col_press));
        %                        tb=-5+273.15; %overriding the above

    case 1
        Z=dat_flt(:,col_alt);

        tol=10;
        iz = find(Z<(H_cb+tol) & Z>(H_cb-tol));
        tb = mean(dat_flt(iz,col_temp))+273.15;
        cb = mean(dat_flt(iz,col_press));
        %                        tb=-5+273.15; %overriding the above
        
    case 2
        %here are starting at a value determined from a background sounding
        %that represents out of cloud air. If use the data for the whole flight then
        %are mixing cloudy and non-cloudy air - thus the cloud base determined from this
        %is likely to be wrong. Need to either average over only cloudy or non-cloudy points
        %Here we use a non-cloudy sounding to see if air parcels rising from this are consistent
        %with the LWC and equivalent potential temperature obtained
        %So, we follow the dry adiabat from this background profile until reach saturation
        T0=interp1(Zsound,Tsound,H_cb)+273.15; %find the values on the background sounding for given H_cb
        P0=interp1(Zsound,Psound,H_cb)*100;
        Q0=interp1(Zsound,qfp_sound,H_cb);        
        pot0 = T0*(1000e2/P0).^0.286;
        
        val_found=0;        
        pvals = P0:-0.1e2:100e2; %reduce the pressure until reach saturation
        for ip=1:length(pvals)
            t_dew=Tdew(Q0,pvals(ip)); %dew point for this vapour MR and pressure
            Tp = pot0./(1000e2./pvals(ip)).^0.286; %find the temperature (assuming constant potemp)
            if abs(Tp-t_dew)<0.1      %if have reached saturation
                val_found=1;
                break
            end
        end
        if val_found==0
            disp('*** cb value not found ***');
            break
        end
        tb = Tp;
        cb = pvals(ip)/100;
        

        
end

%tb=interp1(Zsound,Tsound,H_cb)+273.15;
%cb=interp1(Zsound,Psound,H_cb);

%tb=-12+273.15;
%cb=650;

fprintf(1,'\nH_CB=%.0f, T_CB=%.1f, P_CB=%.1f\n',H_cb,tb-273.15,cb);


                        pressure=dat_flt(inds,col_press)*1e2;
                        %really could just create a vector of monotonic pressure values


                        clear xline Tad
                        for iad=1:length(pressure)
                            [xline(iad),Tad(iad)]=adLWC_PaulLawson_simple(cb*1e2,tb,pressure(iad));
                        end
                        
                        TK=Tad+273.15;
                        pot=(TK).*(1000e2./pressure').^0.286;
                        totw = SatVapPress(tb,'goff','liq',cb*100,1)/f;
                        rho = density(pressure',Tad+273.15);
                        qv = totw-xline.*rho/1000; %calculate the vapour content of the moist adiabatic air
                        equiv = ( (TK + 2.453e6*qv/1004).*(1e5./pressure').^0.286 );
                        
                        %AMS definition - is pretty similar
                        qsat=SatVapPress(TK,'goff','liq',pressure',1)/f;   
                        RH = qv./qsat;
                        equiv2 = TK.*(1e5./pressure').^0.286 .* RH.^(-qv*4.615e2/1005.7)...
                            .* exp( 2.501e6*qv./TK/1005.7 );
                        
                        
                        switch vertical_coord
                            case 'temperature'
                                yline = dat_flt(inds,col_temp); %if plotting with temp as vertical
                            case 'pressure'
                                yline = dat_flt(inds,col_press); %if plotting with temp as vertical
                            case 'height'
                                yline = dat_flt(inds,col_alt)/1e3;
                        end
                        
                    case 'choose cb get temp from flight data'
                        %                tb=273.15-3;
                        %                tb=273.15-4.4;
                        %                tb=273.15-3;
                        cb=880;
                        cb=872;
                        cb=910;  %flight 122
                        %                cb=980;  %flight 120
                        cb=970;  %flight 113
                        cb=630;  %flight 105
                        cb=825;  %flight 99 upper layer
                        %                pressure=[cb:-10:500]*1e2;
                        %                         cb=940;
                        %                         cb=915;  %flight 101
                        %                         cb=715;  %flight 102
                        %                         cb=700;  %flight 102

                        pressure=dat_flt(:,col_press)*1e2;

                        ipress = find(pressure<(cb+2.5)*1e2 & pressure>(cb-2.5)*1e2);
                        temp = mean(dat_flt(ipress,col_temp));
                        %                         temp=-15.4;
                        tb=273.15+temp;



                        clear xline
                        for iad=1:length(pressure)
                            [xline(iad),Tad(iad)]=adLWC_PaulLawson_simple(cb*1e2,tb,pressure(iad));
                        end

                        yline = dat_flt(:,col_alt)/1e3;


                    case 'assume temperature gradient with height'
                        %choose a profile within the data to extend to ground
                        times=[14.81 14.835];
                        if length(times)==0
                            inds=1:length(dat_flt(:,1));
                        else
                            [time_0,time_1]=findheight_nearest(dat_flt(:,1)/1000/3600,times(1),times(2));
                            inds=time_0:time_1;
                        end

                        temp_gradient = (dat_flt(inds(1),col_temp)-dat_flt(inds(end),col_temp))/...
                            (dat_flt(inds(1),col_alt)-dat_flt(inds(end),col_alt));

                        cbz=100; %choose a cloud base height
                        cbz=0; %choose a cloud base height
                        %                         cbz=200; %choose a cloud base height

                        [H_lower ilower]=min(dat_flt(inds,col_alt));
                        [H_upper iupper]=max(dat_flt(inds,col_alt));
                        %                         T_lower=dat_flt(inds(end),col_temp)+273.16;
                        %                         T_upper=dat_flt(inds(1),col_temp)+273.16;

                        T_lower=dat_flt(inds(ilower),col_temp)+273.16;
                        T_upper=dat_flt(inds(iupper),col_temp)+273.16;
                        P_lower=dat_flt(inds(ilower),col_press)*1e2;

                        z2=[0:1:H_upper];
                        temp2=[T_upper+(H_upper-z2')*temp_gradient]';

                        HSPAN=[H_lower 0]; %pressure range for integration
                        [H,P] = ODE45(@hydrostatic,HSPAN,P_lower,[],z2,temp2); %enter the initial height for the given the first value of PSPAN

                        %function [F]=hydrostatic(z,p,zarr,Tarr)
                        %solves hydrostatic balance equation for pressure from starting pressure p0 and height z0 (m) and temperature profile T (in K)

                        cb = interp1(H,P,cbz);
                        tb = interp1(z2,temp2,cbz);

                        %now get the pressure vs. height for the whole profile
                        times=[14.7075 14.835];
                        if length(times)==0
                            inds=1:length(dat_flt(:,1));
                        else
                            [time_0,time_1]=findheight_nearest(dat_flt(:,1)/1000/3600,times(1),times(2));
                            inds=time_0:time_1;
                        end

                        p_profile = dat_flt(inds,col_press)*1e2;
                        h_profile = dat_flt(inds,col_alt);

                        ip=find(p_profile>=P(1)); %P(1) is the pressure at the bottom of the profile
                        p_profile(ip)=[];
                        h_profile(ip)=[];
                        ih=find(h_profile<=H(1));
                        h_profile(ih)=[];
                        p_profile(ih)=[];



                        p_combined = [p_profile; P];
                        h_combined = [h_profile; H];

                        [p_combined isort]=sort(p_combined);
                        h_combined = h_combined(isort);

                        [p_combined,iuni]=unique(p_combined);
                        h_combined = h_combined(iuni);

                        pressure = [cb:-1:min(p_combined)];
                        yline = interp1(p_combined,h_combined,pressure)/1e3;

                        clear xline
                        for iad=1:length(pressure)
                            [xline(iad),Tad(iad)]=adLWC_PaulLawson_simple(cb,tb,pressure(iad));
                        end



                end



        end
        end


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% switch scatter_case %%%
    case 'LWC PADS'
        plot_type='Dan LWC vs hotwire LWC';
        plot_type='PADS LWC vs hotwire LWC';
        %        plot_type='PADS LWC vs Dan LWC';
        %        plot_type='PADS hotwire LWC vs Dan hotwire LWC';

        switch plot_type
            case 'Dan LWC vs hotwire LWC'


                %CAS LWC calculated using my script from our CAS
                xdat(1).x = LWC_dist_PACS_cutoff;

                %hotwire LWC - as calculated from my script - think the spreadsheet ones may be wrong as they
                %don't get the same answer as Darrel Baumgardner's formulae
                ydat(1).y = interp1(PACS_TAS_time,LWC',PACS_time);
                ydat(1).y(ydat(1).y<-0.05)=NaN;

                title_name = ['PADS flight 433'];
                xlabel_name = ['LWC CAS size dist ' num2str(LWC_PACS_cutoff(1)) '-' num2str(LWC_PACS_cutoff(2)) ' \mum (g m^{-3})'];
                ylabel_name = 'LWC hotwire (g m^{-3})';

                ixlims=1;
                iylims=1;

                xlims = [0 2.2];
                ylims = [0 2.2];

                ione_to_one_line=1; %draw a one-to-one line

                icolour = 0;
                %                colour_data = CAS_mean_diameter;

                iadd_line=1;

                nanx = isnan(xdat(1).x);
                nany = isnan(ydat(1).y');
                not_nan = find(nanx==0 & nany==0);

                P = polyfit(xdat(1).x(not_nan),ydat(1).y(not_nan)',1);

                xline = [0:xlims(2)/50:xlims(2)];
                yline = P(1)*xline + P(2);

                if iadd_line==1
                    title_name = [title_name ' ' num2str(P(1)) 'x+' num2str(P(2))];
                end


            case 'PADS LWC vs hotwire LWC'

                %CAS LWC
                %                xdat(1).x = LWC_dist_cas;
                xdat(1).x = (1/1)*PACS_LWC_CAS;

                %hotwire LWC
                ydat(1).y = interp1(PACS_TAS_time,PACS_LWC,PACS_time);
                ydat(1).y(ydat(1).y<-0.05)=NaN;

                title_name = ['PADS flight 433'];
                xlabel_name = ['LWC CAS size dist ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum (g m^{-3})'];
                ylabel_name = 'LWC hotwire (g m^{-3})';

                ixlims=1;
                iylims=1;

                xlims = [0 2.2];
                ylims = [0 2.2];

                ione_to_one_line=1; %draw a one-to-one line

                icolour = 0;
                %                colour_data = CAS_mean_diameter;

                iadd_line=1;

                nanx = isnan(xdat(1).x);
                nany = isnan(ydat(1).y);
                not_nan = find(nanx==0 & nany==0);

                P = polyfit(xdat(1).x(not_nan),ydat(1).y(not_nan),1);

                xline = [0:xlims(2)/50:xlims(2)];
                yline = P(1)*xline + P(2);

                if iadd_line==1
                    title_name = [title_name ' ' num2str(P(1)) 'x+' num2str(P(2))];
                end

            case 'PADS LWC vs Dan LWC'
                %for comparing my calculation of LWC with that from the PADS spreadsheet
                %my LWC is about 4% higher - maybe to do with the method for calculating
                %the mean diameter for each bin - I use the mid-point but maybe then
                %do something different - e.g. assume constant gradient within bins?

                %PADS LWC
                %                xdat(1).x = LWC_dist_cas;
                xdat(1).x = (1/1)*PACS_LWC_CAS;

                %Dan LWC - using the CAS counts and the cas_sample_volume_and_stats.m script
                ydat(1).y = LWC_dist_PACS_cutoff';

                title_name = ['PADS flight 433'];
                xlabel_name = ['LWC PADS CAS size dist (g m^{-3})'];
                ylabel_name = ['LWC CAS size dist ' num2str(CAS_LWC_cut_off_sizes(1)) '-' num2str(CAS_LWC_cut_off_sizes(2)) ' \mum (g m^{-3})'];

                ixlims=1;
                iylims=1;

                xlims = [0 2.2];
                ylims = [0 2.2];

                ione_to_one_line=1; %draw a one-to-one line

                icolour = 0;
                %                colour_data = CAS_mean_diameter;

                iadd_line=1;

                nanx = isnan(xdat(1).x);
                nany = isnan(ydat(1).y);
                not_nan = find(nanx==0 & nany==0);

                P = polyfit(xdat(1).x(not_nan),ydat(1).y(not_nan),1);

                xline = [0:xlims(2)/50:xlims(2)];
                yline = P(1)*xline + P(2);

                if iadd_line==1
                    title_name = [title_name ' ' num2str(P(1)) 'x+' num2str(P(2))];
                end

            case 'PADS hotwire LWC vs Dan hotwire LWC'
                %for comparing my calculation of the hotwire LWC with that from the PADS spreadsheet
                %my LWC is about 14.47% higher - probably because of the different FACT
                %is 13.7% higher when using contstant boiling point of 373.16 K.

                %from the spreadsheet
                xdat(1).x = PACS_LWC;

                %calculated by my script
                ydat(1).y = LWC_calculated;

                title_name = ['PADS flight 433'];
                xlabel_name = 'PADS LWC hotwire (g m^{-3})';
                ylabel_name = 'Dan''s calculated LWC hotwire (g m^{-3})';

                ixlims=1;
                iylims=1;

                xlims = [0 2.2];
                ylims = [0 2.2];

                ione_to_one_line=1; %draw a one-to-one line

                icolour = 0;
                %                colour_data = CAS_mean_diameter;

                iadd_line=1;

                nanx = isnan(xdat(1).x);
                nany = isnan(ydat(1).y);
                not_nan = find(nanx==0 & nany==0);

                P = polyfit(xdat(1).x(not_nan),ydat(1).y(not_nan),1);

                xline = [0:xlims(2)/50:xlims(2)];
                yline = P(1)*xline + P(2);

                if iadd_line==1
                    title_name = [title_name ' ' num2str(P(1)) 'x+' num2str(P(2))];
                end

                figname=plot_type;




        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

if i_override_vertical_coord==1
    switch vertical_coord
        case 'height'
            ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),time_base_scatter)'/1e3;
            ylabel_name = 'Altitude (km)';
%            iylims=0;
%            ylims = [2.5 3.5];
        case 'temperature'
            ydat(1).y = interp1(dat_flt(:,1)/1e3,temperature_data,time_base_scatter)';
            ylabel_name = 'Temperature (^oC)';
            iy_reverse=1;
            iylims=1;
            ylims = [-16 -7];
            ylims = [-18 1]; %flight 99              
%            ylims = [-7 1]; %flight 100  
%            ylims = [-23 1]; %flight 101              
        case 'pressure'
            ydat(1).y = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
            ylabel_name = 'Pressure (hPa)';
            iy_reverse=1;
        case 'potemp'
            ylabel_name = 'Potential Temperature (K)';
            T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),time_base_scatter)';
            P=interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),time_base_scatter)';
            ydat(1).y = T.*(1000./P).^0.286;
            iy_reverse=0;


    end
    
end

if ilimit_data==1
    xdat(1).x=xdat(1).x(inds_limit);
    ydat(1).y=ydat(1).y(inds_limit);
    if icolour==1
        colour_data=colour_data(inds_limit);    
    end
    
    disp('***** WARNING - limiting the data to a sub-set !!! ********');
end

scrsz=get(0,'ScreenSize');
%posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.13];
%posit=[9 50 scrsz(3)/1.9 scrsz(4)/2.1];  %used 22nd Jan, 2009
posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];

figname=[title_name ' scatter plot ' xlabel_name ' vs ' ylabel_name];
figname=[title_name ' scatter plot '];

hf=figure('name',figname,'Position',posit);

if ijoin_points==1
    plot(xdat(1).x,ydat(1).y,'b-','linewidth',2);
    hold on
end

if icolour==1
    hscatter=scatter(xdat(1).x,ydat(1).y,100,colour_data,'o','filled');
    colorbar;
else
%    scatter(xdat(1).x,ydat(1).y,10,'k');
   eval( ['scatter(xdat(1).x,ydat(1).y,' scatter_properties_string ');'] );
end

hold on;


if iadd_line==1
    plot(xline,yline,'k');
end

if ixlims==1
    set(gca,'xlim',xlims);
end
if iylims==1
    set(gca,'ylim',ylims);
end

if ione_to_one_line==1
    ylim_scatter=get(gca,'ylim');
    xlim_scatter=get(gca,'xlim');
    if(strcmp(get(gca,'ydir'),'reverse'))
        ytemp=ylim_scatter(1);
        ylim_scatter(1)=ylim_scatter(2);
        ylim_scatter(2)=ytemp;
    end
    if(strcmp(get(gca,'xdir'),'reverse'))
        xtemp=xlim_scatter(1);
        xlim_scatter(1)=xlim_scatter(2);
        xlim_scatter(2)=xtemp;
    end
    
    L1=min(xlim_scatter(1),ylim_scatter(1));
    L2=min(xlim_scatter(2),ylim_scatter(2));    
    
    line([L1 L2],[L1 L2]);  %    line([0 10],[0 10]);
end

if iclims==1
    set(gca,'clim',clims);
end

if iset_cmap==1
    colormap(cmap);
end

if iy_reverse==1
    set(gca,'ydir','reverse');
end

inan=find(isnan(xdat(1).x)==0 & isnan(ydat(1).y)==0 );

XX = xdat(1).x(inan);
YY = ydat(1).y(inan);

corrXY = corr(XX,YY)
rmse = sqrt ( mean ( (XX-YY).^2 ) )

ylabel(ylabel_name,'fontsize',fsize);
xlabel(xlabel_name,'fontsize',fsize);
title_name = [title_name ' corr=' num2str(corrXY) ', rmse=' num2str(rmse)];
titlenamw=textwrap({title_name},100);
title(titlenamw,'fontsize',fsize);
grid;
set(gca,'fontsize',fsize);

nmax=200;
if length(figname)>nmax
    savename=figname(1:nmax);
else
    savename=figname;
end

savename=[savedir savename];


if iadd_nums_above==1
    add_numbers_above_timeseries  %add numbers above the points of a timeseries
end



if iplot_error_bar==1
    minX=9e99; maxX=-9e99;
for idat=1:length(xdat)
    minX = min([min(xdat(idat).x) minX]);
    maxX = max([max(xdat(idat).x) minX]);
end
spanX = maxX-minX;


%    hE=herrorbar(xdat(1).x, ydat(1).y, error_barH(1).dat);
     hEH = errorbarYY('horiz',xdat(1).x,ydat(1).y,error_barH(1).dat,gca,'b','none',1,0.01);
%    unplot
%    hE=errorbar(xdat(1).x, ydat(1).y, error_barV(1).dat);
%    unplot
     hEV = errorbarYY('vert',xdat(1).x,ydat(1).y,error_barV(1).dat,gca,'b','none',1,0.04,spanX);
end


if ihighlightX==1
    for i=1:length(highlight_pointsX)
        [minval imin] = min(abs(highlight_pointsX(i) - xdat(1).x));
        imin
        plot(xdat(1).x(imin),ydat(1).y(imin),'ks','markerfacecolor','k');
    end
elseif ihighlightY==1
    for i=1:length(highlight_pointsY)
        [minval imin] = min(abs(highlight_pointsY(i) - ydat(1).y));
        plot(xdat(1).x(imin),ydat(1).y(imin),'ks','markerfacecolor','k');
    end
end



if iaxis_square==1
    axis square
end