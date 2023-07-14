amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/HADSST_matched_to_amsre_ssts_global_1990-Sep2020.mat';
load(amsre_matfile,'gcm_Plat2D_AMSRE','gcm_Plon2D_AMSRE','year_amsre','month_amsre','day_amsre','sst_amsre_smooth');
%Gets the sst_amsre_time3 SSTs for the times in daynum_timeseries3 from sst_amsre_smooth
clear ioverride_block_amsre
amsre_block_average_time



%cf = Cloud_Fraction_Liquid.data;
%WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
%tau = Cloud_Optical_Thickness_Liquid_Mean.data;
%                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres

%try

% cf_time=Cloud_Fraction_Liquid.timeseries(:,:);
% NP_time=Cloud_Fraction_Liquid_Pixel_Counts.timeseries(:,:)./cf_time;
%
% tau_time = Cloud_Optical_Thickness_Liquid_Mean.timeseries(:,:);
% reff_time = Cloud_Effective_Radius_Liquid_Mean.timeseries(:,:)*1e-6; %convert to metres
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
%
% [N_time,H,W,k,Q,cw]=MODIS_N_H_func(tau_time,reff_time,Wflag,0);

switch mod_data_type
    case {'timeseries3','timeseries3 lambert'}
        
        daynum_timeseries3_MODIS = daynum_timeseries3;

        Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
        %setting last argument to 'H', 'W', or 'N'
        %tells it what to output - only doing one at a time to save memory
        [H_time3]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'H');
        [W_time3]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'W');
        [N_time3]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');        

%        [cw_time3]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'cw');        


if exist('Cloud_Effective_Radius_37_Liquid_Mean')
    Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
    %setting last argument to 'H', 'W', or 'N'
    %tells it what to output - only doing one at a time to save memory
    [H_time3_37]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_37_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'H');
    [W_time3_37]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_37_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'W');
    [N_time3_37]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_37_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');
end

%Cloud_Optical_Thickness_Liquid_Mean_Uncertainty and Re uncertainties are
%in percent
        make_un=1;
        if make_un==1 & exist('Cloud_Optical_Thickness_Liquid_Mean_Uncertainty') & exist('Cloud_Effective_Radius_Liquid_Mean_Uncertainty')
            N_un_time3 = N_time3/100 .* sqrt( (0.5*Cloud_Optical_Thickness_Liquid_Mean_Uncertainty.timeseries3).^2 + (-5/2*Cloud_Effective_Radius_Liquid_Mean_Uncertainty.timeseries3).^2 );
        end
        
        if exist('Cloud_Optical_Thickness_Liquid_QA_Mean')
            [N_QA_time3]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_QA_Mean.timeseries3,Cloud_Effective_Radius_Liquid_QA_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');
        end
        
        %calculating the normalised (i.e. divided by mean W) std dev of W
        %*** assuming no covariance ***
         if exist('Cloud_Effective_Radius_Liquid_Standard_Deviation')
             stdW_time3 = sqrt(  (Cloud_Optical_Thickness_Liquid_Standard_Deviation.timeseries3./Cloud_Optical_Thickness_Liquid_Mean.timeseries3).^2 + (Cloud_Effective_Radius_Liquid_Standard_Deviation.timeseries3./Cloud_Effective_Radius_Liquid_Mean.timeseries3).^2 );
         end
         
        if exist('Cloud_Water_Path_Liquid_Standard_Deviation')
            %convert back to kg/m2 for consistency with W_time3
            Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3 = Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3/1e3;
        end
        
        ihomog=0;        
        %convert back to kg/m2 for consistency with W_time3
        if exist('Cloud_Water_Path_Liquid') & isnan(maxALL(Cloud_Water_Path_Liquid.timeseries3))==0
            Cloud_Water_Path_Liquid.timeseries3 = Cloud_Water_Path_Liquid.timeseries3/1e3;
            homog_time3_W = (Cloud_Water_Path_Liquid.timeseries3./Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3).^2;
%            homog_time3 = homog_time3_W; ihomog=1;
            homog_cutoff = 1000;
            homog_time3_W(homog_time3_W>homog_cutoff) = homog_cutoff;
            homog_time3 = homog_time3_W;
            ihomog=1;
        end
        if exist('Cloud_Water_Path_Liquid_Log_Mean') & isnan(maxALL(Cloud_Water_Path_Liquid_Log_Mean.timeseries3))==0
            Cloud_Water_Path_Liquid_Log_Mean.timeseries3 = Cloud_Water_Path_Liquid_Log_Mean.timeseries3 - log(1e3);            
             %for Cahalan's homogeneity factor = exp(mean(lnW))./ mean(W)   
             homog_time3_logW_W = exp(Cloud_Water_Path_Liquid_Log_Mean.timeseries3) ./ Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3;
%             homog_time3 = homog_time3_logW; ihomog=1;
        end
         
        if exist('Cloud_Water_Path_Liquid_Standard_Deviation') & isnan(maxALL(Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3))==0
            homog_time3_meanW = (W_time3./Cloud_Water_Path_Liquid_Standard_Deviation.timeseries3).^2;                       
            homog_cutoff = 1000;
            homog_time3_meanW(homog_time3_meanW>homog_cutoff) = homog_cutoff;
             if ihomog==0
                homog_time3=homog_time3_meanW;
             end
        else
           homog_time3  = 1e8*ones(size(Cloud_Top_Temperature_Day_Mean.timeseries3)); 
        end
        
        if exist('Cloud_Optical_Thickness_Liquid_Log_Mean') & isnan(maxALL(Cloud_Optical_Thickness_Liquid_Log_Mean.timeseries3))==0
            homog_tau_Cahalan = exp(Cloud_Optical_Thickness_Liquid_Log_Mean.timeseries3)./Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
        end   %so this is homogeneity - values closer to one mean more homgeneity.
        
        if exist('sst_amsre_time3') & isnan(maxALL(sst_amsre_time3))==0
            CTH.timeseries3 = (273.15 + sst_amsre_time3 - Cloud_Top_Temperature_Day_Mean.timeseries3 - 2.35) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
            CTH_adjust = (273.15 + sst_amsre_time3 - Cloud_Top_Temperature_Day_Mean.timeseries3 + 2.5) / 0.0069 /1e3; %CTH from Zuidema (2009) in km            
        elseif exist('sst_amsre_time3')
            disp('***WARNING - SST data all NaN - setting CTH to 1e8 ***')
%            CTH.timeseries3 = zeros(size(sst_amsre_time3));
            CTH.timeseries3 = 1e8*ones(size(Cloud_Top_Temperature_Day_Mean.timeseries3));            
        else
            disp('***WARNING - no SST data - setting CTH to 1e8 ***')
%            CTH.timeseries3 = zeros(size(sst_amsre_time3));
            CTH.timeseries3 = 1e8*ones(size(Cloud_Top_Temperature_Day_Mean.timeseries3));                        
            
        end
        
        if exist('sst_amsre_time3') & isnan(maxALL(sst_amsre_time3))==0 & exist('Cloud_Top_Temperature_Day_Minimum')
            CTH_max.timeseries3 = (273.15 + sst_amsre_time3 - Cloud_Top_Temperature_Day_Minimum.timeseries3 - 2.35) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
            imax_CTH=1;
        elseif exist('sst_amsre_time3') & exist('Cloud_Top_Temperature_Day_Minimum')
            disp('***WARNING - SST data all NaN - setting CTH to 1e8 ***')
%            CTH.timeseries3 = zeros(size(sst_amsre_time3));
            CTH_max.timeseries3 = 1e8*ones(size(Cloud_Top_Temperature_Day_Mean.timeseries3));  
            imax_CTH=0;
        else
            disp('***WARNING - either no SST data or no min CTT data- setting CTH to 1e8 ***')
%            CTH.timeseries3 = zeros(size(sst_amsre_time3));
            CTH_max.timeseries3 = 1e8*ones(size(Cloud_Top_Temperature_Day_Mean.timeseries3));  
             imax_CTH=0;
        end
        
        if exist('sst_amsre_time3') & isnan(maxALL(sst_amsre_time3))==0 & exist('Cloud_Top_Temperature_Day_ice_liq_Mean')
            CTH_all.timeseries3 = (273.15 + sst_amsre_time3 - Cloud_Top_Temperature_Day_ice_liq_Mean.timeseries3 - 2.35) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
        elseif exist('sst_amsre_time3')
            disp('***WARNING - no Cloud_Top_Temperature_Day_ice_liq_Mean data - setting CTH_all to zero ***')
            disp('... but CTH is ok (just CTH_all affected)');
            CTH_all.timeseries3 = zeros(size(sst_amsre_time3));
        end
        
        if exist('homog_time3') & isnan(maxALL(homog_time3))==1
            disp('*** WARNING - homog_time3 is all NaNs ***');
        end
        
        
        if imax_CTH==1 & exist('Cloud_Top_Pressure_Day_Minimum')
            CTH_max_hybrid.timeseries3 = CTH_max.timeseries3;
            CTP_min_hybrid.timeseries3 = NaN*ones(size(Cloud_Top_Pressure_Day_Minimum.timeseries3)); %set as NaN and populate later
            i=find(CTH_max.timeseries3>2); %find points where max CTH from CTT method is > x km
%            i=find(Cloud_Top_Pressure_Day_Minimum.timeseries3<800); %find points where min CTP is < threshold
            CTH_max_hybrid.timeseries3(i) = NaN;
            CTP_min_hybrid.timeseries3(i) = Cloud_Top_Pressure_Day_Minimum.timeseries3(i);   
            
            %So, where the CTH is > x km the CTH_max_hybrid array  is NaN,
            %but the CTP_min_hybrid array is populated and vice versa
            %Then the screening will allow matches based on either of these
            
            file_SA = '/home/disk/eos1/d.grosvenor/standard_atmos.mat'; load(file_SA);
            CTH_hybrid_std_atmos_CTP = CTH_max_hybrid;
            
            %Get an approx CTH from the min CTP 
            CTH_hybrid_std_atmos_CTP.timeseries3(i) = interp1(SA.P/1e2,SA.z,Cloud_Top_Pressure_Day_Minimum.timeseries3(i));
        end
        
        %timeseries3
%        cf_time3=Cloud_Fraction_Liquid.timeseries3;
%        NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3./cf_time3;       


        if exist('Solar_Zenith_Mean')
%            SZA_time3=Solar_Zenith_Mean.timeseries3;
        end
        if exist('Sensor_Zenith_Mean')
%            sensZA_time3 = Sensor_Zenith_Mean.timeseries3;
        end
        if exist('Solar_Zenith_Maximum')
%            maxSZA_time3 = Solar_Zenith_Maximum.timeseries3;
        end
        if exist('Solar_Zenith_Minimum')
%            minSZA_time3 = Solar_Zenith_Minimum.timeseries3;
        end
        if exist('Sensor_Zenith_Maximum')
%            max_sensZA_time3 = Sensor_Zenith_Maximum.timeseries3;
        end

%        tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
%        reff_time3 = Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6; %convert to metres
%        CTT_time3 = Cloud_Top_Temperature_Day_Mean.timeseries3;
       

%        [N_time3,H_time3,W_time3,k,Q,cw_time3]=MODIS_N_H_func(tau_time3,reff_time3,Wflag,0,CTT_time3);

        %catch
        %    disp('*** warning - missed some calculations due to errors ***');
        %end
        
        if exist('Sensor_Azimuth_Mean') & exist('Solar_Azimuth_Mean')
%             Relative_Azimuth_Mean.timeseries3=(Sensor_Azimuth_Mean.timeseries3-Solar_Azimuth_Mean.timeseries3);
%             Relative_Azimuth_Mean.timeseries3=abs(Relative_Azimuth_Mean.timeseries3);
%             i180=find(Relative_Azimuth_Mean.timeseries3>180);
%             Relative_Azimuth_Mean.timeseries3(i180)=360-Relative_Azimuth_Mean.timeseries3(i180);
            
            %Relative azimuth is defined as zero when the sensor is looking
            %into the Sun - as if there as no scattering (just
%             %transmission). Backscatter is then at 180 degrees. So do the
%             difference of sensor and solar and subtract 180 so that if
%             the difference is 180 (forward scatter) then will get relAZ=0
%NOTE - %with Joint files only get the relative AZ, but have kept the variable
%Solar_Azimuth_Angle for consistency with full L2 and have set the sensorAZ to be the
%relative AZ. But are the MODIS Relative_Azimuth_Angles with the 180
%already removed? I.e. do they follow the above convention? Check using L2?
%Yes, if calculate relAz as below using L2 then get the same answer as that
%given by Relative_Azimuth from the jointL2 product.

if maxALL(Solar_Azimuth_Mean.timeseries3)==0 %if have no solarAz then the relAz should be in the SensorAz
    Relative_Azimuth_Mean.timeseries3=Sensor_Azimuth_Mean.timeseries3;
else %otherwise calculate it.
    Relative_Azimuth_Mean.timeseries3=calc_relAz(Sensor_Azimuth_Mean.timeseries3,Solar_Azimuth_Mean.timeseries3);
end

if exist('Cloud_Effective_Radius_16_Liquid_Mean')
%maximum difference between 1.6, 2.1 and 3.7 mum bands
reffs=cat(4,Cloud_Effective_Radius_16_Liquid_Mean.timeseries3,Cloud_Effective_Radius_Liquid_Mean.timeseries3);
reffs=cat(4,reffs,Cloud_Effective_Radius_37_Liquid_Mean.timeseries3);
dreffs = max(reffs,[],4) - min(reffs,[],4);
rel_dreffs = 100* dreffs ./ (mean(reffs,4));
end




        end

    case 'daily'


        %calculate N and H from MOD06 optical thickness, effective radius and cloud fraction

        use_QA_vals=0;

        if use_QA_vals==1
            %using MOD06 only here

            fprintf(1,'\n**** USING QA means for tau and reff ***\n');
            title_info = 'QA Means';


            tau_time3 = Cloud_Optical_Thickness_Liquid_QA_Mean.data;
            %tau(tau>50)=NaN; %remove large tau values
            %fprintf(1,'\n**** WARNING - am removing data when tau >50 ***\n');

            reff_time3 = Cloud_Effective_Radius_Liquid_QA_Mean.data*1e-6; %convert to metres
            WMOD=Cloud_Water_Path_Liquid_QA_Mean.data/1000; %convert to kg/m2
        else
            %using MOD06 only here

            tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.data;
            %tau(tau>50)=NaN; %remove large tau values
            %fprintf(1,'\n**** WARNING - am removing data when tau >50 ***\n');

            reff_time3 = Cloud_Effective_Radius_Liquid_Mean.data*1e-6; %convert to metres
            %    WMOD=Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
            WMOD=NaN;
        end
        
        
         if exist('Solar_Zenith_Mean')
            SZA_time3=Solar_Zenith_Mean.data;
        end
        if exist('Sensor_Zenith_Mean')
            sensZA_time3 = Sensor_Zenith_Mean.data;
        end
        if exist('Solar_Zenith_Maximum')
            maxSZA_time3 = Solar_Zenith_Maximum.data;
        end
        if exist('Solar_Zenith_Minimum')
            minSZA_time3 = Solar_Zenith_Minimum.data;
        end
        if exist('Sensor_Zenith_Maximum')
            max_sensZA_time3 = Sensor_Zenith_Maximum.data;
        end
        
        


            
        sangle = Scattering_Angle_Mean.data;
        cf_time3 = Cloud_Fraction_Liquid.data;
        NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.data./cf_time3;         
        CTT_time3 = Cloud_Top_Temperature_Day_Mean.data;

        if ~exist('set_MODIS_NH_flags') | set_MODIS_NH_flags==0

            Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
            %Wflag='MODIS'; %use the MODIS LWP

        else
            clear set_MODIS_NH_flags
        end

        [N,H,W,k,Q,cw]=MODIS_N_H_func(tau,reff,Wflag,WMOD,CTT_time3);



end


