%select the timeseries type required
% timeseries_case=''; %for timeseries
% timeseries_case='3'; %for timeseries3

x_axis_vals = 'sensor_ZA'; 
x_axis_vals = 'Nd L2 timeseries'; 
x_axis_vals = 'Nd L2 swath'; 
%x_axis_vals = 'Reff'; 
%x_axis_vals = 'Tau'; 
%x_axis_vals = 'Reff L2 swath';


y_axis_vals='solar_ZA'; data_type='timeseries';
%y_axis_vals = 'sensor_ZA'; data_type='timeseries';
%y_axis_vals = 'Reff'; data_type='timeseries';
%y_axis_vals = 'Tau'; data_type='timeseries';
y_axis_vals = 'Latitude'; data_type='timeseries';
y_axis_vals = 'Nd L2 swath'; data_type='swath';
%y_axis_vals = 'Nd Uncertainty L2 swath'; data_type='swath';
%y_axis_vals = 'Tau L2 swath'; data_type='swath';
%y_axis_vals = 'Reff L2 swath'; data_type='swath';
%y_axis_vals = 'Reff diff L2 swath 6-7'; data_type='swath';
%y_axis_vals = 'Reff diff L2 swath 20-7'; data_type='swath';



%select the number of PDF bins
nXpdf=200;
%nXpdf=1000;
nYpdf=100;

%nXpdf=10;
%nYpdf=10;


%set some default values (leaving these unchanged will result in the bins
%going from the min to the max value in nXpdf bins
minXbins=-9e99;
maxXbins=-9e99;

minYbins=-9e99;
maxYbins=-9e99;


%select the LAT and LON values to make up the
%points included
% LAT_val = [-0.5000];
% LAT_val = [-60.5000];
% %LAT_val = [-59.5000];
% %LAT_val = [-89.5:1:89.5];
% LAT_val = [-30.5:1:30.5];
% %LAT_val = [-12.5000];
% % LAT_val = [72];  %72
% % %LAT_val = [62];
% % LAT_val = [0:2:90];
% % 
% % LAT_val = [72 73];     
% % 
LAT_val = [72.5 73.5];



% 
 
% 
% LON_val = [-179.5:179.5];
% % LON_val = -156.75;
% % LON_val = -153.0000;
% % LON_val = -149.0000;  %
% % LON_val = [-156.75 -153 -149];
% % % 
% % 
% % LON_val = [-152.5 -150];
% % LON_val = [-150 -147.5];
% % 
LON_val = [-149.5 -148.5];

%ilat=1;
%ilon=1;

switch data_type
    case 'timeseries'
        LAT=timLAT;
        LON=timLON;
        
        ilat=findheight_nearest(LAT,LAT_val);
        ilon=findheight_nearest(LON,LON_val);
    case 'swath'
        filtering_data_L2
%will need to find the points between the values for swaths as the required
%area will not be a rectangle
%These should be set in filtering_data_L2
%         row_L2_inds = [4:1349];
%         col_L2_inds = [4:2029];

end




% thresh=150;
% thresh=0.3;
% thresh=0;
% 
% cf = Cloud_Fraction_Liquid.data;
% WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
% tau = Cloud_Optical_Thickness_Liquid_Mean.data;
% %                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres
% 
% 
% 
% cf_time=Cloud_Fraction_Liquid.timeseries(ilat,:);
% NP_time=Cloud_Fraction_Liquid_Pixel_Counts.timeseries(ilat,:)./cf_time;
% 
% tau_time = Cloud_Optical_Thickness_Liquid_Mean.timeseries(ilat,:);
% reff_time = Cloud_Effective_Radius_Liquid_Mean.timeseries(ilat,:)*1e-6; %convert to metres
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
% 
% [N_time,H,W,k,Q,cw]=MODIS_N_H_func(tau_time,reff_time,Wflag,0);
% 
% 
% %timeseries3
% cf_time3=Cloud_Fraction_Liquid.timeseries3(ilat,ilon,:,:);
% NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3(ilat,ilon,:,:)./cf_time3;
% sensZA_time3 = Sensor_Zenith_Mean.timeseries3(ilat,ilon,:,:);
% 
% tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,:,:);
% reff_time3 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,:,:)*1e-6; %convert to metres
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
% 
% [N_time3,H,W,k,Q,cw]=MODIS_N_H_func(tau_time3,reff_time3,Wflag,0);


%%%%
extra_title_info = [' for day ' modis_day_str ', Y' modis_year_str];

switch data_type
    case 'timeseries'
        LAT2 = LAT(ilat);
        LON2 = LON(ilon);
    case 'swath'
        LAT2 = LAT_val;
        LON2 = LON_val;
end


if length(LAT2)==1
    LAT_str=num2str(LAT2);
else
    LAT_str=[num2str(LAT2(1)) ' to ' num2str(LAT2(end))];
end

if length(LON2)==1
    LON_str=num2str(LON2(1));
else
    LON_str=[num2str(LON2(1)) ' to ' num2str(LON2(end))];
end




switch x_axis_vals
    case 'Nd L2 swath other channel'
        xlabelstr = 'Nd (cm^{-3})';

        X = N(row_L2_inds,col_L2_inds);
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
%         minXbins=0;
%         maxXbins=500;
    
    case 'Reff L2 swath'
        xlabelstr = 'Re (\mum)';
        X = re(row_L2_inds,col_L2_inds);
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
    case 'Nd L2 swath'
        xlabelstr = 'Nd (cm^{-3})';

        X = N(row_L2_inds,col_L2_inds);
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        
%         minXbins=0;
%         maxXbins=500;
%         maxXbins=200;         
        
     case 'Reff'
        xlabelstr = 'Re (\mum)';
        X = squeeze(re_5km_tim(ilat,ilon,:));
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        scantime_matlab_tim_pdf = squeeze(scantime_matlab_tim(ilat,ilon,:));
    case 'Tau'
        xlabelstr = 'Optical Depth';
        X = squeeze(tau_5km_tim(ilat,ilon,:));
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        scantime_matlab_tim_pdf = squeeze(scantime_matlab_tim(ilat,ilon,:));    
    case 'Nd L2 timeseries'
        xlabelstr = 'Nd (cm^{-3})';
        Wflag='calc';
        [N,H,W,k,Q,cw]=MODIS_N_H_func(tau_5km_tim,re_5km_tim*1e-6,Wflag,NaN);
        X = squeeze(N(ilat,ilon,:));
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        scantime_matlab_tim_pdf = squeeze(scantime_matlab_tim(ilat,ilon,:));
    case 'sensor_ZA'
        
        xlabelstr = 'Sensor Zenith Angle';
        X = squeeze(sensor_zenith_tim(ilat,ilon,:));
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        scantime_matlab_tim_pdf = squeeze(scantime_matlab_tim(ilat,ilon,:));
        
    case 'solar_ZA'
        xlabelstr = 'Solar Zenith Angle';
        X = squeeze(solar_zenith_tim(ilat,ilon,:));
        ihtot = [1:prod(size(X))]; thresh_str='xxx'; %all the data
        scantime_matlab_tim_pdf = squeeze(scantime_matlab_tim(ilat,ilon,:));
        
    case 'Mean SZA timeseries3'
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:)))]; thresh_str='xxx'; %all the data
        xlabelstr = 'Mean SZA';
        X = Solar_Zenith_Mean.timeseries3(ilat,ilon,:);    
      
    case 'Nd from grid vals'
        xlabelstr = 'N_d (cm^{-3})';
        ihtot = [1:prod(size(N_histo_mean))]; thresh_str='xxx'; %all the data
        xdat(1).x = N;

    case 'Nd from grid vals timeseries'
        xlabelstr = 'N_d (cm^{-3}) from grid values timeseries';
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data
        xdat(1).x = N_time;

    case 'Nd from histogram vals timeseries'
        xlabelstr = 'N_d (cm^{-3}) from histogram timeseries';
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data

        [histo_output]=...
            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,ihtot) );
        xdat(1).x=histo_output.N_histo_mean;
        %                        xdat(1).x=histo_output.N_histo_std;

    case 'Nd from histogram vals timeseries3'
        xlabelstr = 'N_d (cm^{-3}) from histogram timeseries3';
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:)))]; thresh_str='xxx'; %all the data

%        [histo_output]=...
%            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );
        %                        xdat(1).x=histo_output.N_histo_mean;
        X = Nd_timeseries.mean(ilat,ilon,:);
        
%        minXbins=0;
%        maxXbins=500;

    case 'Reff'
        
        xlabelstr = 'Reff (\mum)';
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:)))]; thresh_str='xxx'; %all the data
        X = reff_time3*1e6;
        
   case 'Tau'
        
        xlabelstr = 'Tau';
        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:)))]; thresh_str='xxx'; %all the data
        X = tau_time3;        


    case 'Std. dev of Nd from histogram timeseries'

        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries(1,:)))]; thresh_str='xxx'; %all the data

        [histo_output]=...
            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,ihtot) );

        xlabelstr = 'Std dev N_d (cm^{-3}) from histogram timeseries';
        xdat(1).x=histo_output.N_histo_std;

        xlabelstr = 'Normalised std dev N_d from histogram timeseries';
        xdat(1).x=histo_output.N_std_norm;

    case 'Std. dev of Nd from histogram timeseries3'

        ihtot = [1:prod(size(Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:)))]; thresh_str='xxx'; %all the data

%        [histo_output]=...
%            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );

        xlabelstr = 'Std dev N_d (cm^{-3}) from histogram timeseries';
        X = Nd_timeseries.std_dev(ilat,ilon,:);
        
        
        xlabelstr = 'Normalised std dev N_d from histogram timeseries';
        X = Nd_timeseries.std_dev_norm(ilat,ilon,:);

        



end






switch y_axis_vals
    case 'Reff diff L2 swath 6-7'
        ylabelstr = 'Reff Difference Band 6 - Band 7 (\mum)';
        Y = re_diff(row_L2_inds,col_L2_inds,1);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str]; 
        
    case 'Reff diff L2 swath 20-7'
        ylabelstr = 'Reff Difference Band 20 - Band 7 (\mum)';
        Y = re_diff(row_L2_inds,col_L2_inds,2);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];     
        
    case 'Tau L2 swath'
        ylabelstr = 'Optical Depth';
        Y = tau(row_L2_inds,col_L2_inds);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];         
        
    case 'Reff L2 swath'
        ylabelstr = 'Reff (\mum)';
        Y = re(row_L2_inds,col_L2_inds);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];  
        
    case 'Nd Uncertainty L2 swath'
        ylabelstr = 'N_d Uncertainty (%)';
        Y = percent_error_Nd(row_L2_inds,col_L2_inds);
         

%        minYbins=0;
%        maxYbins=110;
        
%        minXbins=0;
%        maxXbins=500;



    case 'Nd L2 swath'
        ylabelstr = 'N_d (cm^{-3})';
        Y = N(row_L2_inds,col_L2_inds);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];    
        
%         minYbins=0;
%         maxYbins=500;
    case 'Latitude'
        ylabelstr = 'Latitude';
        Y = squeeze( repmat(timLAT(ilat),[size(tau_5km_tim,3) 1])' );
        Y = repmat(Y,[1 1 length(ilon)]);
        Y = permute(Y,[1 3 2]);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
     case 'Reff'
        ylabelstr = 'Re (\mum)';
        Y = squeeze(re_5km_tim(ilat,ilon,:));
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Tau'
        ylabelstr = 'Optical Depth';
        Y = squeeze(tau_5km_tim(ilat,ilon,:));
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];  
     case 'sensor_ZA'    
        ylabelstr = 'Sensor Zenith Angle';
        Y = squeeze(sensor_zenith_tim(ilat,ilon,:));
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        
    case 'solar_ZA'
        ylabelstr = 'Solar Zenith Angle';
        Y = squeeze(solar_zenith_tim(ilat,ilon,:));
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        
        iylim=1; 
        ylims=[0 90];  
        
    case 'WMOD'
        ylabelstr = 'W times five sixths (kg m^{-2})';
        Y = WMOD(:); %MOD06

    case 'CF'
        ylabelstr = 'Cloud Fraction';
        Y = cf(:); %MOD06 CF

    case 'Scattering angle'
        ylabelstr = 'Scattering angle (degrees)';
        Y = sangle(:);

    case 'SZA std dev timeseries'
        ylabelstr = 'Std dev of SZA';
        Y = Solar_Zenith_Standard_Deviation.timeseries(ilat,:);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Mean SZA timeseries'
        ylabelstr = 'Mean SZA';
        Y = Solar_Zenith_Mean.timeseries(ilat,:);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Mean SZA timeseries3'
        ylabelstr = 'Mean SZA';
        Y = Solar_Zenith_Mean.timeseries3(ilat,ilon,:);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
    case 'Mean Sensor Zenith Angle timeseries3'
        ylabelstr = 'Mean Sensor ZA';
        Y = Sensor_Zenith_Mean.timeseries3(ilat,ilon,:);
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        
        
    case 'SZA std dev timeseries2'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'Std dev of SZA';
        Y = Solar_Zenith_Standard_Deviation.timeseries(ilat,:);
        %                        xdat(1).x = N_time(ihtot);

    case 'SZA std dev timeseries3'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'Std dev of SZA';
        Y= Solar_Zenith_Standard_Deviation.timeseries3(ilat,ilon,:);
        %                        xdat(1).x = N_time(ihtot);

    case 'Nd from histogram vals timeseries'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'N_d from histo timeseries (cm^{-3})';

        [histo_output]=...
            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries(:,:,ilat,:) );
        Y=histo_output.N_histo_mean;

    case 'Nd from histogram vals timeseries3'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'N_d from histo timeseries3 (cm^{-3})';

%        [histo_output]=...
%            histo_mean_calc_MODIS_run(Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius,'tau-reff',Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.timeseries3(:,:,ilat,ilon,:) );
%        Y=histo_output.N_histo_mean;

         Y = Nd_timeseries.mean(ilat,ilon,:);

    case 'LWP std dev from histogram vals timeseries3'
        extra_title_info = [', LAT=' LAT_str ', LON=' LON_str];
        ylabelstr = 'Std dev LWP from histo timeseries3 (normalised)';

        Y = W_timeseries.std_dev_norm(ilat,ilon,:);
         
end

%%%%%%%%   now X and Y data are in place. Restrict the data if required


% set_MODIS_NH_flags=1; %gets reset every time
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
% MODIS_N_H_calc %runs the routine to calculate Nd
% N_H_calc_histo
% 
% thresh=0.8;
% 
% [MLAT2d,MLON2d]=meshgrid(MLON,MLAT); %these are the mid points e.g. MLAT=89.5 means 89-90 I think



%ihtot = find(abs(MLAT2d)<=30 & cf>=thresh);  thresh_str='LAT LTE 30 & CF GTE ';
%ihtot = find(abs(Solar_Zenith_Mean.data)>=70 & cf>=thresh);  thresh_str='SZA GTE 70 & CF GTE ';
%                        ihtot = find(cf>=thresh);  thresh_str='CF GTE';




%%%%%  choosing only certain points (lat-lon grid)
%                ihtot = find(totN>nthresh);  thresh_str='No. pixels'; %only plot for data with more a threshold no. pixels in total
%                ihtot = [1:prod(size(N_histo_mean))]; %all the data
%                ihtot = find(cf>=thresh);  thresh_str='CF'; %only plot for data with more a threshold value (of CF in this case)
%                 ihtot = find(WMOD>=thresh);  thresh_str='W times
%                 five-sixths';


%%%%%  choosing only certain points (timeseries)
%                ihtot = find(NP_time>50 & cf_time>0.8);
%                thresh_str='NP.GT.50 AND CF.GT.0.8';

thresh_NP=50;
thresh_CF=0.8;
thresh_SensZA=40;
thresh_tau = 100;
%thresh_tau = 20;
thresh_perror = 100;
thresh_perror_above = 98;
thresh_reff = 10;

switch data_type
    case 'timeseries'

        %thresh_lower_date = datenum('09-Oct-2004 19:00');
        %thresh_upper_date = datenum('13-Oct-2004 01:00');

        %times dhours either side of the 4 MPACE flights
        dhours=0;  %hours. datenum is in days so divide by 24 below
        dhours=1e99;  %hours. datenum is in days so divide by 24 below
        %set to 1e99 for all overpasses
        thresh_lower_date = datenum('09-Oct-2004 20:20') - dhours/24;
        thresh_upper_date = datenum('09-Oct-2004 22:07') + dhours/24;

        thresh_lower_date2 = datenum('10-Oct-2004 00:10')- dhours/24;
        thresh_upper_date2 = datenum('10-Oct-2004 03:10')+ dhours/24;

        thresh_lower_date3 = datenum('10-Oct-2004 21:30')- dhours/24;
        thresh_upper_date3 = datenum('10-Oct-2004 22:46')+ dhours/24;

        thresh_lower_date4 = datenum('12-Oct-2004 23:17')- dhours/24;
        thresh_upper_date4 = datenum('13-Oct-2004 00:00')+ dhours/24;

        %%%%%  choosing only certain points (timeseries3)- comment all out for all
        %%%%%  poiints
        %ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF)];
        %ihtot = find(NP_time3>thresh_NP);thresh_str=['NP.GT.' num2str(thresh_NP)];
        %ihtot = find(cf_time3>thresh_CF); thresh_str=['CF.GT.' num2str(thresh_CF)];
        %ihtot = find(NP_time3>thresh_NP & cf_time3>thresh_CF &
        %sensZA_time3<thresh_SensZA); thresh_str=['NP.GT.' num2str(thresh_NP) ' AND CF.GT.' num2str(thresh_CF) ' AND SensZA.GT.' num2str(thresh_SensZA)];
        %ihtot = find(scantime_matlab_tim>=thresh_lower_date & scantime_matlab_tim<=thresh_upper_date); thresh_str=['date.GTE.' datestr(thresh_lower_date,31) ' AND date.LTE.' datestr(thresh_upper_date,31)];

        ihtot_time = find( (scantime_matlab_tim_pdf>=thresh_lower_date & scantime_matlab_tim_pdf<=thresh_upper_date) | (scantime_matlab_tim_pdf>=thresh_lower_date2 & scantime_matlab_tim_pdf<=thresh_upper_date2) | (scantime_matlab_tim_pdf>=thresh_lower_date3 & scantime_matlab_tim_pdf<=thresh_upper_date3) | (scantime_matlab_tim_pdf>=thresh_lower_date4 & scantime_matlab_tim_pdf<=thresh_upper_date4) );
        if dhours>1e98
            thresh_str=['all overpasses'];
        else
            thresh_str=[num2str(dhours) ' either side of MPACE flights'];
        end
        ihtot=intersect(ihtot,ihtot_time);
        %intersect finds the values common to both arrays

        %%%%  cloud fraction threshold
        ihtot_cf = find(cf_tim(ilat,ilon,:)>thresh_CF); thresh_str=[thresh_str ' CF.GT.' num2str(thresh_CF)];
        ihtot=intersect(ihtot,ihtot_cf);
        %  intersect finds the values common to both arrays - so this only picks
        %  indices that were designated in the ihtot array above, but which also
        %  satisfy the CF condition

        %%%%  tau threshold
        %ihtot_tau = find(tau_5km_tim(ilat,ilon,:)<thresh_tau); thresh_str=[thresh_str ' tau.LT.' num2str(thresh_tau)];
        %ihtot=intersect(ihtot,ihtot_tau);

    case 'swath'
        %restrict by lat lon
%         thresh_str='';
%         
        ihtot_region = find(Plat>=min(LAT_val) & Plat<=max(LAT_val) & Plon>=min(LON_val) & Plon<=max(LON_val) );
        ihtot=intersect(ihtot,ihtot_region); thresh_str=''; extra_title_info = [', LAT=' LAT_str ', LON=' LON_str]; 
        
%        ihtot=intersect(ihtot,ipos_L2); thresh_str=''; extra_title_info = [', LAT=' LAT_str ', LON=' LON_str]; 
% 
% 
         ihtot_perror = find( percent_error_Nd(row_L2_inds,col_L2_inds) < thresh_perror );                
         ihtot = intersect(ihtot,ihtot_perror); thresh_str=[thresh_str ' % error.LT.' num2str(thresh_perror)];
%         
%         ihtot_perror_above = find( percent_error_Nd(row_L2_inds,col_L2_inds) > thresh_perror_above );                
% %        ihtot=intersect(ihtot,ihtot_perror_above); thresh_str=[thresh_str ' % error.GT.' num2str(thresh_perror_above)];
%         
         ihtot_liquid = find(phase_flag(row_L2_inds,col_L2_inds)==2); 
         ihtot=intersect(ihtot,ihtot_liquid); thresh_str=[thresh_str ' liquid  cloud'];
% 
%         ihtot_ocean = find(surface_flag(row_L2_inds,col_L2_inds)==0); 
% %        ihtot=intersect(ihtot,ihtot_ocean); thresh_str=[thresh_str ' OCEAN only'];
%         
        ihtot_tau = find(tau(row_L2_inds,col_L2_inds)<thresh_tau); 
        ihtot=intersect(ihtot,ihtot_tau); thresh_str=[thresh_str ' tau.LT.' num2str(thresh_tau)];
% 
%          ihtot_reff = find(re(row_L2_inds,col_L2_inds)>thresh_reff);
% %        ihtot=intersect(ihtot,ihtot_reff); thresh_str=[thresh_str ' Reff.GT.' num2str(thresh_reff)];
% 
%         ihtot_tau_bounds = find(tau_bounds(row_L2_inds,col_L2_inds)==2); 
% %        ihtot=intersect(ihtot,ihtot_tau_bounds); thresh_str=[thresh_str ' Tau within bounds '];
%         
end

%remove some of the points if not required
X=X(ihtot);
Y=Y(ihtot);

%plot_type
short_plot_name=[x_axis_vals ' vs ' y_axis_vals ];
tit(1).tit=[short_plot_name extra_title_info];
savename=[savedir tit(1).tit];








if minXbins<-8.9e99
    Xbins=make_PDF_bins(X,nXpdf); %if not set limits then use the default of minALL(X):dX:maxALL(X)
else
    Xbins=make_PDF_bins(X,nXpdf,minXbins,maxXbins);
end

if minYbins<-8.9e99
    Ybins=make_PDF_bins(Y,nYpdf);
else
    Ybins=make_PDF_bins(Y,nYpdf,minYbins,maxYbins);
end





