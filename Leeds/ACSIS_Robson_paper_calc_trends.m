yr_start_trend=1850;        
%yr_start_trend=2003;
yr_end_trend=2014;

iscreen_sig=1; %Whether to screen for signficance
marker_size=1; %For non-signficant points
%marker_size=1; %For non-signficant points
iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid

p_conf = 95; % Confidence limit (%) for the trend significance
nthresh_days = 3;
%nthresh_days = 0;




UKESM_Nd_case = 'to ztop';
UKESM_Nd_case = 'to 3.2km';



icoarse_grain=0;

time_round='';
time_format_str='';
icontour_DRIVER=0;
isave_plot=0;
iplot_wind_arrows=0;

cont_col_str_DRIVER='k';
  

%% Trend analysis (linear least squares fit) - model
    yr_start=yr_start_trend; yr_end=yr_end_trend;
    yr_start_trend_used=yr_start_trend; yr_end_trend_used=yr_end_trend;    
    istart=find(dat_ukesm.years_ukesm_1d==yr_start);
    iend=find(dat_ukesm.years_ukesm_1d==yr_end);    
    x = [yr_start:yr_end]';
%     y = Nd_annual(istart:iend,100,100);
%     [Nd_trend,bint,residuals,rint,stats] = regress(y,x);
%     
%     [coeffs,bint,residuals,rint,stats] = trend_linear_fit_nd_data(Nd_annual(istart:iend,:,:),1);
    
   [coeffs,t_trend] = trend_linear_fit_nd_data(x,dat_ukesm.dat_annual(istart:iend,:,:),1); 
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    %t_thresh = tinv(p_conf/100,n_dof); %find the t value needed for 95% confidence using a one-tailed t-test
        %N.B. - here is the function for a 2-tailed test :-
            % E.g. 
            % t=4; v=10;
            % tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));   % 2-tailed t-distribution function
            % tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2; %1-tailed t-distribution function
            % OR just :-  
            % tail2P = 2*tcdf(-abs(t),v);
            % tail1P = tcdf(-abs(t),v);
            % Test :-            
            % T2 = [1-tdist2T(t,v)  tail2P]; %give the same answer
            % T1 = [1-tdist1T(t,v)  tail1P]; %give the same answer
    %itrend_not_sig=find(abs(t_trend)<t_thresh); 
    
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend),n_dof);
    T1 = T2/2;    
    itrend_not_sig = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig = find((1-T1)<=p_conf/100); %1-tailed
   


%% Obs data - calc annual means, etc.

switch var_ukesm
    
    case 'calipso_low_cloud_amount';
        obs_str='CALIPSO';
        
        
        %Load Calipso data - currently 2007-2017
        % Load CF data using :-
        %read_calipso_monthly_IPSL_2007_2017
        years_obs = years_requested;                
        
        
        
        
       
        
    case 'Nd_cf_weighted_UKESM';
        obs_str='MODIS';
        
        
        % MODIS data
        % Using the data given to Jane - screened for sea-ice etc.
        cf_screen_str = 'CF>80';
        cf_screen_str = 'CF>0';
        
        str_2137='21';
        str_2137='37';
        
        res_str='1deg';
        res_str='1km';
        
        
        switch cf_screen_str
            case 'CF>80'
                file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/';
                dataset_str = 'SZA_LT_65_CF_GT_80_CTH_LT_3.2km_screened_for_seaice__2week_max';
                dataset_str = 'SZA_LT_65_CF_GT_80_CTH_LT_3.2km';
                
            case 'CF>0'
                file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.0_meanCTT_173_meanCTH_3.2km_SZA_65/';
                dataset_str = 'SZA_LT_65_CF_GT_0_CTH_LT_3.2km';
        end
        
        
        
        
        switch str_2137
            case '21'
                str_label_2137='2.1 um';
            case '37'
                str_label_2137='3.7 um';
        end
        
        years_MODIS2=[2003:2014];        
        clear mon_me_MODIS2 mon_me_MODIS2_Ndatap
        for iy=1:length(years_MODIS2)
            year_str = num2str(years_MODIS2(iy));
            filename = [file_dir 'Nd_monthly_' str_2137 '_' res_str '_' year_str '_' dataset_str '.mat.nc'];
            nc=netcdf(filename);
            if iy==1
                lat=nc{'lat'}(:);
                lon=nc{'lon'}(:);
                [gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE]=meshgrid(lon,lat);
            end
            
            mon_me_MODIS2{iy} = nc{['Nd_' res_str '_mean']}(:);
            mon_me_MODIS2_Ndatap{iy} = nc{['Nd_' res_str '_Ndatap']}(:);
            inan = find(mon_me_MODIS2_Ndatap{iy} < nthresh_days);
            mon_me_MODIS2{iy}(inan)=NaN;
            %mon_me_filter{iy}(inan)=NaN;
        end
        
     
        
        clear Nd_annual_MODIS
        for iy=1:size(Nd_MODIS,1)            
            Nd_annual_MODIS(iy,:,:) = meanNoNan(mon_me_MODIS2{iy},3);
        end
        
        years_obs = years_MODIS2;
        
end

%% Trend analysis (linear least squares fit) - MODIS
    yr_start_MODIS=2003; yr_end_MODIS=2014;
    istart_MODIS=find(years_MODIS2==yr_start_MODIS);
    iend_MODIS=find(years_MODIS2==yr_end_MODIS);                
    x = [yr_start_MODIS:yr_end_MODIS]';
    %The following takes a litte while...
   [coeffs_MODIS,t_trend_MODIS] = trend_linear_fit_nd_data(x,Nd_annual_MODIS(istart_MODIS:iend_MODIS,:,:),1); 
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    %t_thresh = tinv(p_conf/100,n_dof); %find the t value needed for 95% confidence using a one-tailed t-test
        %N.B. - here is the function for a 2-tailed test :-
            % E.g. 
            % t=4; v=10;
            % tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));   % 2-tailed t-distribution function
            % tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2; %1-tailed t-distribution function
            % OR just :-  
            % tail2P = 2*tcdf(-abs(t),v);
            % tail1P = tcdf(-abs(t),v);
            % Test :-            
            % T2 = [1-tdist2T(t,v)  tail2P]; %give the same answer
            % T1 = [1-tdist1T(t,v)  tail1P]; %give the same answer
    %itrend_not_sig_MODIS = find(abs(t_trend_MODIS)<=t_thresh); %make NaN for now, but can put a dot on, etc.
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend_MODIS),n_dof);
    T1 = T2/2;    
    itrend_not_sig_MODIS = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig_MODIS = find((1-T1)<=p_conf/100); %1-tailed
    
    