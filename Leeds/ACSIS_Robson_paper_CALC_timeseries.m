function [dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std,...
    trend_dat_box,trend_dat_box_ens,trend_ens,dat_annual_box_ukesm_ens,me_t_PI,N_t_PI,std_t_PI]=ACSIS_Robson_paper_CALC_timeseries(dat_ukesm,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
    yr_start_trend_box,yr_end_trend_box,p_conf,ibox)

%% calculate some timeseries in the regional box for the model
if iscell(LON_val) %For situation where specify more than one region for long values in order to cross
    %the date line for example.
    ilat = find(dat_ukesm.gcm_Plat2D_UM(:,1)>LAT_val(1) & dat_ukesm.gcm_Plat2D_UM(:,1)<LAT_val(2));
    ilon=[];
    for iregions=1:length(LON_val)        
        ilon_reg = find(dat_ukesm.gcm_Plon2D_UM(1,:)>LON_val{iregions}(1) & dat_ukesm.gcm_Plon2D_UM(1,:)<LON_val{iregions}(2));
        ilon = [ilon ilon_reg];        
    end
else
    ilat = find(dat_ukesm.gcm_Plat2D_UM(:,1)>LAT_val(1) & dat_ukesm.gcm_Plat2D_UM(:,1)<LAT_val(2));
    ilon = find(dat_ukesm.gcm_Plon2D_UM(1,:)>LON_val(1) & dat_ukesm.gcm_Plon2D_UM(1,:)<LON_val(2));
end
%     clear dat_time_mean
%     for iy=1:size(Nd_ukesm,1)
%         dat_time_mean(iy,:) = meanNoNan(meanNoNan(dat.dat_ukesm(iy,:,ilat,ilon),4),2);        
%     end

    if iscreen_land==1
        load_type = 'merged netCDF';
        var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
        um_case='u-bf666'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
              
        
        land_mask_mean = squeeze(ones(size(dat_ukesm.dat_annual(1,:,:))));
        land_mask_mean(dat_global.dat==1)=NaN;
                
        land_mask_ens = repmat(land_mask_mean,[1 1 size(dat_ukesm.dat_annual_ens,1)]);
        land_mask_ens = permute(land_mask_ens,[3 1 2]);
        
    else
        land_mask_mean = squeeze(ones(size(dat_ukesm.dat_annual(1,:,:))));
        land_mask_ens = repmat(land_mask_mean,[1 1 size(dat_ukesm.dat_annual_ens,1)]);
        land_mask_ens = permute(land_mask_ens,[3 1 2]);
                     
    end
    
  
    clear dat_annual_box_ukesm dat_annual_box_ukesm_ens dat_annual_box_ukesm_ens dat_annual_box_ukesm_ens_std
    for it=1:size(dat_ukesm.dat_annual,1)
        dat_mean = squeeze(dat_ukesm.dat_annual(it,:,:)) .* land_mask_mean;        
        dat_mean_ens = squeeze(dat_ukesm.dat_annual_ens(:,it,:,:)) .* land_mask_ens;
        
        dat_tmp = dat_mean(ilat,ilon);
        area_tmp = gcm_area_UM(ilat,ilon);                
        
        %weight by area of gridbox
        dat_annual_box_ukesm(it) = meanNoNan(dat_tmp(:),1,'',0,1,area_tmp(:)); 
         % function [me,nnums,stdev]=meanNoNan(dat,ndim,op,isqueeze,iweight,weights)
        
        
        area_rep = repmat(gcm_area_UM(ilat,ilon),[1 1 size(dat_mean_ens,1)]);
        area_rep = permute(area_rep,[3 1 2]);
        dat_tmp2 = dat_mean_ens(:,ilat,ilon);  
        
        %weight by area of gridbox
        dat_tmp = meanNoNan(dat_tmp2(:,:),2,'',0,1,area_rep(:,:));
        dat_annual_box_ukesm_ens_std(it) = std(dat_tmp); %std dev across the ensemble
        dat_annual_box_ukesm_ens(it,:) = dat_tmp; %Store the area averages for all 9 members of the ensemble.
    end
    
    land_mask_mean_all_t = repmat(land_mask_mean,[1 1 size(dat_PI.dat_annual,1)]);
    land_mask_mean_all_t = permute(land_mask_mean_all_t,[3 1 2]);
    
    
    dat_tmp = dat_PI.dat_annual(:,ilat,ilon) .* land_mask_mean_all_t(:,ilat,ilon);
    area_tmp = gcm_area_UM(ilat,ilon);
    area_rep = repmat(area_tmp,[1 1 size(dat_tmp,1)]);
    area_rep = permute(area_rep,[3 1 2]);
    dat_PI_box = meanNoNan(dat_tmp(:,:),2,'',0,1,area_rep(:,:)); 
    
    
    [me_t_PI,N_t_PI,std_t_PI] = meanNoNan(dat_PI_box,1);
    
    
    
    
    
%% Calculate the MODEL trend for the boxed region
for it_trend=1:length(yr_start_trend_box)
%     yr_start=yr_start_trend_box(it_trend); yr_end=yr_end_trend_box(it_trend);
%     yr_start_trend_used_box=yr_start; yr_end_trend_used_box=yr_end;    
%     istart=find(dat_ukesm.years_ukesm_1d==yr_start);
%     iend=find(dat_ukesm.years_ukesm_1d==yr_end);    
%     x = [yr_start:yr_end]';
%     y = dat_annual_box_ukesm(istart:iend);
% %     y = Nd_annual(istart:iend,100,100);
% %     [Nd_trend,bint,residuals,rint,stats] = regress(y,x);
% %     
% %     [coeffs,bint,residuals,rint,stats] = trend_linear_fit_nd_data(Nd_annual(istart:iend,:,:),1);
%     
%    [coeffs,t_trend,bint,stats,p_Durbin,d_Durbin] = trend_linear_fit_nd_data(x,y,2); 
%    
%    ylin=coeffs(1)+coeffs(2).*x; %The straight line for the linear trend
%    
%    %t-test threshold for 95% confidence
%     n_dof = length(x)-2; %number of degrees of freedom
%     %t_thresh = tinv(p_conf/100,n_dof); %find the t value needed for 95% confidence using a one-tailed t-test
%         %N.B. - here is the function for a 2-tailed test :-
%             % E.g. 
%             % t=4; v=10;
%             % tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));   % 2-tailed t-distribution function
%             % tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2; %1-tailed t-distribution function
%             % OR just :-  
%             % tail2P = 2*tcdf(-abs(t),v);
%             % tail1P = tcdf(-abs(t),v);
%             % Test :-            
%             % T2 = [1-tdist2T(t,v)  tail2P]; %give the same answer
%             % T1 = [1-tdist1T(t,v)  tail1P]; %give the same answer
%     %itrend_not_sig=find(abs(t_trend)<t_thresh); 
%     
%     % Significance of our t values :-
%     T2 = 2*tcdf(-abs(t_trend),n_dof);
%     T1 = T2/2;    
%     itrend_not_sig = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
%     %itrend_not_sig = find((1-T1)<=p_conf/100); %1-tailed
% 
%     trend_dat_box{ibox,it_trend}.coeffs = coeffs;
%     trend_dat_box{ibox,it_trend}.bint = bint;
%     trend_dat_box{ibox,it_trend}.stats = stats;
%     trend_dat_box{ibox,it_trend}.p_Durbin = p_Durbin;
%     trend_dat_box{ibox,it_trend}.d_Durbin = d_Durbin;        
%     trend_dat_box{ibox,it_trend}.t_trend = t_trend;
%     trend_dat_box{ibox,it_trend}.ylin = ylin;
%     trend_dat_box{ibox,it_trend}.x = x;
%     trend_dat_box{ibox,it_trend}.y = y;
%     trend_dat_box{ibox,it_trend}.itrend_not_sig = itrend_not_sig;
%     trend_dat_box{ibox,it_trend}.T1 = T1;
%     trend_dat_box{ibox,it_trend}.T2 = T2;
    
    [trend_dat_box{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, ...
            it_trend, dat_ukesm.years_ukesm_1d, dat_annual_box_ukesm, p_conf);
    
end
    

%% Calc trends for individual ensemble members
for it_trend=1:length(yr_start_trend_box)
    clear trend_ens
    for iens=1:size(dat_annual_box_ukesm_ens,2)
        [trend_dat_box_ens{ibox,it_trend,iens}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, ...
            it_trend, dat_ukesm.years_ukesm_1d, dat_annual_box_ukesm_ens(:,iens)', p_conf);
        
        trend_ens(iens) = trend_dat_box_ens{ibox,it_trend,iens}.coeffs(2);
    end
    %min, max, mean, std
    
    [me,N,st] = meanNoNan(trend_ens,2);
    
    trend_dat_box_stats{ibox,it_trend}.mean=me;
    trend_dat_box_stats{ibox,it_trend}.std=st;
    trend_dat_box_stats{ibox,it_trend}.min=min(trend_ens);
    trend_dat_box_stats{ibox,it_trend}.max=max(trend_ens);
    
end
