function [trend_dat_box,istart,iend] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, it_trend, years_ukesm_1d, dat_annual_box_ukesm, ...
    p_conf)

sdat=size(dat_annual_box_ukesm);

%for it_trend=1:length(yr_start_trend_box)
    yr_start=yr_start_trend_box(it_trend); yr_end=yr_end_trend_box(it_trend);
    yr_start_trend_used_box=yr_start; yr_end_trend_used_box=yr_end;    
    istart=find(years_ukesm_1d==yr_start);
    iend=find(years_ukesm_1d==yr_end);    
    x = [yr_start:yr_end]';
    if length(sdat)==3
        y = dat_annual_box_ukesm(istart:iend,:,:);
        dim=1;
    else
        y = dat_annual_box_ukesm(istart:iend);
        dim=2;
    end
    
%    if length(y)==0    
        
%    else
        [coeffs,t_trend,bint,stats,p_Durbin,d_Durbin,t_trend2,uncer2,uncer,Neff] = trend_linear_fit_nd_data(x,y,dim);
%    end
    
   ylin=coeffs(1)+coeffs(2).*x; %The straight line for the linear trend
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    n_dof2 = Neff-2; %dof for adjusted sample number (based on lag-1 auto-corr)
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
    
    % Significance of our t values using effective sample size based on lag-1 autocorrelation :-
    T2_2 = 2*tcdf(-abs(t_trend2),n_dof2); %works ok with n_dof2 as a matrix (or singular)
    T1_2 = T2_2/2;    
    itrend_not_sig2 = find((1-T2_2)<=p_conf/100);

    trend_dat_box.coeffs = coeffs;
    trend_dat_box.t_trend = t_trend;
    trend_dat_box.bint = bint;
    trend_dat_box.stats = stats;
    trend_dat_box.p_Durbin = p_Durbin;
    trend_dat_box.d_Durbin = d_Durbin;   
    trend_dat_box.ylin = ylin;
    trend_dat_box.x = x;
    trend_dat_box.y = y;
    trend_dat_box.itrend_not_sig = itrend_not_sig;
    trend_dat_box.T1 = T1;
    trend_dat_box.T2 = T2;
    trend_dat_box.uncer = uncer; %2 times the standard error for the fit
    %Auto-con lag-1 adjusted values :-
    trend_dat_box.itrend_not_sig2 = itrend_not_sig2;
    trend_dat_box.T1_2 = T1_2;
    trend_dat_box.T2_2 = T2_2;
    trend_dat_box.T1_max = max(T1,T1_2);
    trend_dat_box.T2_max = max(T2,T2_2);
    trend_dat_box.uncer2 = uncer2;
    trend_dat_box.uncer_max = max(uncer,uncer2);
    
    if length(istart)==0
        istart=NaN;
    end
    if length(iend)==0
        iend=NaN;
    end
    
%end

