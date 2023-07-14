function trend_dat_box_obs2 = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box2,years_obs2,yr_start,yr_end,p_conf,it_trend,ibox)


    istart=find(years_obs2==yr_start);
    iend=find(years_obs2==yr_end);    
    x = [yr_start:yr_end]';
    y = obs_annual_box2(istart:iend);

    %test whether the significance is affected by subtracting the mean -
    % -- makes no difference --
    %me_lin = meanNoNan(obs_annual_box(istart:iend),2);
    %y = obs_annual_box(istart:iend) - me_lin;
    
if length(y)>2
    
   [coeffs,t_trend] = trend_linear_fit_nd_data(x,y,2); 
   
   ylin=coeffs(1)+coeffs(2).*x; %The straight line for the linear trend
   
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

    trend_dat_box_obs2{ibox,it_trend}.coeffs = coeffs;
    trend_dat_box_obs2{ibox,it_trend}.t_trend = t_trend;
    trend_dat_box_obs2{ibox,it_trend}.ylin = ylin;
    trend_dat_box_obs2{ibox,it_trend}.x = x;
    trend_dat_box_obs2{ibox,it_trend}.itrend_not_sig = itrend_not_sig;
    trend_dat_box_obs2{ibox,it_trend}.T1 = T1;
    trend_dat_box_obs2{ibox,it_trend}.T2 = T2;
    
end
