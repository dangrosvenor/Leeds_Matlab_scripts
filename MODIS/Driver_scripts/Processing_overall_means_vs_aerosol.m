
% Plot domain and time mean LWP and RWP vs Nd
% Data is loaded from a .mat files for the timeseries, which are generated
% by UM_process_runs_for_timeseries and stored in each directory

isave_plot_overall=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

% --- Add run cases below to UM_case_select_runs
UM_cases = '12th Nov case ACI, as of May 2016';
UM_cases = '12th Nov case ACI, as of May 2016 with sed runs'; %N.B - messes up the bar chart for SW etc.
            %Also changes idat_micro to only plot line for the
            %sedimentation runs
UM_cases = '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM';

UM_case_select_RUN; %external scrip to select the case and put all the variables into current memory space


% save name for the figure generated here :-
savename_ACI = '/home/disk/eos1/d.grosvenor/modis_work/plots/UM/VOCALS_Processing_ACI';

%LWP_overall_file = '/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/LWP_overall_20160707T103638.mat';
%SW_overall_file = '/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/SW_values_for_ACI_20160809T092540.mat';
%Also add one that uses the daytime only.
%SW=load(LWP_overall_file);

%back these up
%marker_style_load = marker_style;
%line_colour_load = line_colour;
%line_pattern_load = line_pattern;


nsub=0; %counter for which plot we are on


% vars to load timeseries for
clear var_list xdat_save ydat_save
i=1;
var_list{i} = 'SW_up_TOA'; i=i+1;
var_list{i} = 'LWP'; i=i+1;
var_list{i} = 'LWP_incloud_20gmsq'; i=i+1;
var_list{i} = 'RWP'; i=i+1;
var_list{i} = 'Nd'; i=i+1;
var_list{i} = 'CF_LWP_20'; i=i+1;
%var_list{i} = 'CF_0pt25_LWP_4km_20'; i=i+1;
%var_list{i} = 'UM_time_in'; i=i+1;

%Tests using these different weightings show that using the different
%weights makes very little difference ('choose' vs 'own' option). Using no
%weights affects the magnitude a little, but not the relationships.
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC) for VOCALS (add this on for LST)
clear time_range_weight
iweight=1; %need to set a default value
weight='none';
%weight='own'; %use the SW_up_TOA weighting from each individual run - probably should not use this
%weight='choose'; iweight=4; %Use the same weighting for all runs (from the iweight run in the list)
%weight='SW_down_TOA'; %Use the SW_down_TOA from xmmz-x for all runs
%weight='Time range'; time_range_weight{1}=[datenum('13-Nov-2008 11:45') datenum('13-Nov-2008 12:15')]-time_shift;    %Use one time range only
%weight='Time range'; time_range_weight{1}=[datenum('13-Nov-2008 10:00') datenum('13-Nov-2008 14:00')]-time_shift;    %Use one time range only
%weight='Time range'; time_range_weight{1}=[datenum('12-Nov-2008 10:00') datenum('12-Nov-2008 14:00')]-time_shift;    %Use one time range only
%    time_range_weight{2}=[datenum('13-Nov-2008 10:00') datenum('13-Nov-2008 14:00')]-time_shift;    %Use one time range only


switch UM_cases
    case '12th Nov case ACI, as of May 2016 with sed runs'
        idat_micro = [2 4 5 6]; %for joining the lines together
      case '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM'
         clear idat_micro
        idat_micro{1} = [1 3 5 7]; %for joining the lines together   
        idat_micro{2} = [2 4 6 8]; %for joining the lines together
    otherwise
        idat_micro = [1:length(labs_UM)]; %for joining the lines together
end



%Load in SW_down as a special case, since only processed it for xmmz-x
SW_down_TOA = load('/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/xmmz-x/xmmzx_SW_down_TOA_.pp.nc_timeseries.mat');
w = SW_down_TOA.timeseries_UM(1:end-1);
SW_down_mean  = sum( SW_down_TOA.timeseries_UM(1:end-1) .* w ) ./ sum(w);
% This value is very high (1090 W/m2) since it is TOA - really need to know
% what is coming in to the cloud level - SW at surface will do I think
% (stash 1-235)
% Am moving this over now

Transmission_SW_surf = load('/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/xmmz-x/xmmzx_Transmission_down_surf_LWP_LT_0pt1_.pp.nc_timeseries.mat');
SW_down_surf = load('/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/xmmz-x/xmmzx_SW_down_surf_LWP_LT_0pt1_.pp.nc_timeseries.mat');



%subplotting settings
subplotting_DRIVER=1;
xsub=length(var_list)-1;
ysub=1;

%Puts all the timeseries into e.g. SW_up_TOA(irun).timeseries
for irun=1:length(fileUM) 
    
    if iscell(dirUM)==1
            dirUM_i = dirUM{irun};
        else
            dirUM_i = dirUM;
        end
    
    
    for ivar=1:length(var_list)
        var_str = var_list{ivar};
        if irun==1
            eval(['clear ' var_str]);
        end
        filename_tim = [remove_character( [dirUM_i fileUM{irun}], 'VAR_NAME',var_str) '_timeseries.mat'];
        %Load the data
        eval_str=[var_str '(irun)=load(filename_tim);'];
        eval(eval_str);        
    end    
end

Processing_overall_LWP_ACI_v2_SW_subplot_v2
Processing_overall_LWP_ACI_v2_LWP_subplot_v2
Processing_overall_LWP_ACI_v2_LWPic_subplot_v2
Processing_overall_LWP_ACI_v2_RWP_subplot_v2
Processing_overall_LWP_ACI_v2_CF_subplot_v2

%Based on xmmz-x only (low cloud fraction run)
SW_down_mean_surf  = sum( SW_down_surf.timeseries_UM .* w ) ./ sum(w); 


savename = [savename_ACI];



%% SW calculations and partitioning

savename_SW = [savename_ACI '_SW_partitioning'];

% Do some calculations with the numbers RE breaking down the SW effect.
% Need a formula for albedo as func of LWP_ic, CF and Nd.
% See Sc_calc_albedo_vs_W_Nd
% tau = (N* W ^(5/2) / B )^(1/3)
% albedo = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)

% Data from the N subplots and M runs is stored in the array ydat_save(M,N)

iSW=1; %SW_out
iW=3; %in-cloud LWP
icf=5; %cf
%Nd is in xdat_save

Nd = xdat_save(:,1);
%sort into ascending order
[Nd,isort] = sort(Nd,1,'ascend');
cf = ydat_save(isort,icf);
W = ydat_save(isort,iW);
SW = ydat_save(isort,iSW); %SW_TOA_out
labs_bar = labs(isort);

dSW = SW(2:end) - SW(1:end-1);


%Need to get mean SW_in from the data
%SW_in = 321.4; %
%SW_in = 852.9; %model SW_Down at surface suggests this for clear-sky, but
%gives computed SW_up_TOA differences that are too large...
%SW_in = 1118.3; %For the single time of 17:00 UTC (noon) 13th Nov

%value calculated from very low aerosol run (xmmz-x)
SW_in = SW_down_mean_surf;

%Clear sky transmisison of the atmosphere (SW_down_surf / SW_down_TOA)
transmission_atmos = 1; %model values suggest around 0.8
transmission_atmos = 0.75; %model values suggest around 0.8


%Estimate surface albedo - can't do this yet since need the cloud free SW_TOA_up
%A_surf_estimate = SW(1)./transmission_atmos ./ SW_down_mean_surf;

% Also need the clear-sky albedo - can estimate from cloud free regions?
A_clear = 0.15; %Guess for now - value for CF quite sensitive to this (ranging up to 30%), so should prob
  %check what this is from the model
A_clear = 0.12; %value estimated from the 12 noon snapshot - N.B. - this uses the estimated transmission since
 %it is the surface albedo assuming no atmosphere above

clear dSW_f dSW_W dSW_N Ac tau A_f A_w A_N SW_f SW_W SW_N bar_lab SW_calc

SW_partition_method = 'linear_extrapolation of gradient';
SW_partition_method = 'one-at-a-time sensitivity';

switch SW_partition_method
    case 'one-at-a-time sensitivity'
         % Do calcs for each run and the next highest aerosol run
        for i=1:length(Nd)
            i0 = i; 
            i1=i0+1;
            i1=min(i1,length(Nd));

            f0 = cf(i0); f1 = cf(i1);
            W0 = W(i0); W1 = W(i1);
            N0 = Nd(i0); N1 = Nd(i1);
            % Needed for the actual SW change for comparison :-
            SW0 = SW(i0); SW1 = SW(i1);
            
%            [Ac(i),tau(i)] = albedo_cloudy_func_Seinfeld(W,Nm3)

%            [temp,temp,temp,Ac_calc(i),tau_calc(i),A_calc(i),SW_calc(i)] = SW_partitioning(f0,W0,N0,1,1,1,SW_in,A_clear,transmission_atmos);
            [Ac_calc(i),tau_calc(i),A_calc(i),SW_calc(i)] = calc_SW(f0,W0,N0,SW_in,A_clear,transmission_atmos);            
            if i<=length(Nd)-1
%                [temp,temp,temp,Ac(i),tau(i),A_f(i),SW_f(i)] = SW_partitioning(f1,W0,N0,1,1,1,SW_in,A_clear,transmission_atmos);
%                [temp,temp,temp,Ac(i),tau(i),A_W(i),SW_W(i)] = SW_partitioning(f0,W1,N0,1,1,1,SW_in,A_clear,transmission_atmos);
%                [temp,temp,temp,Ac(i),tau(i),A_N(i),SW_N(i)] = SW_partitioning(f0,W0,N1,1,1,1,SW_in,A_clear,transmission_atmos);
                
                [Ac(i),tau(i),A_f(i),SW_f(i)] = calc_SW(f1,W0,N0,SW_in,A_clear,transmission_atmos);
                [Ac(i),tau(i),A_W(i),SW_W(i)] = calc_SW(f0,W1,N0,SW_in,A_clear,transmission_atmos);
                [Ac(i),tau(i),A_N(i),SW_N(i)] = calc_SW(f0,W0,N1,SW_in,A_clear,transmission_atmos);
                
%                bar_lab{i} = {[labs_bar(i1).l ' vs ']; labs_bar(i0).l};
                bar_lab{i} = [labs_bar(i1).l ' vs\newline' labs_bar(i0).l]; 
            end
            


        end
        

        
        
        clear bar_dat
        bar_dat([1 3 5],1) = SW_f - SW_calc(1:length(Nd)-1);
        bar_dat([1 3 5],2) = SW_W - SW_calc(1:length(Nd)-1);
        bar_dat([1 3 5],3) = SW_N - SW_calc(1:length(Nd)-1);
        bar_dat([1 3 5],4) = 0; %dummy bar for a different colour
        
        %Add the bars for the actual SW change
        bar_dat(2,4) = dSW(1);
        bar_dat(2,1:3) = 0;
        bar_dat(4,4) = dSW(2);
        bar_dat(4,1:3) = 0;
        bar_dat(6,4) = dSW(3);
        bar_dat(6,1:3) = 0;

        figure
        bar(bar_dat,'stacked'); %will give N bars each with M components for bar_dat[N,M]
        title(['SW_{in}=' num2str(SW_in) ', A_{clear}=' num2str(A_clear)]);
        ylabel('\DeltaSW_{up TOA} (W m^{-2})');
        
        %Sort out bar labelling
        pos = get(gca,'xtick'); %get locations
        mid = (pos(2:end)+pos(1:end-1) )/2; %Consolidate into half as many at mid-points
        set(gca,'xtick',mid(1:2:end));
%         clear bar_lab_row01
%         for i=1:length(bar_lab)
%             bar_lab_row01{i} = bar_lab{i}{1};
%             bar_lab_row02{i} = bar_lab{i}{2};
%         end
%        set(gca,'xticklabel',bar_lab_row01);
        [hx,hy] = format_ticks(gca,bar_lab);
        set(hx,'Fontsize',18);
        
        
        %add a legend
        legend({'CF','LWP_{ic}','N_d','Actual'});
        
        %Add labels for the percentage for each bar
        fsize_percents = 14;
        tots = sum(abs(bar_dat),2); %Do percentages of the sum of the absolute changes since
          %some changes can be negative
        tots = repmat(tots,[1 size(bar_dat,2)]);
        percents = 100 * bar_dat ./ tots;
        
        %Get positions to put the labels
        yvals = cumsum(bar_dat,2); %the height of the top of each bar
        yvals2 = zeros([size(bar_dat,1) size(bar_dat,2)+1]);
        yvals2(:,2:end) = yvals;
        yvals_mid = (yvals2(:,2:end) + yvals2(:,1:end-1) ) /2;
        
        cols={'c','c','','k','k','k','k','k','k'}; %The order of this runs across the bars
          %and then up to the next stack, etc.
        iskip=[3]; %[3 3 3]; %i and j indices of bars to avoid labelling for (if very small)
        jskip=[1]; %[1 3 3];
        NI=3;
        NJ=3;
        for i=1:NI %loop over the vertical bar segments
            for j=1:NJ %loop over the bars

                linind = sub2ind([NI NJ],i,j);
                linind2 = sub2ind([NI NJ],iskip,jskip);
                ifind=find(linind2==linind);

                if length(ifind)==0
                    ii = (i-1)*2+1;
                    x = ii-0.25;
                    y = yvals_mid(ii,j);
                    numstr = num2str(percents(ii,j),'%2.1f');
                    text(x,y,[numstr ' %'],'color',cols{linind},'fontsize',fsize_percents);
                end
            end
        end
        
        fontsize_figure(gcf,gca,18);  
        
        
        
        fprintf(1,'\n SW_calc ./ SW = %f\n',SW_calc ./ SW');  %ratios from the other script
                
        
    case 'linear_extrapolation of gradient'

        % Do calcs for each run and the next highest aerosol run
        for i=1:length(Nd)-1
            i0 = i; i1=i+1;
            f0 = cf(i0); f1 = cf(i1);
            W0 = W(i0); W1 = W(i1);
            N0 = Nd(i0); N1 = Nd(i1);
            % Needed for the actual SW chagne for comparison :-
            SW0 = SW(i0); SW1 = SW(i1);

            [dSW_f(i),dSW_W(i),dSW_N(i),Ac(i),tau(i)] = SW_partitioning(f0,W0,N0,f1,W1,N1,SW_in,A_clear);
        end

        clear bar_dat
        bar_dat(:,1) = dSW_f./dSW'
        bar_dat(:,2) = dSW_W./dSW'
        bar_dat(:,3) = dSW_N./dSW'

        figure
        bar(bar_dat,'stacked'); %will give N bars each with M components for bar_dat[N,M]
        title(['SW_{in}=' num2str(SW_in) ', A_{clear}=' num2str(A_clear)]);

        dSW_rel = SW(2:end)./SW(1:end-1) - 1
        dcf_rel = cf(2:end)./cf(1:end-1) - 1
        dW_rel = W(2:end)./W(1:end-1) - 1
        dN_rel = Nd(2:end)./Nd(1:end-1) - 1

        clear bar_dat2
        bar_dat2(:,1) = [dSW_f(1); dSW(1); dSW_f(2); dSW(2); dSW_f(3); dSW(3)];
        bar_dat2(:,2) = [dSW_W(1); 0; dSW_W(2); 0; dSW_W(1); 0];
        bar_dat2(:,3) = [dSW_N(1); 0; dSW_N(2); 0; dSW_N(3); 0];

        figure
        bar(bar_dat2,'stacked');
        title(['SW_{in}=' num2str(SW_in) ', A_{clear}=' num2str(A_clear)]);

end


%% end
