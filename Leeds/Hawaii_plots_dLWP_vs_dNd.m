clear var_test


%% Plot dLWP for dNd bins with LWP, SO2, longitude and time filters too
% The idea was to try and restrict the region to the volcanic plume rather
% than including regions where there was no influence at all.

%plot_type = 'dLWP vs dNd';
%plot_type = 'dNd vs dNaccum';
%plot_type = 'dCFproxy_LWP vs dNd';
%plot_type = 'dLWP vs SO2';
%plot_type = 'LWP_PI vs SO2';
%plot_type = 'LWP_PD vs SO2';
%plot_type = 'Nd_PI vs SO2';
%plot_type = 'Nd_PD vs SO2';
%plot_type = 'dNd vs SO2';
%plot_type = 'LWP_PD vs Nd_PD';
%plot_type = 'LWP_PI vs Nd_PD';
%plot_type = 'LWP_PI vs Nd_PI';
%plot_type = 'dLWP vs Nd_PD';


SO2_filter_type = 'delta';
SO2_filter_type = 'Volc ON';

%SO2 thresholds
thresh_SO2 = -1;
thresh_SO2 = 1e-5;
%thresh_SO2 = 1e-4;

%Longitude thresholds
%LON_filter = [-160 -155];
%LON_filter = [-165 -160];   
%LON_filter = [-170 -165];
%LON_filter = [-175 -170];
%LON_filter = [-361 361];



%LWP thresholds
%thresh_LWP = -1;
%thresh_LWP = 50;

%Default
prc_vals_for_bins = [0:1:100];

%prc_vals_for_bins = [0:2.5:100];


Nd_type_str='';
switch plot_type
    case 'dLWP vs dNd'
        var_test = LWP_PD_ALL - LWP_PI_ALL;
        
        Nd_type = 'CASIM';
        %Nd_type = 'UKCA';
        %Nd_type = 'vs Rain Rate 0.05';
        %Nd_type = 'vs Rain Rate 0.01';
        
        xlab_str = ['\DeltaN_{d ' Nd_type '} (cm^{-3})'];
        
        switch Nd_type
            case 'CASIM'
                var_filter = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
            case 'UKCA'
                var_filter = (Nd_UKCA_PD_ALL - Nd_UKCA_PI_ALL)/1e6;
                Nd_type_str = '_UKCA_Nd';
            case 'vs Rain Rate 0.05'
                var_filter = RainRate_0pt05_PI_ALL;
                Nd_type_str = '_vs_Rain_Rate_0pt05';
                xlab_str = ['Cloud Base Rain Rate in Volc OFF run (kg m^{-2} s^{-1})'];
            case 'vs Rain Rate 0.01'
                var_filter = RainRate_0pt01_PI_ALL;
                Nd_type_str = '_vs_Rain_Rate_0pt01';
                xlab_str = ['Cloud Base Rain Rate in Volc OFF run (kg m^{-2} s^{-1})'];
        end
        
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['\DeltaLWP ' y_units_str];
        
    case 'dLWP vs SO2'
        var_test = LWP_PD_ALL - LWP_PI_ALL;
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc ON run (kg m^{-2})'];
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['\DeltaLWP (Volcano ON minus OFF) ' y_units_str];
        
    case 'LWP_PD vs SO2'
        var_test = LWP_PD_ALL;
        
        
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc ON run (kg m^{-2})'];
        
        
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['LWP Volcano ON ' y_units_str];
        
    case 'LWP_PI vs SO2'
        var_test = LWP_PI_ALL;
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc ON run (kg m^{-2})'];
        
        
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['LWP Volcano OFF ' y_units_str];
        
    case 'dLWP vs SO2'
        var_test = LWP_PD_ALL - LWP_PI_ALL;
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc ON run (kg m^{-2})'];
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['\DeltaLWP (Volcano ON minus OFF) ' y_units_str];
        
    case 'Nd_PD vs SO2'
        var_test = Nd_PD_ALL/1e6;
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc ON run (kg m^{-2})'];
        
        y_units_str = '(cm^{-3})';
        ylab_str = ['N_d (Volcano ON) ' y_units_str];
        
        
    case 'Nd_PI vs SO2'
        var_test = Nd_PI_ALL/1e6;
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc OFF run (kg m^{-2})'];
        
        y_units_str = '(cm^{-3})';
        ylab_str = ['N_d (Volcano OFF) ' y_units_str];  
        
    case 'dNd vs SO2'
        var_test = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
        
        var_filter = SO2_col_PD_ALL;
        Nd_type_str = '_vs_SO2';
        xlab_str = ['Column SO_2 in Volc OFF run (kg m^{-2})'];
        
        y_units_str = '(cm^{-3})';
        ylab_str = ['\DeltaN_d (Volcano ON minus OFF) ' y_units_str];     
        
    case 'LWP_PD vs Nd_PD'
        var_test = LWP_PD_ALL;               
        var_filter = (Nd_PD_ALL)/1e6;
        xlab_str = ['N_{d ' Nd_type ' volc ON} (cm^{-3})'];  
        
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['LWP Volcano ON ' y_units_str];  
        
    case 'LWP_PI vs Nd_PD'
        var_test = LWP_PI_ALL;               
        var_filter = (Nd_PD_ALL)/1e6;
        xlab_str = ['N_{d ' Nd_type ' volc ON} (cm^{-3})'];  
        
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['LWP Volcano OFF ' y_units_str];  
        
    case 'LWP_PI vs Nd_PI'
        var_test = LWP_PI_ALL;               
        var_filter = (Nd_PI_ALL)/1e6;
        xlab_str = ['N_{d ' Nd_type ' volc OFF} (cm^{-3})'];  
        
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['LWP Volcano OFF ' y_units_str];  
        
        
  case 'dLWP vs Nd_PD'
        var_test = LWP_PD_ALL - LWP_PI_ALL;               
        var_filter = (Nd_PD_ALL)/1e6;
        xlab_str = ['N_{d ' Nd_type ' volc ON} (cm^{-3})'];          
        
        y_units_str = '(g m^{-2})';
        ylab_str = ['\DeltaLWP ' y_units_str];          
        
    case 'dNd vs dNaccum'
        
        var_test = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
        var_filter = accum_number_z3000_PD_ALL - accum_number_z3000_PI_ALL;
        
        xlab_str = ['\DeltaN_{accum} (m^{-2})'];
        y_units_str = '(cm^{-3})';
        ylab_str = ['\DeltaN_{d ' Nd_type '} ' y_units_str];
        
        
    case 'dCFproxy_LWP vs dNd'
        thresh_LWP_CF = 5;
        
        var_test_01 = LWP_PI_ALL;
        var_test_02 = LWP_PD_ALL;
        
        CF_PI = zeros(size(LWP_PI_ALL));
        CF_PD = zeros(size(LWP_PD_ALL));
        
        %inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);
        %var_test_01(inan)=NaN;
        %var_test_02(inan)=NaN;
        
        clear N_PI N_PD leg_str
        for i=1:length(thresh_LWP_CF)
            %N_PI(i) = length(find(var_test_01 >= thresh_LWP_CF(i)));
            %N_PD(i) = length(find(var_test_02 >= thresh_LWP_CF(i)));
            CF_PI( find(LWP_PI_ALL >= thresh_LWP_CF(i)) ) = 1;
            CF_PD( find(LWP_PD_ALL >= thresh_LWP_CF(i)) ) = 1;
        end
        
        %Ntot = length(LWP_PI_ALL(:)) - length(inan);
        %CF_PI = N_PI/Ntot; %leg_str{1} = 'Volcano OFF';
        %CF_PD = N_PD/Ntot; %leg_str{2} = 'Volcano ON';
        
        var_test = CF_PD - CF_PI;
        
        Nd_type = 'CASIM';
        %Nd_type = 'UKCA';
        %Nd_type = 'vs Rain Rate 0.05';
        %Nd_type = 'vs Rain Rate 0.01';
        
        xlab_str = ['\DeltaN_{d ' Nd_type '} (cm^{-3})'];
        
        switch Nd_type
            case 'CASIM'
                var_filter = (Nd_PD_ALL - Nd_PI_ALL)/1e6;
            case 'UKCA'
                var_filter = (Nd_UKCA_PD_ALL - Nd_UKCA_PI_ALL)/1e6;
                Nd_type_str = '_UKCA_Nd';
            case 'vs Rain Rate 0.05'
                var_filter = RainRate_0pt05_PI_ALL;
                Nd_type_str = '_vs_Rain_Rate_0pt05';
                xlab_str = ['Cloud Base Rain Rate in Volc OFF run (kg m^{-2} s^{-1})'];
            case 'vs Rain Rate 0.01'
                var_filter = RainRate_0pt01_PI_ALL;
                Nd_type_str = '_vs_Rain_Rate_0pt01';
                xlab_str = ['Cloud Base Rain Rate in Volc OFF run (kg m^{-2} s^{-1})'];
        end
        
        
        y_units_str = '';
        ylab_str = ['\DeltaCF_{LWP>5} ' y_units_str];
        
    case 'dCF_subgrid vs dNd'        
        var_test = f1_orig - f0_orig;
        var_filter = N1_orig - N0_orig;
        
        
        xlab_str = ['\DeltaN_{d} (cm^{-3})'];
        y_units_str = '';
        ylab_str = ['\DeltaCF_{subgridCF} ' y_units_str];  
        
    case 'dLWP_subgrid vs dNd'
        var_test = W1_orig - W0_orig;
        var_filter = N1_orig - N0_orig;
        
        
        xlab_str = ['\DeltaN_{d} (cm^{-3})'];
        y_units_str = '';
        ylab_str = ['\DeltaLWP_{ic subgridCF} ' y_units_str];          
        
    case 'LWP_subgrid PI vs dNd'
        var_test = W0_orig;
        var_filter = N1_orig - N0_orig;        
        
        xlab_str = ['\DeltaN_{d} (cm^{-3})'];
        y_units_str = 'g m^{-2}';
        ylab_str = ['Volc OFF LWP_{ic subgridCF} ' y_units_str];  
        
    case 'LWP_subgrid PD vs dNd'
        var_test = W1_orig;
        var_filter = N1_orig - N0_orig;        
        
        xlab_str = ['\DeltaN_{d} (cm^{-3})'];
        y_units_str = 'g m^{-2}';
        ylab_str = ['Volc OFF LWP_{ic subgridCF} ' y_units_str];  
        
     case 'LWP_subgrid PI vs Nd PI'
        var_test = W0_orig;
        var_filter = N0_orig;        
        
        xlab_str = ['Volc OFF N_{d} (cm^{-3})'];
        y_units_str = '(g m^{-2})';
        ylab_str = ['Volc OFF LWP_{ic subgridCF} ' y_units_str];  
        
     case 'LWP_subgrid PD vs Nd PD'
        var_test = W1_orig;
        var_filter = N1_orig;        
        
        xlab_str = ['Volc ON N_{d} (cm^{-3})'];
        y_units_str = '(g m^{-2})';
        ylab_str = ['Volc ON LWP_{ic subgridCF} ' y_units_str]; 
        
    case 'CF_subgrid PD vs Nd PD'
        var_test = f1_orig;
        var_filter = N1_orig;        
        
        %var_test = f1;
        %var_filter = N1;        
        
        xlab_str = ['N_{d Volc ON} (cm^{-3})'];
        y_units_str = '';
        ylab_str = ['Cloud Fraction_{subgrid Volc ON} ' y_units_str];  
        
        prc_vals_for_bins = [0:5:100];
        
end

%% Filtering.

%LWP
inan = find(LWP_PI_ALL<thresh_LWP | LWP_PD_ALL<thresh_LWP);
var_test(inan)=NaN; var_filter(inan)=NaN; %Remove the filter values too to ensure even bin
%sampling when using percentiles

switch SO2_filter_type
    case 'delta'
        inan = find(SO2_col_PD_ALL - SO2_col_PI_ALL < thresh_SO2);  
        so2_str = '\DeltaSO2';
    case 'Volc ON'
        inan = find(SO2_col_PD_ALL < thresh_SO2);
        so2_str = 'volc on';
end
var_test(inan)=NaN; var_filter(inan)=NaN;


%CF
iCF_filter = 0;
if iCF_filter==1
    thresh_CF = 0.95; inan = find(low_CF_PD_ALL<thresh_CF | low_CF_PI_ALL<thresh_CF);
    thresh_CF = -1; inan = find(low_CF_PD_ALL<thresh_CF | low_CF_PI_ALL<thresh_CF);
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    CF_filter_str = [', CF threshold = ' num2str(thresh_CF)];
else
    CF_filter_str = '';
end



lon_rep = repmat(gcm_Plon2D_UM,[1 1 size(var_test,3)]);
inan=find(lon_rep<LON_filter(1) | lon_rep>=LON_filter(2));
var_test(inan)=NaN; var_filter(inan)=NaN;

%Time filtering
itime_filter = 0;
if itime_filter==1
    time_t0 = time_out-time_out(1); %time in days
    time_rep = repmat(time_t0(:),[1 size(var_test,1) size(var_test,2)]);
    time_rep = permute(time_rep,[2 3 1]);
    
    thresh_time = [0 1]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
    %thresh_time = [1 2]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
    %thresh_time = [2 3]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
    thresh_time = [3 4]; inan = find(time_rep < thresh_time(1) | time_rep >= thresh_time(2));
    var_test(inan)=NaN; var_filter(inan)=NaN;
    
    days_str = ['days=' num2str(thresh_time(1)) ' to ' num2str(thresh_time(2))];
else
    days_str='';
end

%Make bins using percentiles

%var_filter_bin_edges = [-50:10:800];
%var_filter_bin_edges = prctile(var_filter(:),[0:0.5:100]);
%var_filter_bin_edges = prctile(var_filter(:),);
%var_filter_bin_edges = prctile(var_filter(:),[0:10:100]);


var_filter_bin_edges = prctile(var_filter(:),prc_vals_for_bins);


title_str=[''];

opts=[];
[y,std_dev,N,ymean_overall,yn_overall,ystd_overall]=Hawaii_binned_by_filter_var(var_test,var_filter,var_filter_bin_edges,title_str,xlab_str,ylab_str,y_units_str,opts);
title_str=['LWP threshold = ' num2str(thresh_LWP) ' g m^{-2}, SO_2 threshold = ' num2str(thresh_SO2,'%1.0e')...
    SO2_filter_type ', ' CF_filter_str ', ' num2str(LON_filter(1)) ' to ' num2str(LON_filter(2)) '^{o}E'...
    ', ' days_str...
    ', minN=' num2str(min(N))];

std_err = std_dev./sqrt(N);
yplot = y;
%yplot(N<=2 | std_err>10 | isnan(std_err)==1)=NaN;
%yplot(N<=200 | std_err>10 | isnan(std_err)==1)=NaN;

figure('color','w');
%plot(thresh_Nd_multi,yvals_timemean,'linewidth',3);
mid_points = 0.5*( var_filter_bin_edges(1:end-1) + var_filter_bin_edges(2:end) );
plot(mid_points,yplot,'bo-','linewidth',3);
xlabel(xlab_str);
ylabel(ylab_str);
title(title_str);
fontsize_figure(gcf,gca,18);
grid on

%%
savefile = [savedir_date run_set ' ' plot_type ' ' title_str Nd_type_str];
save(savefile,'mid_points','yplot','LON_filter','thresh_LWP','thresh_SO2','thresh_CF','xlab_str',...
    'y_units_str','ylab_str','-V7.3');

%    loadfile = [savedir_date 'Blending option=3, volcano starting 12 UTC, no orography, adjusted emission height LWP threshold = 50 g m^{-2}, SO_2 threshold = 1e-05, -160 to -155^{o}E, , minN=1372.mat'];
%    dat = load(loadfile);

% %% Restrict region to threshold DeltaNd - filtered using time mean
%     thresh_Nd = 100; %per cc
%     dNd_map = (Nd_PD_map - Nd_PI_map)/1e6;
%     %dNd_map2 = repmat(dNd_map
%
%
%     iplume = find(dNd_map >= thresh_Nd);  %1d linear indices for 2D map.
%     [ix,iy]=ind2sub(size(dNd_map),iplume);
%     ix2=repmat(ix,[nT 1]);
%     iy2=repmat(iy,[nT 1]);
%     iz=[1:nT]';
%     iz2=(repmat(iz,[1 length(ix)]))';
%     iz2=iz2(:);
%     ii = sub2ind(size(dLWP),ix2,iy2,iz2);
%
%     dLWP_plume = NaN*ones(size(dLWP));
%     dLWP_plume(ii) = dLWP(ii);
%     dom_mean_dLWP_plume = meanNoNan(meanNoNan(dLWP_plume,1),1);
%
%     dNd_map2 = NaN*ones(size(dNd_map));
%     dNd_map2(iplume)=dNd_map(iplume);
%     qpcolor(dNd_map2);
%     caxis([0 400]);
%     title('\DeltaN_d');
%
%     figure
%     plot(time,dom_mean_dLWP_plume,'linewidth',3);
%     datetick('x','dd');
%
%
%     %legend(leg_strs,'location','northwest');
%     xlabel('Time');
%     %ylabel('SW surface forcing (W m^{-2})');
%     ylabel('\DeltaLWP (g m^{-2})');
%     set(gca,'ylim',[-30 50])
%     fontsize_figure(gcf,gca,18);
%     grid on
%     title(['\DeltaN_d thresh=' num2str(thresh_Nd) ' cm^{-3}, filtered with time-av N_d']);


