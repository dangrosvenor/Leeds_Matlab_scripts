% Box for regional averages
offset_lon = 0; %for text in box
offset_lat = 0;

inset_ax_pos = [0.1636    0.7186    0.2500    0.2500];
inset_ax_pos_cf = [0.1636    0.7186    0.2500    0.2500];
inset_ax_pos_sw = [0.1636    0.7186    0.2500    0.2500];

ylims_main = [45 230];
ylims_main_Nd = [45 230];
ylims_main_cf = [0.25 0.55];
ylims_main_totcf = [0.62 0.75];
ylims_main_totcf = [0.63 0.74];
ylims_main_sw = [90 130];
ylims_inset_sw = [90 130];

iplot_legend=1;
iplot_legend_cf=1;
iplot_legend_sw=1;

sigfig_str = '%.3f'; %Default

box_region_str='unlabelled';
switch box_region
    case 'a'
        LAT_val = [30 42]; LON_val =[-68 -52]; %
    case '1'
        LAT_val = [30 40]; LON_val =[-73 -63]; %box in plots sent to Ken
        %box in plots sent to Ken
        box_region_str = 'US outflow';
        box_region_num_str = '1';
        
        ylims_main_sw = [90 130];
        ylims_main_cf = [0.26 0.56];
        
        ylims_inset=[112 161];
        ylims_inset_cf=[0.31 0.43];        
        ylims_inset_sw = [89 117];
        
        inset_ax_pos_cf = [0.1636    0.6813    0.2500    0.2500];
        
    case '2'
        LAT_val = [30 45]; LON_val =[-73 -40]; %wider box region in NW Atlantic
        box_region_str = '2';
    case '3'
        LAT_val = [30 40]; LON_val =[-20 -10]; %Region SW of Spain down African coast
        box_region_str = 'Europe outflow';
        box_region_num_str = '3';
        inset_ax_pos = [0.1636    0.5634    0.2500    0.2500];
        inset_ax_pos_cf = [0.1636    0.7342    0.2500    0.2500];
        
        ylims_main_sw = [65 105];
        ylims_main_cf = [0.26 0.56];
        
        ylims_inset=[90 162];
        ylims_inset_cf=[0.26 0.53];
        ylims_inset_sw = [71 91];
        
        iplot_legend_cf=0;
        
        
        
    case '4'
        LAT_val = [18 60]; LON_val = [-75 0]; % "main" part of NA - need to exclude land here
        %box_region_str = 'N Atlantic';
        box_region_str = 'N. Atlantic basin';
        box_region_num_str = '';
        offset_lon = 10; %for text in box
        
        ylims_main_cf = [0.37 0.45];       
        %ylims_main_sw = [85 110];
        ylims_main_sw = [86 103];
        ylims_main_Nd = [50 130];
        
        ylims_inset=[81 105];
        ylims_inset_cf=[0.37 0.44];
        ylims_inset_sw = [86 96];                
        
        inset_ax_pos = [0.1636    0.5792    0.2500    0.2500];
        inset_ax_pos_cf = [0.1636    0.6876    0.2500    0.2500];
        
        iplot_legend_cf=1;
        
    case '5'
        LAT_val = [50 60]; LON_val = [-10 3]; % UK region
        box_region_str = 'N Atlantic';
        box_region_str = '5';    
        
    case '6'
        LAT_val = [50 60]; LON_val = [-20 -10]; % UK region
        box_region_str = 'N Atlantic';
        box_region_str = '6';
        
    case '7'
        LAT_val = [50 60]; LON_val = [-40 -20]; % UK region
        box_region_str = 'N Atlantic';
        box_region_str = '7';        
        
    case '8'
        LAT_val = [30 40]; LON_val = [-45 -35]; % UK region
        box_region_str = 'Mid-Atlantic';
        %box_region_str = '8';  
        %offset_lon = 3; %for text in box
        
        ylims_main_sw = [75 115];
        ylims_main_cf = [0.26 0.56];
        ylims_main_Nd = [40 120];
        
        ylims_inset=[64 93];
        ylims_inset_cf=[0.26 0.44];
        ylims_inset_sw = [73 100];
        
        inset_ax_pos = [0.1636    0.5292    0.2500    0.2500];
        inset_ax_pos_cf = [0.1636    0.6255    0.2500    0.2500];
        box_region_num_str = '2';
        
        
        
     case '9'
        LAT_val = [50 60]; LON_val = [-60 -40]; %Canada region
        box_region_str = 'N Atlantic';
        box_region_str = '9';  
        
    case '10'
        LAT_val = [22 42]; LON_val = [-80 -43]; %
        box_region_str = 'US outflow SW trend';
        box_region_str = '10';
        
     case '11'
        LAT_val = [24 42]; LON_val = [-80 -53]; %
        box_region_str = 'US outflow SW trend';
        box_region_str = '11';
        
        if iplot_individual_ens==1
            ylims_main_sw = [85 105];
        else
            ylims_main_sw = [86 98];
            ylims_main_sw = [87 104]; %1900 onwards plot
        end    
           
        
end
%for label for plot_box_on_map
%box_text = box_region_str;
box_text = box_region_num_str;


switch var_ukesm
    case {'calipso_low_cloud_amount','calipso_total_cloud_amount'}
        sigfig_str = '%.1e';
end
