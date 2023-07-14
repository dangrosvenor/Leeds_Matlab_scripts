% Box for regional averages

if ~exist('iplot_individual_ens')
    iplot_individual_ens=0;
end

loc='NorthWest'; %default

offset_lon = 0; %for text in box
offset_lat = 0;

inset_ax_pos = [0.1636    0.7186    0.2500    0.2500];
inset_ax_pos_cf = [0.1636    0.7186    0.2500    0.2500];
inset_ax_pos_sw = [0.1636    0.7186    0.2500    0.2500];

ylims_main = [45 230];
ylims_main_Nd = [45 230];
ylims_main_cf = [0.25 0.55];
ylims_main_sw = [90 130];
ylims_main_SO2 = [0 1.2];
ylims_main_ts = [291 293]; %surface temperature
ylims_main_lwp = [70 85];
ylims_main_lwpic = [100 120];

ylims_inset_sw = [90 130];

 ylims_main_totcf = [0.595 0.78];

iplot_legend=1;
iplot_legend_cf=1;
iplot_legend_sw=1;

box_region_num_str = ''; %default
box_region_str='unlabelled';
switch box_region
    case '0'
        LAT_val = [-1e9 1e9]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Global';
        ylims_main_sw = [93 100];
    case '00'
        LAT_val = [-60 60]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Global 60S to 60N';
        switch land_ocean
            case 'land only'
                ylims_main_sw = [109 118];
                ylims_main_ts = [289.8 292.5];
            case 'ocean only'
                ylims_main_sw = [92 100];
                ylims_main_ts = [292 294]; %surface temperature
            case 'land+ocean'                
                ylims_main_sw = [97 104];
        end
    case '01'
        LAT_val = [-55 60]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Global -55 to 60';
        ylims_main_sw = [92 98];
    case '02'
        LAT_val = [-50 55]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Global -50 to 55';
        ylims_main_sw = [90 97];    
    case '03'
        LAT_val = [-60 -40]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Southern Ocean 40 to 60^{o} S';
        ylims_main_sw = [107 117];      
        loc='SouthWest';
    case '04'
        LAT_val = [-80 -60]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Antarctic sea-ice region';
        ylims_main_sw = [107 117];      
        %loc='SouthWest';     
    case '05'
        LAT_val = [60 1e9]; LON_val =[-1e9 1e9]; %
        box_region_str = 'Arctic region'; 
        ylims_main_sw = [107 117];      
        %loc='SouthWest';  
        
    case '06'
        LAT_val = [-40 0]; LON_val ={[-1e9 -70],[150 1e9]}; %Will need to deal with wraparound if want to go further west than -180.
        box_region_str = 'South Pacific region';
        ylims_main_sw = [107 117];      
        %loc='SouthWest';  
        
    case '07'
        LAT_val = [-40 0]; LON_val =[30 115]; %
        box_region_str = 'Indian Ocean region';
        ylims_main_sw = [107 117];      
        %loc='SouthWest'; 
        
    case '08'
        LAT_val = [-40 0]; LON_val =[-55 20]; %
        box_region_str = 'South Atlantic region';
        ylims_main_sw = [107 117];      
        %loc='SouthWest';  
        
    case '09'
        LAT_val = [0 60]; LON_val =[-80 -10]; %
        box_region_str = 'North Atlantic region';
        ylims_main_sw = [107 117];      
        %loc='SouthWest';    
        
    case '10'
        LAT_val = [0 60]; LON_val ={[-1e9 -95],[125 1e9]}; %
        box_region_str = 'North Pacific region';
        ylims_main_sw = [107 117];
        %loc='SouthWest';
        
    case 'L01'
        LAT_val = [-55 15]; LON_val =[-85 -35]; %
        box_region_str = 'South America region';
        ylims_main_sw = [107 117];
        %loc='SouthWest';  
        
    case 'L02'
        LAT_val = [15 75]; LON_val =[-170 -55]; %
        box_region_str = 'North America region';
        ylims_main_sw = [107 117];
        %loc='SouthWest';
        
     case 'L03'
        LAT_val = [-35 35]; LON_val =[-20 60]; %
        box_region_str = 'Africa region';
        ylims_main_sw = [107 117];
        %loc='SouthWest'; 
        
     case 'L04'
        LAT_val = [35 80]; LON_val =[-10 1e9]; %
        box_region_str = 'Europe and Northern Asia';
        ylims_main_sw = [107 117];
        %loc='SouthWest'; 
        
    case 'L05'
        LAT_val = [0 35]; LON_val =[60 140]; %
        box_region_str = 'Southern Asia';
        ylims_main_sw = [107 117];
        %loc='SouthWest';         
        
        
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
        %ylims_main_totcf = [0.63 0.78];
        ylims_main_totcf = [0.61 0.71];
        %ylims_main_lwpic = [103 122];
        ylims_main_lwpic = [106 122];
        ylims_main_lwp = [71 83];
        %ylims_main_sw = [85 110];
        %ylims_main_Nd = [50 130];
        ylims_main_Nd = [50 130];
        ylims_main_rsutcs = [38 45];
        ylims_main_od550tot = [0.05 0.25];
        
        %ylims_main_ts = [289.5 292.5];
        ylims_main_ts = [291.5 294];
        
        if iplot_individual_ens==1
            ylims_main_sw = [86 102];
        else
            switch model_str   
                case 'DAMIP_hist-aer'
                    ylims_main_sw = [84 102];
                case 'DAMIP_hist-GHG' %For the offline flux calculation plots
                    ylims_main_sw = [81 96];
                case 'DAMIP_hist-nat'
                    ylims_main_sw = [84 96];                    
                otherwise
                    ylims_main_sw = [86 98];
                    ylims_main_sw = [87 104]; %1900 onwards plot
                    %ylims_main_sw = [87 104]; %1900 onwards plot
                    ylims_main_sw = [84 100];
            end
            ylims_main_swtoa_calc = [92 105]; 
        end
        if iadd_DAMIP==1
            ylims_main_sw = [-5.5 12]; 
            %ylims_main_totcf = [-4 7]*1e-2; %perturbation plot
            ylims_main_totcf =[-0.045 0.06];
            %ylims_main_lwpic = [-5.5 7.5];
            ylims_main_lwpic = [-6 9];
            %ylims_main_lwp = [-4 6];
            ylims_main_lwp = [-6 6];
            %ylims_main_Nd = [-3 40];
            ylims_main_Nd = [-3 43];
            ylims_main_rsutcs = [-1 7];
            ylims_main_ts = [-1.7 1.7]; 
            ylims_main_od550tot = [-0.03 0.1];
        end
        if exist('iadd_AerChemMIP') & iadd_AerChemMIP==1
            ylims_main_sw = [-7 11]; 
            ylims_main_ts = [-1.6 1.9]; 
            ylims_main_rsutcs = [-1.0 9.0]; 
            ylims_main_lwpic = [-12 18]; 
            ylims_main_totcf = [-0.07 0.1]; 
            ylims_main_Nd = [-10 60];
            ylims_main_od550tot = [-0.05 0.1];
            
        end
        
        if exist('i_CMIP6_multimodel') & i_CMIP6_multimodel==1
            ylims_main_Nd = [0 230];
            ylims_main_Nd = [0 310];
            ylims_main_reff = [4.5 7.5];
        end
        
        
        
        
        
        ylims_inset=[81 105];
        ylims_inset_cf=[0.37 0.44];
        ylims_inset_sw = [86 96];                
        
        
        inset_ax_pos = [0.1636    0.5792    0.2500    0.2500];
        inset_ax_pos_cf = [0.1636    0.6876    0.2500    0.2500];
        
        iplot_legend_cf=0;
        
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
        box_region_str = 'US outflow';
        box_region_str = '10';
        
        if iplot_individual_ens==1
            ylims_main_sw = [85 105];
        else
            ylims_main_sw = [86 98];
            ylims_main_sw = [87 104]; %1900 onwards plot
        end
        
    case '11'
        LAT_val = [26 42]; LON_val = [-80 -56]; %Also used for Matt C natural labs paper
        box_region_str = 'US outflow';
        %box_region_str = '10';
        
        %if iplot_individual_ens==1
            switch season
                case 'DJF'
                    ylims_main_sw = [75 120];
                case 'MAM'
                    ylims_main_sw = [100 155];
                case 'JJA'
                    ylims_main_sw = [90 125];
                case 'SON'
                    ylims_main_sw = [74 110];
                otherwise
                    ylims_main_sw = [90 116];
                    ylims_main_totcf = [0.595 0.78];
                    ylims_main_lwp = [80 110];
                    inset_ax_pos = [0.1636    0.57    0.2500    0.2500];
                    loc = 'NorthEast';

            end
       % else
       %     ylims_main_sw = [86 98];
       %     ylims_main_sw = [87 104]; %1900 onwards plot
       % end
       
     case '12'
        LAT_val = [26 50]; LON_val = [-100 -70]; %
        box_region_str = 'US mainland';
        
        ylims_main_SO2 = [0 1.2]; %converted from kg/m2/s to tonnes/km2/yr
        
        
    case '13'
        LAT_val = [18 50]; LON_val = [-75 0]; % "main" part of NA - need to exclude land here
        %box_region_str = 'N Atlantic';
        box_region_str = 'N. Atlantic basin to 50N';
        box_region_num_str = '';
        offset_lon = 10; %for text in box
        
        ylims_main_cf = [0.37 0.45];
        %ylims_main_sw = [85 110];
        if iplot_individual_ens==1
            ylims_main_sw = [86 102];
        else
            ylims_main_sw = [86 98];
            ylims_main_sw = [87 104]; %1900 onwards plot
        end
        ylims_main_Nd = [50 130];
        
        ylims_inset=[81 105];
        ylims_inset_cf=[0.37 0.44];
        ylims_inset_sw = [86 96];                
        
        
        inset_ax_pos = [0.1636    0.5792    0.2500    0.2500];
        inset_ax_pos_cf = [0.1636    0.6876    0.2500    0.2500];
        
        iplot_legend_cf=0;
        
        
    case '14'
        LAT_val = [20 45]; LON_val =[-30 -7.5]; %
        %LAT_val = [35 40]; LON_val =[-30 -12.5]; %
        box_region_str = 'Eastern N. Atlantic';
        ylims_main_sw = [75 102];         
            
                   
    case '15'
        LAT_val = [-40 0]; LON_val =[-35 10]; %
        %LAT_val = [35 40]; LON_val =[-30 -12.5]; %
        box_region_str = 'Southern Atlantic';
        ylims_main_sw = [75 102];                  
        
    case '16'
        LAT_val = [0 20]; LON_val =[-35 10]; %
        box_region_str = 'Sahara';
        
    case '17'
        LAT_val = [-30 -10]; LON_val =[-110 -80]; %
        box_region_str = 'SE Pacific Zhou2016';  
        ylims_main_ts = [293 297]; %surface temperature
        
    case '18'
        LAT_val = [36 54]; LON_val = [-52 -10]; %
        box_region_str = 'Northern NA paper01'; %region used for 2020 1st ACSIS paper.
        
    case '19'
        LAT_val = [18 36]; LON_val = [-69 -20]; %
        box_region_str = 'Southern NA paper01';  
        
    case '20'
        LAT_val = [-10 80]; LON_val = [-85 30]; % Full NA region used for 2020 1st ACSIS paper.
        box_region_str = 'N Atlantic paper01';    
        
    case '21'
        LAT_val = [18 40]; LON_val = [-75 0]; % As for 1st SW paper March, 2021, but up to 40N.
        box_region_str = 'South N. Atlantic SW paper'; 
        if iadd_DAMIP==1
            ylims_main_sw = [-5.5 9]; 
            ylims_main_totcf = [-5 5]*1e-2; %perturbation plot
            ylims_main_lwpic = [-5.5 7.5];
            ylims_main_lwp = [-4 6];
            ylims_main_Nd = [-3 50];
        end
        
    case '22'
        LAT_val = [47 60]; LON_val = [-60 -45]; % 
        box_region_str = 'Sea-ice region Newfoundland SW paper'; 
        if iadd_DAMIP==1
            ylims_main_sw = [-5.5 9]; 
            ylims_main_totcf = [-5 5]*1e-2; %perturbation plot
            ylims_main_lwpic = [-5.5 7.5];
            ylims_main_lwp = [-4 6];
            ylims_main_Nd = [-3 50];
        end     
        
    case '23'
        LAT_val = [30 40]; LON_val = [-80 -65]; % 
        box_region_str = 'US high clear-sky SWTOA bias region SW paper'; 
        if iadd_DAMIP==1
            ylims_main_sw = [-5.5 9]; 
            ylims_main_totcf = [-5 5]*1e-2; %perturbation plot
            ylims_main_lwpic = [-5.5 7.5];
            ylims_main_lwp = [-4 6];
            ylims_main_Nd = [-3 50];
        end          

        
    case '24'
        LAT_val = [32.5 50]; LON_val = [-115 -80]; %Jon's US domain for AMOC work July 2021
        %box_region_str = 'N Atlantic';
        box_region_str = 'US Jon Robson';
        box_region_num_str = '';
        offset_lon = 10; %for text in box
        
        ylims_main_cf = [0.37 0.45];
        %ylims_main_totcf = [0.63 0.78];
        ylims_main_totcf = [0.61 0.71];
        %ylims_main_lwpic = [103 122];
        ylims_main_lwpic = [106 122];
        ylims_main_lwp = [71 83];
        %ylims_main_sw = [85 110];
        %ylims_main_Nd = [50 130];
        ylims_main_Nd = [50 130];
        ylims_main_rsutcs = [38 45];
        
        if iplot_individual_ens==1
            ylims_main_sw = [86 102];
        else
            switch model_str
                case 'DAMIP_hist-aer'
                    ylims_main_sw = [84 102];
                case 'DAMIP_hist-GHG' %For the offline flux calculation plots
                    ylims_main_sw = [81 96];
                case 'DAMIP_hist-nat'
                    ylims_main_sw = [84 96];                    
                otherwise
                    ylims_main_sw = [86 98];
                    ylims_main_sw = [87 104]; %1900 onwards plot
                    %ylims_main_sw = [87 104]; %1900 onwards plot
                    ylims_main_sw = [84 100];
            end
            ylims_main_swtoa_calc = [92 105]; 
        end
        if iadd_DAMIP==1
            ylims_main_sw = [-5.5 9]; 
            %ylims_main_totcf = [-4 7]*1e-2; %perturbation plot
            ylims_main_totcf =[-0.045 0.06];
            %ylims_main_lwpic = [-5.5 7.5];
            ylims_main_lwpic = [-6 9];
            %ylims_main_lwp = [-4 6];
            ylims_main_lwp = [-6 6];
            %ylims_main_Nd = [-3 40];
            ylims_main_Nd = [-3 43];
            ylims_main_rsutcs = [-1 7];
        end
        
        if i_CMIP6_multimodel==1
            ylims_main_Nd = [0 230];
            ylims_main_Nd = [0 310];
            ylims_main_reff = [4.5 7.5];
        end
        
        ylims_main_ts = [289.5 292.5];
        
        
        
        ylims_inset=[81 105];
        ylims_inset_cf=[0.37 0.44];
        ylims_inset_sw = [86 96];                
        
        
        inset_ax_pos = [0.1636    0.5792    0.2500    0.2500];
        inset_ax_pos_cf = [0.1636    0.6876    0.2500    0.2500];
        
        iplot_legend_cf=0;
       


        
end
%for label for plot_box_on_map
%box_text = box_region_str;
box_text = box_region_num_str;
