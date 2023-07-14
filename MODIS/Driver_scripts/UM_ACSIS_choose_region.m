function [LAT_val_DRIVER2, LON_val_DRIVER2, region_shortname] = UM_ACSIS_choose_region(region_choice)

switch region_choice
    case 'North Atlantic'
        %LAT_val_DRIVER2 = [-10 60]; LON_val_DRIVER2 = [-85 0];
        LAT_val_DRIVER2 = [-10 80]; LON_val_DRIVER2 = [-85 30]; %Looking further north to show Iceland
        region_shortname = 'NA';
    case 'North Atlantic SW paper'       
        LAT_val_DRIVER2 = [18 60]; LON_val_DRIVER2 = [-75 0]; %Need to exclude land too...
        region_shortname = 'NASW';   
     case 'Global'
        %LAT_val_DRIVER2 = [-10 60]; LON_val_DRIVER2 = [-85 0];
        LAT_val_DRIVER2 = [-1e9 1e9]; LON_val_DRIVER2 = [-1e9 1e9]; %Looking further north to show Iceland
        region_shortname = 'Global';    
    case 'North Atlantic further north'
        LAT_val_DRIVER2 = [-10 80]; LON_val_DRIVER2 = [-85 30];
        region_shortname = 'NA';        
    case 'Northern Hemisphere Yawen'
        LAT_val_DRIVER2 = [-15 55]; LON_val_DRIVER2 = [-150 40];
        region_shortname = 'NA';        
    case 'Northern Hemisphere Yawen, further north'
        LAT_val_DRIVER2 = [-15 85]; LON_val_DRIVER2 = [-150 40];
        region_shortname = 'NA';                
    case 'Southern NA' %Southern part of N. Atlantic where CF changes dominate forcing 
        %LAT_val_DRIVER2 = [25 42]; LON_val_DRIVER2 = [-69 -20]; %Region of mainly Nd contribution - v10.3 model
        LAT_val_DRIVER2 = [18 36]; LON_val_DRIVER2 = [-69 -20]; %Region of mainly Nd contribution - UKESM model
        region_shortname = 'SNA';
    case 'Northern NA' %Northern part of N. Atlantic where Nd changes dominate forcing 
        %LAT_val_DRIVER2 = [43 60]; LON_val_DRIVER2 = [-52 -10]; %Region of mainly Nd contribution - v10.3 model 
        LAT_val_DRIVER2 = [36 54]; LON_val_DRIVER2 = [-52 -10]; %Region of mainly Nd contribution - UKESM model
        region_shortname = 'NNA';
    case 'Rosenfeld VOCALS' %Pacific region west of VOCALS
        LAT_val_DRIVER2 = [-40 15]; LON_val_DRIVER2 = [-120 -70]; 
        region_shortname = 'RosVOC';
    case 'Rosenfeld ALL' 
        LAT_val_DRIVER2 = [-40 0]; LON_val_DRIVER2 = [-1e9 1e9];
        region_shortname = 'RosALL';
    case 'VOCALS CPT'
        LAT_val_DRIVER2 = [-30 -10]; LON_val_DRIVER2 = [-100 -80];
        region_shortname = 'VOC';
    case 'VOCALS coastal'
        LAT_val_DRIVER2 = [-24 -18]; LON_val_DRIVER2 = [-80 -72];  
        region_shortname = 'VOCcoast';
    case 'East of US SW trend'
         %LAT_val = [22 42]; LON_val = [-80 -43]; 
         % Data already cut down, so select all here
         LAT_val_DRIVER2 = [-1e9 1e9]; LON_val_DRIVER2 = [-1e9 1e9]; 
         region_shortname = 'USswTrend';

    case 'US continent Jon Robson AMOC'
         LAT_val_DRIVER2 = [32.5 50]; LON_val_DRIVER2 = [-115 -80]; 
         region_shortname = 'USamocJonR';

    otherwise
        region_shortname = 'unspec';
        LAT_val_DRIVER2 = [-1e9 1e9]; LON_val_DRIVER2 = [-1e9 1e9]; %Global - or just set irestrict_domain_DRIVER=0
        %LAT_val_DRIVER2 = [30 45]; LON_val_DRIVER2 = [-60 -10]; %N Atlantic region of high negative forcing ocean only
        %LAT_val_DRIVER2 = [30 50]; LON_val_DRIVER2 = [-60 -10]; %N Atlantic region of high negative forcing ocean only
        %LAT_val_DRIVER2 = [40 50]; LON_val_DRIVER2 = [-60 -10]; %Smaller forcing region since is narrower for runs with varying SSTs
        %LAT_val_DRIVER2 = [33 50]; LON_val_DRIVER2 = [-60 -10]; %Extending further south for runs with varying SSTs and wind only nudging (u-ay837 and 838)
        %LAT_val_DRIVER2 = [5 32]; LON_val_DRIVER2 = [-40 -10]; %Region of large cloud fraction change near Africa
        %LAT_val_DRIVER2 = [40 60]; LON_val_DRIVER2 = [-50 -10]; %Region of mainly Nd contribution.
        
        
end
