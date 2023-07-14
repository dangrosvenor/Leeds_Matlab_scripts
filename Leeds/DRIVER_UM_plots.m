% 12th nov time difference to LST :-
% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)


%LWP
plot_case_UM = '26thOct2008 POC case - timeseries of LWP from Moose NetCDF output (25Nov2014)';

%% For internal seminar
    %revised location for box (seemed to be plotting for the NE corner for
    %the AGU plots.... and different runs (ones shown in AGU talk, since
    %previous timeseries were for older runs)

plot_case_UM = '26thOct2008 POC case - timeseries of LWP from Moose NetCDF output (08May2015)';     
plot_case_UM = '26thOct2008 POC case - PDFs of CF from Moose NetCDF output (08May2015)'
%plot_case_UM = '26thOct2008 POC case - PDFs of LWP from Moose NetCDF output (08May2015)'    

plot_case_UM = '12thNov2008 Boutle case - PDFs of LWP from Moose NetCDF output (09Jul2015)'  

%% Timeseries
%RWP
%plot_case_UM = '26thOct2008 POC case - timeseries of RWP from Moose NetCDF output (25Nov2014)';
%plot_case_UM = 'Boutle_case_12Nov2008_RWP_timser_20150727.m';
%plot_case_UM = 'Boutle_case_12Nov2008_RWP_timser_20161121.m'; %Updated version that uses the UM_case_select_runs

%LWP
%plot_case_UM = '12thNov2008 Boutle case - PDFs of ship LWP vs model (10 min averages) (13Aug2015)';
%plot_case_UM = '12thNov2008 Boutle case - timeseries of LWP from Moose NetCDF output (29Jul2015)';


%% Maps
% plot_case_UM = '12thNov2008 Boutle case - maps of LWP from Moose NetCDF output (26Nov2014)';    
% plot_case_UM = '12thNov2008 Boutle case - map of Werner parameter vs ship PDFs';
% plot_case_UM = 'Boutle_case_12Nov2008_Nd_maps_GOES_20151112';
% %plot_case_UM = 'Boutle_case_12Nov2008_Nd_maps_20151112';
% plot_case_UM = 'Boutle_case_12Nov2008_CF_maps_20151112'; %Used these for Iceland cases too
% plot_case_UM = 'Boutle_case_12Nov2008_SWLW_maps_CERES_20151112';
% plot_case_UM = 'Boutle_case_12Nov2008_LWP_maps_20141126T041216'; %LWP map for UM
% plot_case_UM = 'Boutle_case_12Nov2008_orog_SST_maps_20141126T041216.m'; %Combined orography and SST map with domain region
% plot_case_UM = 'Boutle_case_12Nov2008_SHIP_LWP_PDFs_20150805T151800.m'; %Maps of Werner metric vs ship LWP PDFs.
% 
%% PDFs
% %plot_case_UM = '12thNov2008 Boutle case - PDFs of Nd (11Nov2015)';
% plot_case_UM = 'Boutle_case_12Nov2008_SWLW_PDFs_20151111';
% plot_case_UM = 'Boutle_case_12Nov2008_SHIP_LWP_PDFs_20150805T151800.m'; %Ship PDFs of LWP vs model

%% Joint histograms
%plot_case_UM = 'Boutle_case_12Nov2008_NdLWP_joint_PDFs_20151111';
%plot_case_UM = 'Boutle_case_12Nov2008_Height_dBZ_joint_PDFs_20151111';

%SWITCH

switch plot_case_UM
    
%% Timeseries    

    case '26thOct2008 POC case - timeseries of LWP from Moose NetCDF output (08May2015)'
        POC_26Oct2008_LWP_timser_20150508          
    case '26thOct2008 POC case - timeseries of LWP from Moose NetCDF output (25Nov2014)'
        POC_26Oct2008_LWP_timser_20141125T032943  
    case '26thOct2008 POC case - timeseries of RWP from Moose NetCDF output (25Nov2014)'
        POC_26Oct2008_RWP_timser_20141125T032943  
    case '12thNov2008 Boutle case - timeseries of LWP from Moose NetCDF output (29Jul2015)'
        Boutle_case_12Nov2008_LWP_timser_20150727    
    case '12thNov2008 Boutle case - timeseries of Nd from Moose NetCDF output (29Jul2015)'        
%        Boutle_case_12Nov2008_Nd_timser_20151112.m
        Boutle_case_12Nov2008_Nd_timser_20161121.m      %New version to allow the use of UM_case_select_runs
 
        %RWP as of July, 2016
%        'Boutle_case_12Nov2008_RWP_timser_20150727.m'
        
% N.B. see --- read_GOES_vocals_netcdf_files_MULTI.m ---
% for reading in muliple files and saving the data in a .mat file        
        
%% Maps        
    case '12thNov2008 Boutle case - maps of LWP from Moose NetCDF output (26Nov2014)'
%        Boutle_case_12Nov2008_LWP_maps_20141126T041216 
        Boutle_case_12Nov2008_LWP_maps_20141126T041216_Fig_v1        
    case '12thNov2008 Boutle case - map of Werner parameter vs ship PDFs'
        Boutle_case_12Nov2008_Werner_map_20151126
        
    % For GOES from the .mat file (top half of UM domain only) :-    
    case 'Boutle_case_12Nov2008_Nd_maps_GOES_20151112' %also LWP maps
        Boutle_case_12Nov2008_Nd_maps_GOES_20151112
        
    %GOES maps for more general cases (whole fields, etc.)    
    case 'POC_26Oct2008_GOES_general_maps_20141125T032943'     
        POC_26Oct2008_GOES_general_maps_20141125T032943
        
    case 'Boutle_case_12Nov2008_Nd_maps_20151112'    
          Boutle_case_12Nov2008_Nd_maps_20151112
          
    case 'Boutle_case_12Nov2008_max_dBZ_maps_20151112'    
          Boutle_case_12Nov2008_max_dBZ_maps_20151112          
        
%% PDFs        
    case '26thOct2008 POC case - PDFs of CF from Moose NetCDF output (08May2015)'
        POC_26Oct2008_CF_0pt25deg_PDFs_20141125T032943        
    case '12thNov2008 - PDFs of CF from Moose NetCDF output (08May2015)'
        POC_26Oct2008_CF_0pt25deg_PDFs_20141125T032943       
    case 'Boutle_12Nov2008_CF_0pt25deg_PDFs_20151215.m'
        Boutle_12Nov2008_CF_0pt25deg_PDFs_20151215.m
    case '26thOct2008 POC case - PDFs of LWP from Moose NetCDF output (08May2015)'
        POC_26Oct2008_LWP_PDFs_20141125T032943
    case '12thNov2008 Boutle case - PDFs of LWP from Moose NetCDF output (09Jul2015)'
%        Boutle_case_12Nov2008_LWP_PDFs_20150709T172500
        Boutle_case_12Nov2008_LWP_PDFs_20160216        
        
        %PDF based on ship data :-
    case '12thNov2008 Boutle case - PDFs of ship LWP vs model (10 min averages) (13Aug2015)'
        Boutle_case_12Nov2008_SHIP_LWP_PDFs_20150805T151800   
        
    case '12thNov2008 Boutle case - PDFs of Nd (11Nov2015)'
        Boutle_case_12Nov2008_Nd_PDFs_20151111.m
        
    case '12thNov2008 Boutle case - PDFs of SWLW (28Jan2016)'
        Boutle_case_12Nov2008_SWLW_PDFs_20151111.m        
        
%% Joint PDFs
    case 'Boutle_case_12Nov2008_NdLWP_joint_PDFs_20151111'
        Boutle_case_12Nov2008_NdLWP_joint_PDFs_20151111
        
%% Or just execture the name given in plot_case_UM
    otherwise
        eval(plot_case_UM)
        
end