isave_highlight=0;
%do the plot in watervap first - then this just scales and saves it

time_highlight_path=[];% 
y_axis_type='';
x_axis_type='';
iytick_relabel=0;    
fsize=14;


clear plot_graph_case
istring=1;
%plot_graph_case{istring}='Temp'; istring=istring+1;
%plot_graph_case{istring}='3D flight track'; istring=istring+1;
%plot_graph_case{istring}='Vertical wind speed'; istring=istring+1;
%plot_graph_case{istring}='Altitude'; istring=istring+1;
%plot_graph_case{istring}='Aircraft T pert'; istring=istring+1;
%plot_graph_case{istring}='Potemp'; istring=istring+1;
%plot_graph_case{istring}='Pressure'; istring=istring+1;
%plot_graph_case{istring}='RHi hygrometer'; istring=istring+1;
%plot_graph_case{istring}='RHi humi'; istring=istring+1;
%plot_graph_case{istring}='RH hygrometer'; istring=istring+1;
%plot_graph_case{istring}='RH humi'; istring=istring+1;
plot_graph_case{istring}='LWC'; istring=istring+1;
%plot_graph_case{istring}='Total large ice number CIP'; istring=istring+1;
%plot_graph_case{istring}='Total number CAS'; istring=istring+1;
%plot_graph_case{istring}='Ice Number (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='Total number size range (Jonny)'; istring=istring+1; %or ice too
%plot_graph_case{istring}='Binned ice number size range (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='N crystals per bin ice number size range (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='Number of ice accepts (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='Round Number (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='Wind direction'; istring=istring+1;
%plot_graph_case{istring}='Wind speed'; istring=istring+1;
%plot_graph_case{istring}='Ice Mass (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='Ice Diameter (Jonny)'; istring=istring+1;
%plot_graph_case{istring}='Airspeed'; istring=istring+1;
%plot_graph_case{istring}='Estimated displacement'; istring=istring+1;
%plot_graph_case{istring}='Number aerosol 610-1030 nm'; istring=istring+1;
%plot_graph_case{istring}='Vertical speed of aircraft'; istring=istring+1;





if ~strcmp(flight_no,'Calibration')
    eval(['time_flt=time_flt' flight_no ';']);
    air_speed_type = 'CIP probe';
    airspeed_constant=0;
    savedir=['Y:\BAS_flights\flight' flight_no '\CAS_plots\']; 
else
    air_speed_type = 'constant';
    airspeed_constant=17;
    switch comp
        case 'uni'
            savedir=['c:/documents and settings/dan/My Documents\logbook\Antarctica\Flights and instruments_Feb2010\Manchester and BAS CAS comparisons\'];
        case 'UWchallenger'
            savedir='/home/disk/eos1/d.grosvenor/Ant_flights/plots/ALL_flights_paper/';
    end
end
        
switch flight_no
        case 'Calibration'
            clear times
        
        case '102'            

        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                
                
                %%%%%%%%%%%%
                 segment=10;  %for highlighting a section only
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[20+20/60+0/3600 21+00/60+0/3600];%
                    case 2
                        time_highlight_path=[21+7/60+0/3600 21+9/60+0/3600];%
                    case 3
                        time_highlight_path=[21+00/60+0/3600 21+7/60+0/3600];%
                    case 4
                        time_highlight_path=[21+34/60+0/3600 21+42/60+0/3600];%
                    case 5
                        time_highlight_path=[21+11/60+0/3600 21+12/60+0/3600];%                        
                    case 6
                        time_highlight_path=[21+18/60+0/3600 21+20/60+0/3600];%   
                    case 7
                        time_highlight_path=[20+20/60+0/3600 20+40/60+0/3600];%  
                    case 8
                        time_highlight_path=[20.2795 20.2934];  %flight 102 50-54 km
                    case 9
                        time_highlight_path=[20.785 20.795];
                    case 10
                        time_highlight_path=[20+54/60+20/3600 20+54/60+40/3600];
                    case 99
                        time_highlight_path=[];%                                                                        
                       
                        
                end
                
             xlims3D=[275 375];
             ylims3D=[300 400];   
             
             xlims3D=[280 340];
             ylims3D=[310 340];   
             
             xlims3D=[260 340];
             ylims3D=[310 340];  
                          
             xlims3D=[260 390];
             ylims3D=[310 400];  

            case 2
  
                times=[20+20/60 20+40/60]/24;   
%                times=[20+40/60 21+0/60]/24;                   
%                times=[21+0/60 21+20/60]/24;   
%                times=[21+30/60 21+42/60]/24;    
%                 times=[20.3791 20.4107]/24; %flight 102 70-80 km
%                times=[20.5122 20.5344]/24;   %flight 102 108-115 km
%                times=[20.6312 20.6514]/24;   %flight 102 146-152 km
%               times=[20.8143 20.8646]/24;   %flight 102 146-152 km
%               times=[21.2306 21.4508]/24;   %flight 102 300-350 km
%               times=[20.2795 20.2934]/24;  %flight 102 50-54 km
%               times=[20.7894 20.7902]/24;  %flight 102 3.22 km   
%               times=[20.75 20.82]/24;  %flight 102 3.22 km   
%               times=[20.5 21]/24;  %flight 102 3.22 km                  
               

                xlims3D=[250 400];
                ylims3D=[275 425];
                                
                xlims3D=[275 375];
                ylims3D=[300 400];
                
                %%%%%%%%%%%%
                 segment=3;
                %%%%%%%%%%%%                
                switch segment
                    case 1
                        time_highlight_path=[21+7/60+0/3600 21+9/60+0/3600];%
                    case 2
                        time_highlight_path=[21+11/60+0/3600 21+20/60+0/3600];%  
                    case 3
                        time_highlight_path=[];%      
                end
                
            case 3
                times=[21+0/60 21+20/60]/24;        
                
                xlims3D=[280 340];
                ylims3D=[310 340];
                
                
                %%%%%%%%%%%%
                 segment=4;
                %%%%%%%%%%%%                
                switch segment
                    case 1
                        time_highlight_path=[21+7/60+0/3600 21+9/60+0/3600];%
                    case 2
                        time_highlight_path=[21+11/60+0/3600 21+12/60+0/3600];%                        
                    case 3
                        time_highlight_path=[21+18/60+0/3600 21+20/60+0/3600];%                                                                        
                    case 4
                        time_highlight_path=[21+30/60+0/3600 21+45/60+0/3600];%                                                                                                

                end
                
                
             case 4
                times=[21+20/60 21+45/60]/24;        
                
                xlims3D=[280 340];
                ylims3D=[310 340];
                
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%                
                switch segment
                    case 1
                        time_highlight_path=[21+34/60+0/3600 21+42/60+0/3600];%                                                                                                

                end   
                
            case 5
                times=[21+30/60 21+45/60]/24;        
                
                xlims3D=[280 340];
                ylims3D=[310 340];
                
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%                
                switch segment
                    case 1
                        time_highlight_path=[21+34/60+0/3600 21+42/60+0/3600];%                                                                                                

                end       
                   
                        
                
                
            otherwise                                                                   
                                                    
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];% 
                time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%
                time_highlight_path=[16+26/60+10/3600 16+28/60+35/3600];%
                                
        end
        
        case '104'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=99;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%    
                    otherwise
                        time_highlight_path=[];
                        
                end
                
                xlims3D=[300 700];
                ylims3D=[200 400];
                
            case 2
                times=[14+0/60 14+20/60]/24;        
                
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+4/60+0/3600 14+12/60+0/3600];%
                    case 2
                        time_highlight_path=[14+13/60+0/3600 14+19/60+0/3600];%  

                end
                   
                        
                
                
            otherwise                                                                   
                                                    
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];% 
                time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%
                time_highlight_path=[16+26/60+10/3600 16+28/60+35/3600];%
                                
        end
        
     case '105'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[300 700];
                ylims3D=[200 400];
                
            case 2
                times=[14+0/60 14+20/60]/24;        
                
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+4/60+0/3600 14+12/60+0/3600];%
                    case 2
                        time_highlight_path=[14+13/60+0/3600 14+19/60+0/3600];%  

                end
                   
                        
                
                
            otherwise                                                                   
                                                    
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];% 
                time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%
                time_highlight_path=[16+26/60+10/3600 16+28/60+35/3600];%
                                
        end    
        
        
     case '108'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[300 700];
                ylims3D=[200 400];
                
   
                                
        end 
        
         case '113'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[0 700];
                ylims3D=[0 400];
                
   
                                
        end     
        
case '117'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[0 700];
                ylims3D=[0 400];
                
   
                                
        end           
        
        case '120'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[0 700];
                ylims3D=[0 400];
                
   
                                
        end
        
    case '122'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[0 700];
                ylims3D=[0 400];
                
   
                                
        end    
        
 case '123'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[0 700];
                ylims3D=[0 400];
                
   
                                
        end    
                
        
        
    case '100'

        itimser='antjan06_flt';
        
         xlims3D=[0 700];
         ylims3D=[0 700];
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=99;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1
                times=[13+35/60 14+10/60]/24;                
                times=[time_flt(1) 14+10/60]/24;
                
                segment=2;
                
                switch segment
                    case 1
                        time_highlight_path=[13+39/60+15/3600 13+45/60+37/3600];%needles on first straight run level at -0.7 degC
                    case 2
                        time_highlight_path=[13+52/60+46/3600 13+53/60+31/3600];%droplets near cloud top on ascent
                        time_highlight_path=[13+56/60+1/3600 13+57/60+40/3600];%droplets near cloud top on descent
                        time_highlight_path=[13+52/60+46/3600 13+58/60+40/3600];%both above two periods
                        
                end
                
            case 2
                times=[14+0/60 14+20/60]/24;        
                
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+4/60+0/3600 14+12/60+0/3600];%
                    case 2
                        time_highlight_path=[14+13/60+0/3600 14+19/60+0/3600];%  

                end
                
             case 3
                 times=[14+20/60 14+40/60]/24;         
                
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+25/60+0/3600 14+38/60+0/3600];% droplet only, ice only and mixed regions.
               
                    case 2
                        time_highlight_path=[14+38/60+0/3600 14+40/60+0/3600];% 

                end   
                
               case 4
                 times=[14+40/60 time_flt(end)]/24;  
                times=[14+40/60 15+2/60+20/3600]/24; 
                times=[14+40/60 15+2/60+20/3600]/24;        
                
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+25/60+0/3600 14+38/60+0/3600];% droplet only, ice only and mixed regions.
               
                    case 2
                        time_highlight_path=[14+38/60+0/3600 14+40/60+0/3600];% 

                end     
                
            case 5
                times=[15+2/60+20/3600 15+21/60+40/3600]/24;                  
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+25/60+0/3600 14+38/60+0/3600];% droplet only, ice only and mixed regions.
               
                    case 2
                        time_highlight_path=[14+38/60+0/3600 14+40/60+0/3600];% 

                end 
                
             
            case 6
               times=[15+21/60+40/3600 time_flt(end)]/24;               
                times=[15+21/60+40/3600 15+55/60+0/3600]/24;                              
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+25/60+0/3600 14+38/60+0/3600];% droplet only, ice only and mixed regions.
               
                    case 2
                        time_highlight_path=[14+38/60+0/3600 14+40/60+0/3600];% 

                end                                     
                        
                
                
            otherwise                                                                   
                                                    
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];% 
                time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%
                time_highlight_path=[16+26/60+10/3600 16+28/60+35/3600];%
                time_highlight_path=[];%
                
                times=[time_flt(1) time_flt(end)]/24;%                
                                
        end
        
        case '101'

        itimser='antjan06_flt';
        
         xlims3D=[0 700];
         ylims3D=[0 700];
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
            switch transect
                case 1


                    time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];% 
                    time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%
                    time_highlight_path=[16+26/60+10/3600 16+28/60+35/3600];%
                    time_highlight_path=[];%

                    times=[time_flt(1) time_flt(end)]/24;%                

            end
        
        
        case '99'            

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                %%%%%%%%%%%%
                 segment=1;
                %%%%%%%%%%%%
                
                switch segment
                    case 1
                        time_highlight_path=[19+0/60+0/3600 19+30/60+0/3600];%
                    case 2
                        time_highlight_path=[19+30/60+0/3600 20+30/60+0/3600];%
                    case 3
                        time_highlight_path=[20+30/60+0/3600 22+0/60+0/3600];%                        
                        
                end
                
                xlims3D=[300 700];
                ylims3D=[200 400];
                
                
%                xlims3D=[300 700];
%                ylims3D=[100 500];
                
            case 2
                times=[14+0/60 14+20/60]/24;        
                
                
                segment=2;
                switch segment
                    case 1
                        time_highlight_path=[14+4/60+0/3600 14+12/60+0/3600];%
                    case 2
                        time_highlight_path=[14+13/60+0/3600 14+19/60+0/3600];%  

                end
                   
                        
                
                
            otherwise                                                                   
                                                    
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];% 
                time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%
                time_highlight_path=[16+26/60+10/3600 16+28/60+35/3600];%
                                
        end
        

        

       end 
        
        
        



%load wrf data first for plotting the terrain on 3D graphs

% if ~exist('is_met_em')
%     is_met_em=0;
% end

        
% Make sure that set i_highlight_path=0 if aren't highlighting a section
% otherwise the plot might not be restricted to the section requested
% in the times array

for iplot_comparison=1:length(plot_graph_case)

    plot_graph_case2 = plot_graph_case{iplot_comparison};
    
    ismooth=[0 0 0 0 0 0];   
    i_set_dateticks=0;
    iadd_nums_above=0;
    
    switch plot_graph_case2
        
        case '3D flight track'
                

                
        itimser='antjan06_flt';                
        man_choose_water_graph=1;
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        iadd_terrain=1; %switch these off for some 3D views        
        iplot_latlon=0;
        i_highlight_path=1;
        iset_xlimits=1; 
        xlimits=times*24;
        flt_graph='Altitude3D';
%        flt_graph='Temperature3D';
        highlight_type='fill';
%        multisaveplot
        waterVapourMay2005        
        i3D=1; %switch to zoom in on a certain area - switch on for top down
        zoom_and_highlight_plot
        
        case 'Pressure'
        
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Pressure';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
       case 'Vertical speed of aircraft'
            man_choose_water_graph=1;
            itimser='antjan06_flt';
            graph=46;
            man_choose_flt_graph=1;
            man_choose_itimser=1;
            i_highlight_path=0;
            iset_xlimits=1;  %if want a close up then will have to set these
                              %and sort out the cloud highlight bit - seems to work through setting
                              %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)                    

            set_ix_distance=1;          
            ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
            xlimits=times*24;
            flt_graph='Vertical speed of aircraft';
            highlight_type='fill';
            waterVapourMay2005        
            zoom_and_highlight_plot
        
       case 'Temp'
            
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=0;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)                    
                          
        set_ix_distance=1;          
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Temp';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
      case 'Estimated displacement'
            
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=0;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)                    
                          
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Estimated displacement';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot  
        
    case 'Number aerosol 610-1030 nm'
         man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=0;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)                    
                          
        set_ix_distance=1;          
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Number aerosol 610-1030 nm';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot  
        
        
        
     case 'Vertical wind speed'
            
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Vertical wind speed';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        
        
   case 'Potemp'
        
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Potemp';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        
    case 'RHi hygrometer'
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='RHi hygrometer';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
    case 'RHi humi'        
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='RHi humi';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
 case 'RH hygrometer'
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='RH hygrometer';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
    case 'RH humi'        
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='RH humi';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot        
        
    case 'Altitude'        
         man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                        
        xlimits=times*24;
        flt_graph='Altitude';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
    case 'Aircraft T pert'
        
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Aircraft T pert';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
    
        case 'LWC'
                        
            
        man_choose_water_graph=1;
        graph=46;
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;       
        time_graph = {'LWC_dist_CAS','LWC_hotwire'};
%        time_graph = {'LWC_hotwire'};
%        time_graph = {'LWC_dist_CAS'};
%        time_graph = {'LWC_dist_CAS','LWC_dist_CAS'};
%        time_graph = {'LWC mode','LWC mode'};
 %       time_graph = {'LWC mode'};
%        instrument={'BAS CAS','MAN CAS'};
        instrument={'BAS CAS','BAS CAS'};
%        instrument={'BAS CAS'};
%        instrument={'MAN CAS'};        
        highlight_type='fill';    
        
        set_ix_distance=1;          
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track 
        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        case 'Airspeed'
                  
         man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Airspeed';
%        flt_graph='Airspeed CAS';        
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        
        case 'Total number CAS'
        
        man_choose_water_graph=1;
        graph=46;
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        ix_distance=1;
        time_graph = {'Total number CAS'};
        instrument={'BAS CAS'};
        highlight_type='fill';        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        case 'Total large ice number CIP'
        
        man_choose_water_graph=1;
        graph=46;
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Total large ice number CIP'};
        highlight_type='fill';        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        case 'Ice Number (Jonny)'
        
        instrument={'BAS CAS'};
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Ice number'};
        highlight_type='fill';   
        
        set_ix_distance=1;          
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        %set to one for distance (not 2) for this case (CAS plots)
        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        case 'Total number size range (Jonny)'
            instrument={'BAS CAS'};
            man_choose_water_graph=1;
            graph=46; %timeseries
            itimser='CAS plots';
            man_choose_itimser=1;
            man_choose_flt_graph=1;
            i_highlight_path=1;     
            time_graph = {'Total number size range (Jonny)'};
            highlight_type='fill';   

            set_ix_distance=1;          
            ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                

            waterVapourMay2005
            zoom_and_highlight_plot
            
        case 'Binned ice number size range (Jonny)'
            instrument={'BAS CAS'};
            man_choose_water_graph=1;
            graph=46; %timeseries
            itimser='CAS plots';
            man_choose_itimser=1;
            man_choose_flt_graph=1;
            i_highlight_path=1;     
            time_graph = {'Binned ice number size range (Jonny)'};
            highlight_type='fill';   

            set_ix_distance=1;          
            ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                

            waterVapourMay2005
            zoom_and_highlight_plot
            
        case 'N crystals per bin ice number size range (Jonny)'
            instrument={'BAS CAS'};
            man_choose_water_graph=1;
            graph=46; %timeseries
            itimser='CAS plots';
            man_choose_itimser=1;
            man_choose_flt_graph=1;
            i_highlight_path=1;     
            time_graph = {'N crystals per bin ice number size range (Jonny)'};
            highlight_type='fill';   

            set_ix_distance=1;          
            ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                

            waterVapourMay2005
            zoom_and_highlight_plot
            
        case 'Number of ice accepts (Jonny)'
            instrument={'BAS CAS'};
            man_choose_water_graph=1;
            graph=46; %timeseries
            itimser='CAS plots';
            man_choose_itimser=1;
            man_choose_flt_graph=1;
            i_highlight_path=1;     
            time_graph = {'Number of ice accepts (Jonny)'};
            highlight_type='fill';   

            set_ix_distance=1;          
            ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                

            waterVapourMay2005
            zoom_and_highlight_plot
            
     case 'Round Number (Jonny)'  %round particles as classified by Jonnny's s/ware
        
        instrument={'BAS CAS'};
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Round Number (Jonny)'};
        highlight_type='fill';   
        
        set_ix_distance=1;          
        ix_distance=2;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        
        case 'Ice Mass (Jonny)'
        
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Ice mass'};
        highlight_type='fill';   
        
        set_ix_distance=1;          
        ix_distance=1;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
                
        waterVapourMay2005
        zoom_and_highlight_plot
        
     case 'Ice Diameter (Jonny)'
        
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Mean ice diameter'};
        highlight_type='fill';   
        
        set_ix_distance=1;          
        ix_distance=1;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
                
        waterVapourMay2005
        zoom_and_highlight_plot
        
        
        
     case 'Wind direction'        
         man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Wind dir';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
      case 'Wind speed'        
         man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit - seems to work through setting
                          %inds in timeseries_dan. ydat is set and then changed after ydat(1).y(inds)
        set_ix_distance=1;          
        ix_distance=0;  %whether to plot x-axis as distance or time  %2 is for plotting as a function of distance along the flight track                                
        xlimits=times*24;
        flt_graph='Wind';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        
        
        case 'Wind direction profile'
 
         highlight_type='';
            ihighlight_cloud=0;
            man_choose_water_graph=1;
            graph=89;
            man_choose_case89_xlab=1;
            i_highlight_path=1;
            xlab = 'Wind direction (degrees)'; 
            waterVapourMay2005
            save_plot
            
        case 'Wind speed'
            
            man_choose_water_graph=1;
            graph=89;
            man_choose_case89_xlab=1;
            i_highlight_path=1;
            xlab = 'Wind speed (m s^{-1})';
            waterVapourMay2005
            save_plot
        
 case 'need to fill in'
        
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Ice mass'};
        highlight_type='fill';        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Mean ice diameter'};
        highlight_type='fill';        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        
                        
        man_choose_plotTimeHeight_graph=1;
        i_highlight_path=1;        
        i577 = 'Particle size dist vs time';
        highlight_type='rectangle';        
        plotTimeHeightVap3
        zoom_and_highlight_plot
        

        

        
        

           
    end
        
end


disp('done');