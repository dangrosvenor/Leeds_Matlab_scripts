isave_highlight=0;
%do the plot in watervap first - then this just scales and saves it


plots_todo = 'flight segments';
%plots_todo = 'profiles';

savedir=['Y:\BAS_flights\flight' flight_no '\CAS_plots\']; 
savedir=['C:\Documents and Settings\G\My Documents\Ohio2010\flight' flight_no '\'];         

eval(['time_flt=time_flt' flight_no ';']);
        
switch flight_no
        case '102'            

        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1            
                times=[time_flt(1) time_flt(end)]/24;
                
                
                
                %%%%%%%%%%%%
                 segment=11;
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
                    case 8
                        time_highlight_path=[20+20/60+0/3600 20+40/60+0/3600];%
                    case 9
                        time_highlight_path=[20+40/60+0/3600 21+0/60+0/3600];%
                    case 10
                        time_highlight_path=[21+0/60+0/3600 21+20/60+0/3600];%

                    case 11
                        time_highlight_path=[21+30/60+0/3600 21+42/60+0/3600];%
    
                    case 99
                        time_highlight_path=[20 20.0001];%                                                                        


                        
                       
                        
                end
                
             xlims3D=[275 375];
             ylims3D=[300 400];   
             
             xlims3D=[280 340];
             ylims3D=[310 340];   
             
             xlims3D=[260 340];
             ylims3D=[310 340];  
                          
             xlims3D=[260 390];
             ylims3D=[310 400];  

             xlims3D=[280 330];
             ylims3D=[315 340];    
             
         case 2
                times=[20+20/60 20+40/60]/24;     
                times=[20+40/60 21+0/60]/24;   
%                times=[21+0/60 21+20/60]/24; 
%                times=[21+30/60 21+42/60]/24;                 
                
                xlims3D=[250 400];
                ylims3D=[275 425];
                                
                xlims3D=[275 375];
                ylims3D=[300 400];
                
                xlims3D=[280 330];
                 ylims3D=[315 340];
                
                %%%%%%%%%%%%
                 segment=99;
                %%%%%%%%%%%%                
                switch segment
                    case 1
                        time_highlight_path=[21+7/60+0/3600 21+9/60+0/3600];%
                    case 2
                        time_highlight_path=[21+11/60+0/3600 21+12/60+0/3600];%  
                    case 3
                        time_highlight_path=[21+18/60+0/3600 21+20/60+0/3600];%                          
					case 4
                        time_highlight_path=[21+30/60+0/3600 21+42/60+0/3600];%                        
                        
                    case 99
                        time_highlight_path=[];                                                                   
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
                times=[time_flt104(1) time_flt104(end)]/24;
                
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
        
    case '100'

        itimser='antjan06_flt';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        transect=6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        switch transect
            case 1
                times=[13+35/60 14+10/60]/24;                
                times=[time_flt100(1) 14+10/60]/24;
                
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
                 times=[14+40/60 time_flt100(end)]/24;  
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
               times=[15+21/60+40/3600 time_flt100(end)]/24;               
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
                                
        end
        
        

        

       end 
        
        
        
        

        
        
        switch plots_todo
            case 'flight segments'
                

        man_choose_water_graph=1;
        itimser='antjan06_flt';                         
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        iadd_terrain=1; %switch these off for some 3D views
        iplot_latlon=0;
        i_highlight_path=1;
        iset_xlimits=1; 
        xlimits=times*24;
        flt_graph='Altitude3D';
        highlight_type='fill';
%        multisaveplot
        waterVapourMay2005        
%        i3D=1; %switch to zoom in on a certain area - switch on for top down
        zoom_and_highlight_plot
        
        clear i3D
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit
        xlimits=times*24;
        flt_graph='Pressure';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        clear i3D        
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=0;
        iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit
        xlimits=times*24;
        flt_graph='Temp';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        clear i3D
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
         iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit
         xlimits=times*24;
        flt_graph='Potemp';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        
        clear i3D
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
         iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit
         xlimits=times*24;
        flt_graph='RHi hygrometer';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        clear i3D
        man_choose_water_graph=1;
        itimser='antjan06_flt';
        graph=46;
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
         iset_xlimits=1;  %if want a close up then will have to set these
                          %and sort out the cloud highlight bit
         xlimits=times*24;
        flt_graph='RH humi';
        highlight_type='fill';
        waterVapourMay2005        
        zoom_and_highlight_plot
        
        clear i3D
        man_choose_water_graph=1;
        itimser='antjan06_flt';        
        graph=46;        
        man_choose_flt_graph=1;
        man_choose_itimser=1;
        i_highlight_path=1;
        xlimits=times*24;        
         iset_xlimits=1;       
        flt_graph='Altitude';
        highlight_type='fill';        
        waterVapourMay2005
        zoom_and_highlight_plot
        
        clear i3D
        man_choose_water_graph=1;
        graph=46;
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;       
        time_graph = {'LWC_dist_CAS','LWC_hotwire'};
        time_graph = {'LWC_hotwire'};
        highlight_type='fill';        
        xlimits=times*24;        
         iset_xlimits=1;       
        waterVapourMay2005
        zoom_and_highlight_plot

        clear i3D        
        man_choose_water_graph=1;
        graph=46;
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Total number CAS'};
        highlight_type='fill';        
		xlimits=times*24;        
         iset_xlimits=1;                      
        waterVapourMay2005
        zoom_and_highlight_plot
        
        clear i3D        
        man_choose_water_graph=1;
        graph=46;
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Total large ice number CIP'};
        highlight_type='fill';        
		xlimits=times*24;        
         iset_xlimits=1;                       
        waterVapourMay2005
        zoom_and_highlight_plot
        
        clear i3D        
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Ice number'};
        highlight_type='fill'; 
		xlimits=times*24;        
         iset_xlimits=1;               
        waterVapourMay2005
        zoom_and_highlight_plot
        
        clear i3D       
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Ice mass'};
        highlight_type='fill'; 
		xlimits=times*24;        
         iset_xlimits=1;               
        waterVapourMay2005
        zoom_and_highlight_plot
        
        clear i3D
        man_choose_water_graph=1;
        graph=46; %timeseries
        itimser='CAS plots';
        man_choose_itimser=1;
        man_choose_flt_graph=1;
        i_highlight_path=1;     
        time_graph = {'Mean ice diameter'};
        highlight_type='fill'; 
		xlimits=times*24;        
         iset_xlimits=1;               
        waterVapourMay2005
        zoom_and_highlight_plot
        
        
        clear i3D                        
        man_choose_plotTimeHeight_graph=1;
        i_highlight_path=1;        
        i577 = 'Particle size dist vs time';
        highlight_type='rectangle';        
		xlimits=times;        
         iset_xlimits=1;                       
        plotTimeHeightVap3
        zoom_and_highlight_plot
        

        

        
        
        case 'profiles'
            highlight_type='';
            ihighlight_cloud=0;
            man_choose_water_graph=1;
            graph=89;
            man_choose_case89_xlab=1;
            i_highlight_path=1;
            xlab = 'Wind direction (degrees)'; 
            waterVapourMay2005
            save_plot
            
            man_choose_water_graph=1;
            graph=89;
            man_choose_case89_xlab=1;
            i_highlight_path=1;
            xlab = 'Wind speed (m s^{-1})';
            waterVapourMay2005
            save_plot
        end


