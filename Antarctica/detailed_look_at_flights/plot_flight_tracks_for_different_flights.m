isave=0;

figure
add_terrain
axis square

clear savedir

eval(['dat_flt=dat_flt' flight_no ';']);




switch flight_no
    case '100'
        
        savedir='Y:\BAS_flights\flight100\CAS_plots\';

        manual_set=1;

        segment=1; %14:00 to 14:20

        
        ihighlight_sections=1; %flag to say whether to highlight some sections of flight path
        

        scatter_var='Time (UTC)';
        scatter_var='Altitude (m)';

        switch segment
            case 1
                time_inds=[time_flt100(1) 14+0/60];
                time_highlight_path=[13+40/60 13+45/60];%needles on first straight run level at -0.7 degC
                time_highlight_path=[13+52/60+46/3600 13+57/60+40/3600];%droplets near cloud top on ascent and descent
                xlims=[200 500];
                ylims=[250 550];


            case 2
                xlims=[200 500];
                ylims=[250 550];
                time_highlight_path=[13+52/60+46/3600 13+53/60+31/3600];%droplets near cloud top on ascent
                %                time_highlight_path=[13+56/60+1/3600 13+57/60+40/3600];%droplets near cloud top on descent
                time_highlight_path=[13+52/60+46/3600 13+57/60+40/3600];%droplets near cloud top on ascent
            case 3
                %        time_inds=[time_flt100(1) time_flt100(end)];
                time_inds=[14+0/60 14+20/60];
                %        time_inds=[14+20/60 14+37/60];

                time_highlight_path=[14+4/60 14+12/60];
                time_highlight_path=[14+13/60+0/3600 14+19/60+0/3600];%both above two periods
                %                xlims=[270 360];
                %                ylims=[370 460];

                xlims=[200 500];
                ylims=[250 550];

            case 4
                %        time_inds=[time_flt100(1) time_flt100(end)];
                time_inds=[14+20/60 15+0/60];
                time_inds=[14+20/60 14+40/60];
                %                time_inds=[14+40/60 15+0/60];
                %        time_inds=[14+20/60 14+37/60];

                %                time_highlight_path=[14+4/60 14+12/60];
                time_highlight_path=[14+25/60+0/3600 14+35/60+0/3600];%
                %               time_highlight_path=[14+38/60+0/3600 14+43/60+0/3600];%


                %                xlims=[270 360];
                %                ylims=[370 460];

                xlims=[200 500];
                ylims=[250 550];

            case 5
                %        time_inds=[time_flt100(1) time_flt100(end)];
                time_inds=[14+20/60 15+0/60];
                time_inds=[14+20/60 14+40/60];
                time_inds=[14+40/60 15+2/60+20/3600];
                %               time_inds=[15+2/60+20/3600 15+21/60+40/3600];

                %                time_highlight_path=[14+4/60 14+12/60];
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];%
                %                time_highlight_path=[15+0/60+0/3600 15+5/60+0/3600];%


                %                xlims=[270 360];
                %                ylims=[370 460];

                xlims=[200 500];
                ylims=[250 550];

            case 6

                time_inds=[15+2/60+20/3600 15+21/60+40/3600];

                %                time_highlight_path=[14+4/60 14+12/60];
                time_highlight_path=[14+40/60+0/3600 14+55/60+0/3600];%
                time_highlight_path=[15+3/60+0/3600 15+5/60+0/3600];%
                time_highlight_path=[15+16/60+30/3600 15+18/60+30/3600];%


                %                xlims=[270 360];
                %                ylims=[370 460];

                xlims=[200 500];
                ylims=[250 550];

            case 7

                time_inds=[15+21/60+40/3600 time_flt100(end)];

                time_highlight_path=[15+26/60+10/3600 15+28/60+35/3600];%

                xlims=[200 500];
                ylims=[250 550];


            case 8

                time_inds=[time_flt100(1) time_flt100(end)];

                time_highlight_path=[15+26/60+10/3600 15+28/60+35/3600];%

                xlims=[200 500];
                ylims=[250 550];

        end
        
        
    case '113'
        savedir='Y:\BAS_flights\flight113\CAS_plots\';
        
        ihighlight_sections=0; %flag to say whether to highlight some sections of flight path        

        manual_set=1;

        segment=1; %14:00 to 14:20
        segment=1; %
                
        scatter_var='Time (UTC)';
        scatter_var='Altitude (m)';

        switch segment
            case 1
                
                time_inds=[time_flt113(1) time_flt113(end)];
%                time_inds=[time_flt113(1) 15.5];

                xlims=[0 850];
                ylims=[0 850];
        end
        
     case '117'
        savedir='Y:\BAS_flights\flight117\CAS_plots\';
        
        ihighlight_sections=0; %flag to say whether to highlight some sections of flight path        

        manual_set=1;

        segment=1; %14:00 to 14:20
        segment=1; %
                
        scatter_var='Time (UTC)';
        scatter_var='Altitude (m)';

        switch segment
            case 1
                
                time_inds=[time_flt117(1) time_flt117(end)];
%                time_inds=[time_flt113(1) 15.5];

                xlims=[300 550];
                ylims=[300 550];
        end  
        
    case '104'
         savedir=['Y:\BAS_flights\flight' flight_no '\CAS_plots\'];
        
        ihighlight_sections=0; %flag to say whether to highlight some sections of flight path        

        manual_set=1;

        segment=1; %14:00 to 14:20
        segment=1; %
                
        scatter_var='Time (UTC)';
        scatter_var='Altitude (m)';

        switch segment
            case 1
                
                time_inds=[time_flt(1) time_flt(end)];
%                time_inds=[time_flt113(1) 15.5];

                xlims=[300 650];
                ylims=[200 550];
        end  
        
        
     case '102'
         savedir=['Y:\BAS_flights\flight' flight_no '\CAS_plots\'];
        
        ihighlight_sections=0; %flag to say whether to highlight some sections of flight path        

        manual_set=1;

        segment=2; %
                
        scatter_var='Time (UTC)';
%        scatter_var='Altitude (m)';

        switch segment
            case 1
                
                time_inds=[time_flt(1) time_flt(end)];
%                time_inds=[time_flt113(1) 15.5];
%                time_inds=[20+20/60+0/3600 21+00/60+0/3600];%
%                time_inds=[21+7/60+0/3600 21+09/60+0/3600];%   
%                time_inds=[21+11/60+0/3600 21+20/60+0/3600];%                   

                xlims=[200 550];
                ylims=[200 550];
                
                
             
                
            case 2
                    
                time_inds=[20+20/60 21+45/60];
                
                ihighlight_sections=0;
                time_highlight_path=[21+32/60+0/3600 21+42/60+0/3600];%
                
                xlims=[200 550];
                ylims=[200 550];
                
                
             xlims=[280 340];
             ylims=[310 340];
        end  

        
        case '99'
         savedir=['Y:\BAS_flights\flight' flight_no '\CAS_plots\'];
        
        ihighlight_sections=0; %flag to say whether to highlight some sections of flight path        

        manual_set=1;

        segment=1; %14:00 to 14:20
        segment=1; %
                
        scatter_var='Time (UTC)';
        scatter_var='Altitude (m)';

        switch segment
            case 1
                
                time_inds=[time_flt(1) time_flt(end)];
%                time_inds=[time_flt113(1) 15.5];

                xlims=[300 650];
                ylims=[200 550];
        end  


end


iplot_return=0;


iset_clims_scatter=1; %so that the colour limits are the time limits and not the terrain


flight_plot_type='MODIS cloud temperature';
flight_plot_type='Flight track segment';

switch flight_plot_type
    case 'MODIS cloud temperature'
        type_of_track = 'line'; %normal line plot
        ihighlight_sections=0; %flag to say whether to highlight some sections of flight path
        plot_flight_path

        MOD_DAT_STR='CT-273.15';
        %                MOD_DAT_STR='CTP';
        %                MOD_DAT_STR='CF';
        %                MOD_DAT_STR='CTP_IR';
        %                MOD_DAT_STR='CEFR';

        eval(['MOD_DAT=' MOD_DAT_STR ';']);

        pcolor(XLAT,YLAT,MOD_DAT);



        shading flat;

        clims=[minALL(MOD_DAT) maxALL(MOD_DAT)];
        %                clims=[350 800];
        set(gca,'clim',clims);

        lb_map=lbmap(256,'brownblue'); %nice colormap for colorblind people
        lb_map=lbmap(256,'bluegray'); %nice colormap for colorblind people  blue,bluegray, redblue
        lb_map=lbmap(256,'redblue'); %nice colormap for colorblind people  blue,bluegray, redblue

        cmap=flipud(lb_map);
        %                cmap=lb_map(1:100,:);

        %                colormap('copper');
        colormap(cmap);
    case 'Flight track segment'
        type_of_track = 'scatter_altitude';
        plot_flight_path
end

colorbar
axis([xlims ylims]);





xinds=findheight(x_grid,xlims);
xinds(1)=max([1 xinds(1)-1]);
xinds(end)=min([length(x_grid) xinds(end)+1]);
xinds=[xinds(1):xinds(2)];

yinds=findheight(y_grid,ylims);
yinds(1)=max([1 yinds(1)-1]);
yinds(end)=min([length(y_grid) yinds(end)+1]);
yinds=[yinds(1):yinds(2)];

hold on
plot_latlon_lines2

xlabel('X (km)');
ylabel('Y (km)');

switch flight_plot_type
    case 'MODIS cloud temperature'
        titlenam=['MODIS plot of ' MOD_DAT_STR];
        savename=[savedir titlenam];
    case 'Flight track segment'
        titlenam=[scatter_var ' from ' num2str(time_inds(1)) ' to ' num2str(time_inds(end)) ' UTC'];
        savename = [savedir 'Flight track in ' titlenam];
end

title(titlenam);




if isave==1
    print(gcf,[savename '.emf'],'-dmeta');
    saveas(gcf,savename,'fig');
end