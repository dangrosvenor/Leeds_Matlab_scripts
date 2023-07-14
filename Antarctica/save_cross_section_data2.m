isave=1; %saves the figure even if isave==0

savedir_plots = '/home/disk/eos1/d.grosvenor/Antarctica/plots/';

XY_pot_cross_data.X_cross = timesTH(1).t(1:pend);
XY_pot_cross_data.Y_cross = zz(1).z;

%savedir = 'Y:\WRF\ecmwf_ml_0.5_nudging\figures_for_paper_Aug2010\';
savedir = '/home/disk/eos1/d.grosvenor/matlab/work/Antarctica/Fohn_upstream_work/';
%savename = [savedir 'pot_slice_streamline_LAT=67.5_z0=950_str_succ=129_0-5km.mat'];

savename2 = remove_character(short_plot_name,' ','_'); savename2 = remove_character(savename2,':','');
savename = [savedir savename2 '_' datestr(now,30) '.mat'];

savename_plot = [savedir_plots savename2 '_' datestr(now,30)];
saveas_ps_fig_emf(gcf,[savename_plot],'',0,1);


switch var_plot
    case 'Wind speed (m s^{-1})'       
        data_type='wind';
        U_cross_AP = pdat(1).p(:,1:pend);
        if isave==1
            save(savename,'U_cross_AP');
        end
    case 'Potential temperature (K)'
        data_type='potemp';
        pot_cross_AP = pdat(1).p(:,1:pend);
        if isave==1
            save(savename,'pot_cross_AP');
        end
    case 'Vertical wind speed (m s^{-1})'
        data_type='vert wind';
         vertwind_AP = pdat(1).p(:,1:pend);
         if isave==1
            save(savename,'vertwind_AP');
        end
    case 'Component horizontal (UV) wind speed (m s^{-1})'
        data_type = 'wind_component';
        Ucomp_cross_AP = pdat(1).p(:,1:pend);
        if isave==1
            save(savename,'Ucomp_cross_AP');
        end
    case 'Wind direction (degrees)'
        data_type = 'wind_dir';
        WindDir_cross_AP = pdat(1).p(:,1:pend);
        if isave==1
            save(savename,'WindDir_cross_AP');
        end
        
    case 'Relative humidity (%)'   
        data_type = 'RH';
        RH_cross_AP = pdat(1).p(:,1:pend);
        if isave==1
            save(savename,'RH_cross_AP');
        end
end

if isave==1
        save(savename,'XY_pot_cross_data','lon_slice','d_line','x_line','y_line','hor_vert','-APPEND');
    %x_line and y_line describe the straight line used to specify the cross
    %section (unless were doing cross sections along streamlines).
    %hor_vert tells us this - 0 = straight line, 3 = along streamline
end



% if isave==1
% 
% %     if exist(savename)==2
% %         app_str = ',''-APPEND''';
% %     else
% %         app_str = '';
% %     end
% 
%     switch data_type
%         case 'wind'
% %            eval_str = ['save(savename,''XY_pot_cross_data'',''U_cross_AP'',''lon_slice'',''d_line''' app_str ');'];           
%             save(savename,'U_cross_AP');
%         case 'wind_component'
%             save(savename,'U_cross_AP');
%         case 'potemp'
% %            eval_str = ['save(savename,''XY_pot_cross_data'',''pot_cross_AP'',''lon_slice'',''d_line'' ');'];
%             save(savename,'pot_cross_AP');
%         case 'vert wind'
% %            eval_str = ['save(savename,''XY_pot_cross_data'',''vertwind_cross_AP'',''lon_slice'',''d_line'' ');'];
%             save(savename,'vertwind_cross_AP');
% 
%     end
% 
% %    eval(eval_str);
%     
%     save(savename,'XY_pot_cross_data','lon_slice','d_line','x_line','y_line','hor_vert','-APPEND');
%     %x_line and y_line describe the straight line used to specify the cross
%     %section (unless were doing cross sections along streamlines).
%     %hor_vert tells us this - 0 = straight line, 3 = along streamline
%     
%     
% end