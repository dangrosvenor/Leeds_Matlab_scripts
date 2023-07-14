
ioverride_plotglobal_loc=1;
ioverride_plotglobal_thresh=1;  %comment out if want to use screenings set in plot_global
%time override should aready be set (ioverride_time_selection)
ioverride_years_time_screen=1; %required to specify the different years
inew_cticks=1;  %colorbar is of the non-linear type
plot_global_maps
if exist('isave') & isave==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',save_notes_filepath);
    %                    saveas_ps_fig_emf(gcf,savename,'',0);
end

%Store the monthly data for all the locations.
season_mean_ALL(:,:,iseason) = P_save;
season_Ndatap_ALL(:,:,iseason) = Npoints;
season_std_ALL(:,:,iseason) = P_std_dev;
