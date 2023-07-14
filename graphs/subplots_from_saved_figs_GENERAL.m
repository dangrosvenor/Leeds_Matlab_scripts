figlab='blah blah';
fig_dir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS/plots_20190903T071109/';
Nsub=1;
Msub=2;

clear fig_load pos_sub pos_cb
ifig=1;

for i=1:12
    clims{i}=[NaN NaN];
end

%Create brown/blue colormap - better for pltos showing negative and
%positive values (centred around zero).
 lb_map=lbmap(16,'brownblue'); %nice colormap for colorblind people       
 i_reverse_cmap=1;
 if i_reverse_cmap==1
     cmap=flipud(lb_map);
 else
     cmap=lb_map;
 end
 
colormap(cmap);

%Couldn't get the colorbar to display properly, so leave this for now. Will
%use images with the colormap already changed.


% ---- Fig. a ------


ifig=1;


% ------------ NEW blue to brown colourscale -------------------------

clims_all = [0 100];
%CAM5 - Day and Night
clims{ifig} = clims_all;
fig_load{ifig} = 'SW_TOA_model_time-mean_calc.fig'; ifig=ifig+1;
clims{ifig} = clims_all;
fig_load{ifig} = 'SW_TOA_calculated_from_obs_mean_values_using_average_CF'; ifig=ifig+1;



ichange_colormap = 1; %flag to overrule the default jet colormap with those specified above in cmaps



run_multi_file ='run_multi_new_GENERAL';
eval(run_multi_file);
%finally, mnake a copy of the above file with the tag to match the plots
%(from saveas...). This will allow the figure to re-created in case
%run_multi... changes at any point. E.g. for different layouts. Will need
%to copy the logged file over the original run_multi* file (after backing it up) 
%eval_str = ['!C:\cygwin\bin\cp "C:\Users\Dan\Documents\MATLAB\work\graphs\' run_multi_file '.m" "C:\Users\Dan\Documents\MATLAB\work\graphs\run_multi_panel_copy_commands_log\' run_multi_file '_' datestr_now '.m"'];
%eval(eval_str);