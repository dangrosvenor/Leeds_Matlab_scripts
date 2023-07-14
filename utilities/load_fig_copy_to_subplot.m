function [hsub_new,Hcb,pos_sub,pos_cb] = load_fig_copy_to_subplot(hf,fig_dir,fig_load,isub,fontsize)

hf_load = open([fig_dir fig_load]);
title(''); %remove the title
hs = subplot_spaceplots(3,2,isub,'parent',hf);
%hs = hA(1);
[hsub_new,Hcb,pos_sub,pos_cb]=fig2subplot(hf_load,hs);
temp = increase_font_size_map_figures_func(hsub_new,fontsize,0);
close(gcf);