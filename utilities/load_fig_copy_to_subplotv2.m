function [hsub_new,Hcb,pos_sub,pos_cb,Hleg] = load_fig_copy_to_subplotv2(hf,fig_dir,fig_load,isub,fontsize,subN,subM,clims,colbar_loc)

filename=[fig_dir fig_load];

hf_load = open(filename);


if isnan(clims(1))==0
    caxis(clims);
end

if ~exist('colbar_loc')
    colbar_loc='';
end

 
 
%title(''); %remove the title
htit=get(gca,'title');
%set(htit,'string','');
hs = subplot_spaceplots(subN,subM,isub,'parent',hf);
%hs = hA(1);
[hsub_new,Hcb,pos_sub,pos_cb,Hleg]=fig2subplot(hf_load,hs,fontsize,colbar_loc);
 
 
temp = increase_font_size_map_figures_func(hsub_new,fontsize,0);
if length(Hcb)==0
    Hcb=NaN;  %to avoid an error when trying to pass the handle out of the function
end
close(gcf);