function [pos_cb] = increase_font_size_map_figures_func(hax,ref_size,Hcb1)
%function [pos_cb] = increase_font_size_map_figures_func(hax,ref_size,Hcb1)
%Supply the colorbar handle, Hcb1, if known
%, otherwise will attempt to get the colorbar handle.
% Set Hcb1 to zero to leave the colorbar alone


%changes the font size for the labels around a map and moves the colorbar
%to compensate
%ref_size = 22; %for single panels
%ref_size = 16; %for 3x2 sub-plots

%change the fontsize of the colorbar (and title)
set(hax,'fontsize',ref_size);

%change the size of the font for the lat lon labels.
h = findobj(hax,'Type','text');
set(h,'Fontsize',ref_size-4);

if nargin<3
    %get the handle to the colorbar
    Hcb1 = find_peer_colorbars_of_an_axes(hax);


    if length(Hcb1)>0
        pos2 = get(Hcb1,'position'); %position of colorbar
        pos = get(hax,'position'); %position of main axis
        dY = pos(4)*0.15; %fraction of the size of the main y-axis to move colorbar down by
        pos_cb = [pos2(1) pos(2)-dY pos2(3) pos2(4)];
        
    else
        disp('No colorbar!');
        pos_cb=[];
    end

elseif Hcb1==0
 pos_cb = []; %do nothing
else
    
pos_cb = get(Hcb1,'position');

end

if Hcb1~=0
set(Hcb1,'position',pos_cb);
set(Hcb1,'xaxislocation','bottom');
end

fname = get(gcf,'FileName');
fname2 = [fname(1:end-4) '_LARGER_FONT'];
%saveas_ps_fig_emf(gcf,fname2,'',0,1);