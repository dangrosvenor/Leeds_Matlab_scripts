%Defaults for the multi-panel plot function

fsize_legend=10;
fontsize=16;
resizing=1;

dX_xlab=0;
dY_xlab=0;

clear clims fig_tit idelete_cb idelete_legend pos_cb_set icbar_horiz pos_legend ...
  idelete_xlab idelete_ylab cbar_label ioverride_xticks ioverride_yticks icol_plot
    
for i=1:100
    clims{i}=[NaN NaN];
    fig_tit{i}='';
    idelete_cb{i}=0; %whether to delete the colorbars for each figure
    idelete_legend{i}=0; %same, but for the legend.
    %position of the colorbars in page relative coords
    pos_cb_set{i}=[0.1+0.1*(i-1) 0.1 0.3 0.08];  
    %colbar_loc{i}=''; %default - keep location as it came
        %Otherwise can get to 'eastoutisde' for a vertical one and
        %southoutwide for a horizontal one.  
        %However, setting colobar location doesn't seem to work....
    icbar_horiz{i}=0; %Use this instead - assume default is vertical. may need to deal with converting from horiz to vert at some point
    dX(i)=0; %amount to move each figure by in X and Y
    dY(i)=0;
    pos_legend{i}=[NaN NaN NaN NaN];
    idelete_xlab{i}=0;
    idelete_ylab{i}=0;
    idelete_xtick_labs{i}=0;
    cbar_label{i}='';
    ioverride_xticks{i}='';
    ioverride_yticks{i}='';
    icol_plot{i}=1; %whether is a color plot, or something else
end