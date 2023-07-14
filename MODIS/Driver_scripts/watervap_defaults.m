
% --- default flags for watervap
    subplotting = 0; clear idirs
    time_highlight_path=[];
    
    iytick_relabel_log=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
    y_axis_type=''; %default
    x_axis_type='';
    i_set_dateticks=0;
    iadd_nums_above=0;
    
    xlims=0;
    fsize=22;
    
    idatetick=0; %flag to put times as HH:MM instead of decimal time
    
    noplot=0;
    
    iwrite_text_dat=0; %
    
    ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
    line_pattern(1).p=NaN;
    line_colour(1).c=NaN;
    marker_style(1).m=NaN;
    line_widths(1).l = NaN;
    
    iovr_leg_line=0; %flag which plots a '-' line in the legend even if linestyle has been set to 'none'
    
    logflag=0;
    idate_ticks_fix=0;
    
    
    iaxis_square=1; %switch to make axis square                                                        
    
    
    
     iytick_relabel_log=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
    y_axis_type=''; %default
    x_axis_type='';
    i_set_dateticks=0;
    iadd_nums_above=0;
    
    xlims=0;
    fsize=22;
    
    idatetick=0; %flag to put times as HH:MM instead of decimal time
    
    noplot=0;
    
    iwrite_text_dat=0; %
    
    time_highlight_path=[];
    
   
    
    ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
    line_pattern(1).p=NaN;
    line_colour(1).c=NaN;
    marker_style(1).m=NaN;
    line_widths(1).l = NaN;
    
    iovr_leg_line=0; %flag which plots a '-' line in the legend even if linestyle has been set to 'none'
    
    
    for idat=1:99
        ismooth_x(idat)=0;
        ismooth_y(idat)=0;
    end
Nsmooth_window=5;
smooth_mode='mean';



add_ground_height=0.62;

lwidth=4; %line width
marksize=10;
%marksize=7;
dual=0;
izlim=0; %flag to say that want y axis limited by zmin and zmax
ixtime=0; %flag set to one when x-axis is time so that times over 24hrs are changed
ititle=1; %flag to give graph a title of titlenam
xloc=0;
add_points=0; %flag to say whether to add labelled points on a plot
no_title=0; %switch to remove the title
%to run from 00:00
ixtick_relabel_log=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
%x_axis_type=''; %default
iplot_3D=0;
ierror_bars='none';
errordatU=[];
errordatL=[];
ihighlight_points=0;
iset_xticks = 0; %set xticks to specified locations
iset_xticklabs=0; %set the xtick labels to specifed strings too
fsize_title = 14; %fontsize for the title



iexecute_script=0; %flag to say whether to execute a script at the end with name script_name
iadd_line=0; %add a line defined by addlineX and addlineY

ixdir=0; %set to -1 to reverse x direction
iydir=0; %same for ydir


secyA=0;
secyB=0;
lab2='';

zmax=25;
zmin=0;
%lor=0;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane

logflag=0;
lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=0; %no. of markers - set to -1 for markers for every point
idirlabel=0; %flag to put directory of run in bottom left corner
savename='';

gridon=1; %switch for grid default =1 - note probs with grid when resizing due to extra ticks.s


for idat=1:99
    ismooth_x_import(idat)=0;
    ismooth_y_import(idat)=0;
end

scrsz=get(0,'ScreenSize');
posit=[scrsz];
    
    
    
                                
                                
                                
% -- Examples of setting linestyles

% ichoose_styles=1; %flag to say whether we want to specifiy the specific line patterns and colours
% istyle=1;
% line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; line_widths(istyle).l = 3; istyle=istyle+1;
% 
% line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; line_widths(istyle).l = 2; istyle=istyle+1;
% line_pattern(istyle).p= '--'; line_colour(istyle).c=[1 0.7 0.7]; marker_style(istyle).m='d'; line_widths(istyle).l = 2; istyle=istyle+1;
% line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
% line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
% 
% line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='s'; line_widths(istyle).l = 3; istyle=istyle+1;
% 
% line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
% line_pattern(istyle).p= '--'; line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='o'; line_widths(istyle).l = 2; istyle=istyle+1;
% line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.7 0.7 0]; marker_style(istyle).m='s'; line_widths(istyle).l = 2; istyle=istyle+1;
% line_pattern(istyle).p= '--'; line_colour(istyle).c=[0.7 0.7 0]; marker_style(istyle).m='s'; line_widths(istyle).l = 2; istyle=istyle+1;


                                
% -- How to move one line to appear on top
%   uistack(h(1).h,'top');
                            

                            
                            
                            

