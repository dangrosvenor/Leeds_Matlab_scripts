graph=1;

scrsz=get(0,'ScreenSize');
%posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.13];
%posit=[9 50 scrsz(3)/1.9 scrsz(4)/2.1];  %used 22nd Jan, 2009
posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];

lwidth=3; %line width
dual=0;
xlims=0;
zlims=0; %flag to say that want y axis limited by zlimits(1) and zlimits(2)
ixtime=0; %flag set to one when x-axis is time so that times over 24hrs are changed
ititle=1; %flag to give graph a title of titlenam
xloc=0;
add_points=0; %flag to say whether to add labelled points on a plot
no_title=0; %switch to remove the title
iaxis_square=1; %switch to make axis square
%to run from 00:00

ixdir=0; %set to -1 to reverse x direction
iydir=0; %same for ydir

f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

secyA=0;
secyB=0;
lab2='';


if ~exist('subplotting'); subplotting=0; end
if ~exist('idir'); idir=1; end
if ~exist('file_type'); file_type=''; end

if ~exist('justplot')
    justplot=0;
end

clear labs xdat ydat diff xpos ypos point_labs
    figname='WRF profile';
    gridon=1; %switch for grid default =1 - note probs with grid when resizing due to extra ticks.s
    noplot=0;
    logflag=0;
    lor=1;  %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
    nmark=0; %no. of markers - set to -1 for markers for every point
    idirlabel=0; %flag to put directory of run in bottom left corner
    xlims=0;
    savename='';

    xlab='';
    ylab='';
    
    i_paper_labels=0;
    