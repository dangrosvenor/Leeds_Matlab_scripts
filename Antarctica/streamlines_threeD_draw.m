figdir='C:\Documents and Settings\dan\My Documents\WRF\plots\plots_for_AP_paper\level4_winds_6UTC_6thJan.fig';
%note - the above figure is wrong as the x and y grids were slightly offset
figdir='C:\Documents and Settings\dan\My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\1_Wind speed (m s^{-1}) at level 4 (~293m above terrain) for 06 Jan 06,00 UTC for ecmwf-ml-0.5-nudging_CORRECT_grid_x_grid_y.fig';
figdir='C:\Documents and Settings\dan\My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\1_Terrain height 06 Jan 06,00 (m) for ecmwf-ml-0.5-nudging_CORRECT_grid_x_grid_y.fig';
figdir='C:\Documents and Settings\dan\My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\1_Terrain height 06 Jan 06,00 (m) for ecmwf-ml-0.5-nudging_pcolor_CORRECT_grid_x_grid_y.fig';
figdir='C:\Documents and Settings\dan\My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC_black_coarse_2km_xy_closeup\1_First level height 06 Jan 06,00 (m) for ecmwf-ml-0.5-nudging_pcolor_CORRECT_grid_x_grid_y.fig';
%WARNING - the above one plots the terrain as the lowest WRF model level (as defined by the -9e9 values in ufine and vfine etc)
figdir='C:\Documents and Settings\dan\My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC_black_coarse_2km_xy_closeup\Terrain height 06 Jan 06,00 (m) for ecmwf-ml-0.5-nudging_pcolor_CORRECT_grid_x_grid_y_ZOOM.fig';
%the above is the full field (not zoomed out)

open(figdir);
%this a good figure to plot the streamlines on


iunplot=0;
if iunplot==1
    for i=1:length(STREAM_z0)
        unplot
    end
end

x_line=[235.5 602.9];
y_line=[413 254.8]; %

dX=diff(x_line);
dY=diff(y_line);

D=sqrt(dX^2+dY^2);

d_perp = 250; %length of the perpendicular line to draw
%d_perp = 300; %length of the perpendicular line to draw
d_perp = 400; %length of the perpendicular line to draw
nstream=200;
%nstream=800;

px1=x_line(1) + (-dY*d_perp/D/2); %perpendiulcar vector drawn as (x_line(1),y_line(1)) + f*(-dY,dX) and
% (x_line(1),y_line(1)) - f*(-dY,dX).
% (-dY,dX) is the vector perpendicular to (dX,dY) and f is just the fraction of the lenght of it to draw
px2=x_line(1) - (-dY*d_perp/D/2);


dpx=px2-px1;
px=[px1:dpx/nstream:px2];

px=px+55;

py1=y_line(1) + (dX*d_perp/D/2);
py2=y_line(1) - (dX*d_perp/D/2);

dpy=py2-py1;
py=[py1:dpy/nstream:py2];

%%%% set starting height here %%%%%%
if ~exist('iset_zstart')
    zstart=950; %height of streamline starting positions in m
else
    clear iset_zstart %if are setting zstart from elsewhere
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




pz=zstart*ones(size(px));
%STREAM_z0=stream3(xstream,ystream,zstream,ufine,vfine,wfine,px*1e3,py*1e3,pz);
%STREAM_z0=stream3(xstream,ystream,zstream,ufine,vfine,wfine,px*1e3,py*1e3,pz,1); %the ,1 at the end selects a stepsize of 1*gridcell
%the default is 0.1 of a cell.
STREAM_z0=stream3(xstream,ystream,zstream,ufine2,vfine2,wfine2,px*1e3,py*1e3,pz,1); %the ,1 at the end selects a stepsize of 1*gridcell

%using a coarser z grid
%STREAM_z0=stream3(xstream3,ystream3,zstream3,ufine3,vfine3,wfine3,px*1e3,py*1e3,pz,1); %the ,1 at the end selects a stepsize of 1*gridcell

terr=nc{'HGT'}(1,:);

col='w';
for i=1:length(STREAM_z0)    
            plot(STREAM_z0{i}(:,1)/1e3,STREAM_z0{i}(:,2)/1e3,col,'linewidth',0.5);
            %diff(STREAM_z0{152}(:,3))./diff(STREAM_z0{152}(:,1));
            
%             dx=diff(STREAM_z0{i}(:,1));
%             istr_end=find(abs(dx)<5);
            
%if have used ufine2 (with NaNs) then can just do this:
            inan=isnan(STREAM_z0{i}(:,1));
            istr_end=find(inan==1)-1;
            
            ndiff=3;
            if length(istr_end)>=1
                if istr_end(1)>ndiff+1
                    %now find the terrain gradient by 2D interpolation of the last two points
                    
                    terrA = interp2(x_grid*1e3,y_grid*1e3,squeeze(terr_fine),STREAM_z0{i}(istr_end(1)-ndiff,1),STREAM_z0{i}(istr_end(1)-ndiff,2));
                    terrB = interp2(x_grid*1e3,y_grid*1e3,squeeze(terr_fine),STREAM_z0{i}(istr_end(1),1),STREAM_z0{i}(istr_end(1),2));

                    if terrA>terrB
                        plot(STREAM_z0{i}(:,1)/1e3,STREAM_z0{i}(:,2)/1e3,'k','linewidth',2);
                    end
                    % %    x=[x; STREAM_z0{i}{1}(:,1)];
                    % %    y=[y; STREAM_z0{i}{1}(:,2)];
                end
            end
end

iset_lims=1;
if iset_lims==1
%    set(gca,'clim',[1200 1600]);
    set(gca,'clim',[1000 2000]);    
    colorbar
%    set(gca,'xlim',[330 490]);
%    set(gca,'ylim',[230 390]);
%    set(gca,'xlim',[250 500]);
%    set(gca,'ylim',[180 430]);
    
%    set(gca,'xlim',[150 500]);
%    set(gca,'ylim',[130 480]);    
end

z0_string=['z_{s0} = ' num2str(zstart) ' m'];
xlims=get(gca,'xlim');
ylims=get(gca,'ylim');
text(xlims(2)+0.05*(xlims(2)-xlims(1)),ylims(2)+0.05*(ylims(2)-ylims(1)),z0_string);
savename = ['Streamline plot for 6UTC 6th Jan for ' z0_string];

iadd_cross_section_streamline=0;
if iadd_cross_section_streamline==1
    plot(STREAM2{1}(:,1),STREAM2{1}(:,2),'r','linewidth',2);
end



disp('Finished plotting streamlines');