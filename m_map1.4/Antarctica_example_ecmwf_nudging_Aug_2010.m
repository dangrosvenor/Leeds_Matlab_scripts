%set the projection type to polar sterographic
%only seems to work well when set latitude to be the poles
%lon=98.0 makes the grid box square for the ecmwf_ml_nudgding runs
%since this is what stand_lon is set to in namelist.wps
m_proj('stereographic','lat',-90,'lon',-98.0,'rad',67);
%xaxislocation property is important for south pole plots as otherwise
%it puts the longitude labels in the middle, which is messy
%setting to 'top' means they are at the northernmost part of the plot
figure
%draw the coastlines - patch does it as a filled patch instead of a line
m_coast('patch',[.6 .6 .6]);
set(gca,'fontsize',24);

m_grid('xaxislocation','top','yaxislocation','middle','tickdir','out','linest','-');
%N.B. m_ungrid removes a grid so that another can be drawn instead


imanual_select_load_case=1;
rundir(1).dir='ecmwf_ml_0.5_nudging'; fileWRF(1).file=['d01'];   %usual ECMWF run (24th April, 2009) 
load_WRF_vars

%now draw the box of the WRF domain using the line plotting equivalent
%for the mapping tool - specify all the pairs of lat lon for each value
%on the box boundary - actuallly, is a square box anyway (is plotting 
%correctly - can test by only plottin half the of boudary values for example
%Must be because was using a sterographic projection in WRF too
m_line(lon2d.var(1,:),lat2d.var(1,:),'color','k','linewidth',2);
m_line(lon2d.var(end,:),lat2d.var(end,:),'color','k','linewidth',2);
m_line(lon2d.var(:,1),lat2d.var(:,1),'color','k','linewidth',2);
m_line(lon2d.var(:,end),lat2d.var(:,end),'color','k','linewidth',2);

H=m_text(lon2d.var(end,1),lat2d.var(end,1)+1.8,'Domain 1');
set(H,'FontSize',16);


imanual_select_load_case=1;
rundir(1).dir='ecmwf_ml_0.5_nudging'; fileWRF(1).file=['d02'];   %usual ECMWF run (24th April, 2009) 
load_WRF_vars

%now draw the box of the WRF domain using the line plotting equivalent
%for the mapping tool - specify all the pairs of lat lon for each value
%on the box boundary - actuallly, is a square box anyway (is plotting 
%correctly - can test by only plottin half the of boudary values for example
%Must be because was using a sterographic projection in WRF too
m_line(lon2d.var(1,:),lat2d.var(1,:),'color','k','linewidth',2);
m_line(lon2d.var(end,:),lat2d.var(end,:),'color','k','linewidth',2);
m_line(lon2d.var(:,1),lat2d.var(:,1),'color','k','linewidth',2);
m_line(lon2d.var(:,end),lat2d.var(:,end),'color','k','linewidth',2);

H=m_text(lon2d.var(end,1)-0.5,lat2d.var(end,1)+1.8,'Domain 2');
set(H,'FontSize',16);



imanual_select_load_case=1;
rundir(1).dir='ecmwf_ml_0.5_nudging'; fileWRF(1).file=['d03'];   %usual ECMWF run (24th April, 2009) 
load_WRF_vars

%now draw the box of the WRF domain using the line plotting equivalent
%for the mapping tool - specify all the pairs of lat lon for each value
%on the box boundary - actuallly, is a square box anyway (is plotting 
%correctly - can test by only plottin half the of boudary values for example
%Must be because was using a sterographic projection in WRF too
m_line(lon2d.var(1,:),lat2d.var(1,:),'color','k','linewidth',2);
m_line(lon2d.var(end,:),lat2d.var(end,:),'color','k','linewidth',2);
m_line(lon2d.var(:,1),lat2d.var(:,1),'color','k','linewidth',2);
m_line(lon2d.var(:,end),lat2d.var(:,end),'color','k','linewidth',2);

H=m_text(lon2d.var(end,1)-2,lat2d.var(end,1)+1.8,'Domain 3');
set(H,'FontSize',16);