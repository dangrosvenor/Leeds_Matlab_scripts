function plot_wind_barbs(p,u,v)
% function plot_wind_barbs(p,u,v)
% plot wind barbs on the tephigram
k_gas=0.286;

alim=axis;
a=get(gca,'position');
axes('position',a);

% the rotated coordinates for iso pressure lines
p_iso=[1025 950 850 750 650 550 450 375 325 275 225 175 125 85 65].*100;
t_p_iso=180:1:350;t_p_iso=t_p_iso';
for i=1:length(p_iso)
    th_p_iso(:,i)=t_p_iso.*(100000./p_iso(i)).^k_gas;
    [t_p_iso_r(:,i),th_p_iso_r(:,i)]=rotate_coords(t_p_iso,th_p_iso(:,i),pi./4);
end
% rotated x we want to find
x_coord=435;
for i=1:length(p_iso)
    y_coord(i)=interp1(t_p_iso_r(:,i),th_p_iso_r(:,i),x_coord);
end
% interpolate to find gridded u, v
u_grid_iso=interp1(p,u,p_iso);
v_grid_iso=interp1(p,v,p_iso);
% plot wind barbs on graph
wind_barbs(x_coord.*ones(size(y_coord)),y_coord,u_grid_iso,v_grid_iso,2,'o');

axis equal;
set(gca,'ylim',[alim(3) alim(4)]);
%set(gca,'xlim',[alim(1) alim(2)]);
a=get(gca,'position');
set(gca,'position',[a(1)+0.35 a(2) a(3) a(4)])
set(gca,'ygrid','off');
set(gca,'xgrid','off');
set(gca,'ygrid','off');
set(gca,'color','none');
set(gca,'xtick',[]);
set(gca,'ytick',[]);box off
set(gca,'visible','off');
set(gca,'plotboxaspectratiomode','auto');
