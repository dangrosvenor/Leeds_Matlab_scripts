function plot_wind_quiver_arrows(u,v,x,y,nx,ny,scale_speed)

%    spx=25;^M
%    spy=15;^M
%    nx = 15;
%    ny = 15;
    
    idx=floor(length(x)/nx);
    idy=floor(length(y)/ny);
    
    xinds=[1:idx:length(y)];
    yinds=[1:idy:length(x)];
    
%    sf=maxALL(u)/maxALL(v);

if nargin==7 %if want to scale vectors - e.g. in order to make arrows proportional to wind speed on different plots
	speed = sqrt(u(xinds).^2 + v(yinds).^2);
	max_sp = maxALL(speed);
	sf = 2 * max_sp/scale_speed; %scale factor to make arrows proportional to the wind speed - sf=2 is about right for getting arrows the right size on the plot. Then scale according to a max wind speed expected (scale_speed m/s)
	quiver(x(xinds),y(yinds),u(xinds,yinds),v(xinds,yinds),sf,'w');

else
	quiver(x(xinds),y(yinds),u(xinds,yinds),v(xinds,yinds),'w');
end

    hold on;




