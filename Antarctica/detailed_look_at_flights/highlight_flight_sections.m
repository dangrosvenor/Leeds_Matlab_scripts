if ihighlight_sections==1
    eval_str=['inds = findheight(time_flt' flight_no ',time_highlight_path);';];
    eval(eval_str);
    inds=inds(1):inds(2);
    
    eval_str=['LAT_plot = dat_flt' flight_no '(inds,2);'];
    eval(eval_str);
    eval_str=['LON_plot = dat_flt' flight_no '(inds,3);'];
    eval(eval_str);

    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_plot,LON_plot,0.1);
    
    X=(ilon-1)*DX;
    Y=(ilat-1)*DY;
    
    inds=round(1:length(X)/5:length(X));
    dy=diff(Y(inds));
    dx=diff(X(inds));
    angle = atan(dx./dy);
    angle = [angle(1) angle];
    
    d=25;
    xdiffs=d*sin(angle);
    ydiffs=d*cos(angle);
    
    xvals_left=X(inds)-xdiffs;
    yvals_left=Y(inds)+ydiffs;    

    xvals_right=X(inds)+xdiffs;
    yvals_right=Y(inds)-ydiffs;
    
    XX=[xvals_right fliplr(xvals_left)];
    YY=[yvals_right fliplr(yvals_left)];
    
%    fill(XX,YY,'g','linestyle','none');
    


%    scatter((ilon-1)*DX,(ilat-1)*DY,length(inds)*4,'g','o','filled');
    
    plot((ilon-1)*DX,(ilat-1)*DY,'wo','markerfacecolor','w','markersize',3);
end