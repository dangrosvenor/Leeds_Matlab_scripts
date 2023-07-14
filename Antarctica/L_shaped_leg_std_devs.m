clear std_leg mean_leg xdat ydat

var='wind speed';
%var='wind direction';
%var='Temp';

switch var
        case 'wind speed'
            col=col_wind;
        case 'wind direction';
            col=col_winddir;
        case 'Temp';
            col=col_temp;    
end

dat=dat_flt(:,col);
    
idat=0;
idat=idat+1;
[it it2]=findheight(time_flt19,20.4,20.57);
xdat(idat).x = dat_flt(it:it2,col_lat);
ydat(idat).y = dat(it:it2);
labs(idat).l = '15 m';
iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
aws_times = time_flt19(it:it2);
aws_time(idat) = aws_times(iaws_time(idat));

idat=idat+1;
[it it2]=findheight(time_flt19,20.97,21.17);
xdat(idat).x = dat_flt(it:it2,col_lat);
ydat(idat).y = dat(it:it2);
labs(idat).l = '152 m';
iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
aws_times = time_flt19(it:it2);
aws_time(idat) = aws_times(iaws_time(idat));


idat=idat+1;
[it it2]=findheight(time_flt19,21.21,21.38);
xdat(idat).x = dat_flt(it:it2,col_lat);
ydat(idat).y = dat(it:it2);
labs(idat).l = '305 m';
iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
aws_times = time_flt19(it:it2);
aws_time(idat) = aws_times(iaws_time(idat));


idat=idat+1;
[it it2]=findheight(time_flt19,21.77,21.96);
xdat(idat).x = dat_flt(it:it2,col_lat);
ydat(idat).y = dat(it:it2);
labs(idat).l = '610 m';
iaws_time(idat) = findheight_nearest(xdat(idat).x,-67.01);
aws_times = time_flt19(it:it2);
aws_time(idat) = aws_times(iaws_time(idat));


xlimits=[xdat(1).x(1) xdat(1).x(end)]; %whole range of the flight
%             izlim=1;
%             zmin=0;
%             zmax=20;
xlab=['Latitude'];

ylab = [ylab ' leg N-S'];



for i=1:length(xdat)
    %    istart=find(xdat(i).x>-67.15); %stop at -66.8 degrees south
    %    iend=find(xdat(i).x<-66.8); %stop at -66.8 degrees south

    inds = find(xdat(i).x>-67.15 & xdat(i).x<-66.8);

    
    switch var
        case {'wind speed','Temp'}
            std_leg(i) = std(ydat(i).y(inds));
            mean_leg(i) = mean(ydat(i).y(inds));
        case 'wind direction'; %don't do the mean, just the range for wind direction
            ydat(i).y=ydat(i).y+180; %need to add 180 to get conventional wind direction
            ydat(i).y(ydat(i).y>360)=ydat(i).y(ydat(i).y>360)-360;

            %for the direction will plot the range rather than the average
            min_leg = min(ydat(i).y);
            max_leg = max(ydat(i).y);
            mean_leg(i) = 0.5*(min_leg+max_leg); %the mid-point
            std_leg(i) = max_leg - mean_leg(i);
            
    end
    
    aws_value(i) = ydat(i).y(iaws_time(i));




end

iplot_AWS=1;

ileg=1; %the leg to plot. 15m is the first


if iplot_AWS==1
    
    plot(6+aws_time(ileg)/24,mean_leg(ileg),'bx','markersize',15);
    plot(6+aws_time(ileg)/24,mean_leg(ileg)-std_leg(ileg),'bs','markersize',10,'markerfacecolor','b');
    plot(6+aws_time(ileg)/24,mean_leg(ileg)+std_leg(ileg),'bs','markersize',10,'markerfacecolor','b');
    
%     switch var
%     case 'wind direction';
            plot(6+aws_time(ileg)/24,aws_value(ileg),'bo','markersize',10);
%     end
    
end

    %find the approx location of the mini descent min height (roughly 120m either side)
    [it,it2]=findheight_nearest(time_flt19,21.995,22.01);
    dat_mini = dat_flt19(it:it2,col);
    heights = dat_flt19(it:it2,col_alt)-65;
    times = dat_flt19(it:it2,1)/1e3/3600; %hours
    iheights = find(heights<=25);
    [min_height i]=min(heights);

    switch var
        case {'wind speed','Temp'}
            mean_dip = mean(dat_mini(iheights));
            std_dip = std(dat_mini(iheights));
            time_dip = times(i);
        case 'wind direction';
            dat_mini=dat_mini+180; %need to add 180 to get conventional wind direction
            dat_mini(dat_mini>360)=dat_mini(dat_mini>360)-360;

            %for the direction will plot the range rather than the average
            min_dip = min(dat_mini);
            max_dip = max(dat_mini);
            mean_dip = 0.5*(min_dip+max_dip);
            std_dip = max_dip - mean_dip;
    end



if iplot_AWS==1    
    plot(6+time_dip/24,mean_dip,'bx','markersize',15);
    plot(6+time_dip/24,mean_dip-std_dip,'bs','markersize',10,'markerfacecolor','b');
    plot(6+time_dip/24,mean_dip+std_dip,'bs','markersize',10,'markerfacecolor','b');
end
