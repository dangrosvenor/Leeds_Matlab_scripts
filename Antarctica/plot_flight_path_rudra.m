just_track=0;
just_labels=1; %flag to say whether to just print the location labels (and not the flight track).
i_use_set_labels=0; %use the location labels from watervap...
i_add_cross_section_line=0;

i_plot_location_labels=1;

marker_type='filled square';
marker_type='plus sign';

type_of_track = 'line'; %normal line plot
%type_of_track = 'scatter_altitude';


clear text %in case there has been a variable named 'text' as this causes trouble when using text command
np = 100;
%ntot = length(dat_flt19(:,2));
ntot = findheight(time_flt19,20.8);
inds = round([1:ntot/np:ntot]);

%LAT_plot = dat_flt19(inds,2);
%LON_plot = dat_flt19(inds,3);

% [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_plot,LON_plot,0.1);
 
 hold on
%  x_vals_flight = (ilon-1)*DX;
%  y_vals_flight = (ilat-1)*DY;
%  inds_flight = inds;

if just_labels==0
    switch type_of_track
        case 'line'
            plot((ilon-1)*DX,(ilat-1)*DY,'w--','linewidth',2.5);
        case 'scatter_altitude'
            scatter((ilon-1)*DX,(ilat-1)*DY,50,dat_flt19(inds,11),'filled');
    end
end

 iplot_return=1;
 if iplot_return==1 & just_labels==0
     ntot2 = length(time_flt19);
     dn=(ntot2-ntot)/np;
     inds = round([ntot:dn:ntot2]);
     
%     [a b]= findheight(time_flt19,20.4,20.75);
%     inds=round(a:dn:b); 

     LAT_plot = dat_flt19(inds,2);
     LON_plot = dat_flt19(inds,3);

     [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_plot,LON_plot,0.1);

     hold on
%      x_vals_flight = [x_vals_flight (ilon-1)*DX];
%      y_vals_flight = [y_vals_flight (ilat-1)*DY];
%      inds_flight = [inds_flight inds];

%switch type_of_track
%    case 'line'
         plot((ilon-1)*DX,(ilat-1)*DY,'w--','linewidth',2.5);
%    case 'scatter_altitude'
%         scatter((ilon-1)*DX,(ilat-1)*DY,50,time_flt19(inds,1));
%end
              
% height_col = dat_flt19(inds,11)/maxALL(dat_flt19(inds,11)); %scale to between 0 and 1
%      for iplot=1:length(ilon)
%          plot((ilon(iplot)-1)*DX,(ilat(iplot)-1)*DY,'o','color',[height_col(iplot) height_col(iplot) height_col(iplot)]);
%      end

 end
 
% dist_flight = cumsum(sqrt(diff(x_vals_flight).^2 + diff(y_vals_flight).^2)); %find distances between each point and the last and sum cumulatively
 %to get the distance from the start
% time_dist_flight = time_flt19(inds_flight(2:end)); %ignore the first index as did differences
 %now have distance vs. time for equally spaced points along the flight path for using to plot e.g. distance vs. a variable to check variabilty
 
 if just_track==1
     return
 end
 
 
 if i_plot_location_labels==1

     %LAT_plot=[-67.2 -67.2 -67.2 -68.0168 lat2d(1).var(140,240) lat2d(1).var(174,290)];  %places during descent plus other locations
     %LON_plot=[-62 -63 -64 -62.4159 lon2d(1).var(140,240) lon2d(1).var(174,290)];
     LAT_plot = LAT;
     LON_plot = LON; %ones left over from watervap...m
     [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_plot,LON_plot,0.1);

     abc='ABCDEFGHIJKLMNOP';

     clear labels
     for i=1:length(ilat)
         if i_use_set_labels==1
             labels(i,:) = location_lab(i).l;   %EAST side
         else
             labels(i,:) = ['' abc(i)];   %EAST side
         end
     end

     %----------------   plotting EAST SIDE POINTS   ---------------------------------------------%
     %plot((ilon-1)*DX,(ilat-1)*DY,'ws','markersize',20);
     %text((ilon-6)*DX,(ilat-1)*DY,labels,'color','k','fontweight','bold');

     switch marker_type
         case 'filled square'
             text((ilon-1)*DX,(ilat-1)*DY,labels,'color','k','fontweight','bold','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[.7 .9 .7]);
         case 'plus sign'
             for i=1:length(ilat)

                 %            labels(i,:) = ['o'];

             end
             text((ilon-5)*DX,(ilat-5)*DY,labels,'color',[0.3 0.3 0.3],'fontsize',13,'fontweight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
             plot((ilon-1)*DX,(ilat-1)*DY,'o','MarkerFaceColor','w','MarkerEdgeColor','k');
             plot((ilon-1)*DX,(ilat-1)*DY,'b+');
     end
     %text((ilon-1)*DX,(ilat-1)*DY,labels,'color','r','fontweight','bold','HorizontalAlignment','center','BackgroundColor',[.7 .9 .7]); %verical middle is the default
     %--------------------------------------------------------------------------------------------%

     LAT_plot=[-67.55 -67.62 -67.55]; %start, middle and end of ascent
     LON_plot=[-68.1 -67.8 -67.5];
     LAT_plot=[-67.55 -68.4]; %start of ascent and point before peninsula where flows come over
     LON_plot=[-68.1 -68.4];
     [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_plot,LON_plot,0.1);



     clear labels
     for i=1:length(ilat)
         labels(i,:) = ['W' abc(i)];
     end

     %----------------   plotting WEST SIDE POINTS   ---------------------------------------------%
     %plot((ilon-1)*DX,(ilat-1)*DY,'ws','markersize',20); %WEST side
     %text((ilon-6)*DX,(ilat-1)*DY,labels,'color','k','fontweight','bold');
     %--------------------------------------------------------------------------------------------%


 end

idraw_AWS=0;
if idraw_AWS==1
    LAT_AWS(1).dat = -68.34;  LON_AWS(1).dat =  -69.01;  short_name(1).name = 'Kirk';
    LAT_AWS(2).dat = -68.09;  LON_AWS(2).dat =  -68.82;  short_name(2).name = 'Dis';
    LAT_AWS(3).dat = -64.78;  LON_AWS(3).dat =  -64.07;  short_name(3).name = 'Bon';
    LAT_AWS(4).dat = -74.79;  LON_AWS(4).dat =  -71.49;  short_name(4).name = 'Sky';
    LAT_AWS(5).dat = -75.91;  LON_AWS(5).dat =  -59.26;  short_name(5).name = 'Lim';
    LAT_AWS(6).dat = -72.21;  LON_AWS(6).dat =  -60.17;  short_name(6).name = 'But';
    LAT_AWS(7).dat = -67.01;  LON_AWS(7).dat =  -61.55;  short_name(7).name = 'Lar';

    % (1) Kirkwood Island - Lat : 68.34S  Long :  69.01W    Elev :   30 M  %wind looks ok
    % (2) Dismal Island - Lat : 68.09S  Long :  68.82W      Elev :   10 M  %bad wind data
    % (3) Bonaparte Point -  Lat : 64.78S  Long :  64.07W   Elev :    8 M  %bad wind data
    % (4) Sky Blu   -   Lat : 74.79S  Long :  71.49W        Elev : 1510 M  %wind looks ok
    % (5) Limbert   -   Lat : 75.91S  Long :  59.26W        Elev :   40 M  %wind looks ok
    % (6) Butler Island  -  Lat : 72.21S  Long :  60.17W    Elev :   91 M  %not sure about wind data - think is ok


    %    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,-67.01,-61.55,0.1);
    for iaws=1:length(short_name)
        [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT_AWS(iaws).dat,LON_AWS(iaws).dat,0.1);

        %----------------   plotting AWS POINTS   ---------------------------------------------%
        %    plot((ilon-1)*DX,(ilat-1)*DY,'ws','markersize',20);
        %    text((ilon-6)*DX,(ilat-1)*DY,'AWS','color','k','fontweight','bold');

        %        text((ilon-6)*DX,(ilat-1)*DY,short_name(iaws).name,'color','k','fontweight','bold');
    end
end

% hline=line([timesTH(1).t(1);timesTH(1).t(pend)],[380 380]);
% set(hline,'color','k');
% 
% hline=line([timesTH(1).t(1);timesTH(1).t(pend)],[280 280]);
% set(hline,'color','k');
if exist('LAT') & exist('LON')
    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1); %restore the original ilat ilon values
end

if i_add_cross_section_line==1
    line(x_line,y_line,'color','k','linewidth',2);
end


 disp('Done plot aircraft path');