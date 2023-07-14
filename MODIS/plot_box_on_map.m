try
    
    if ~exist('ioverride_box_colour')==1 | ioverride_box_colour==0
        col_str='c-';  %doesn't like white for some reason, so use cyan
%col_str='k-';

        %lat_box = [72 75];
        %lon_box = [0 50.0];

        %lat_box = [72 75];
        %lon_box = [-3 48];

        %LAT_val = [-24.5 -15.44]; LON_val = [-86.93 -77.08];
        %LAT_val = [-20.25 -20]; LON_val = [-81.25 -81];


        %LON_val = [75:1:88]; LAT_val=[16:1:28];

        LAT_val = [-22.70 -17.28]; LON_val =[-78.93 -73.08]; %12th Nov UM FULL domain
        %LAT_val = [-20.5 -17.5]; LON_val = [-78.75 -73.25]; %12th Nov UM top half of domain

        LAT_val = [32.5 50]; LON_val =[-115 -80]; %12th Nov UM FULL domain
        
        itext_in_box=0;
        box_lwidth=2;
        
        irotated_pole_box=0;
        
        imap=1;

    end
    
if irotated_pole_box==1
    
   %Use the unrotated regular grid to define the box edges and then rotate
   %them to get lat,lon coords
   lat_bot = nest.lat(1)*ones([length(nest.lon) 1]);
   lat_top = nest.lat(end)*ones([length(nest.lon) 1]);
   lon_left = nest.lon(1)*ones([length(nest.lat) 1]);
   lon_right = nest.lon(end)*ones([length(nest.lat) 1]);
      
   [lat_bot2,lon_bot2]=em2gm(lat_bot,nest.lon,pole_lat,pole_lon);
   [lat_top2,lon_top2]=em2gm(lat_top,nest.lon,pole_lat,pole_lon);   
   [lat_left2,lon_left2]=em2gm(nest.lat,lon_left,pole_lat,pole_lon);
   [lat_right2,lon_right2]=em2gm(nest.lat,lon_right,pole_lat,pole_lon);   
   
   %Need to flip the top and left to keep the direction anti-clockwise
   xvec = [lon_bot2; lon_right2; flipud(lon_top2); flipud(lon_left2)];
   yvec = [lat_bot2; lat_right2; flipud(lat_top2); flipud(lat_left2)];   
   
   
else
    if iscell(LON_val)
        nbox=length(LON_val);
    else
        nbox=1;
    end
    
    for ibox=1:nbox        
        if iscell(LON_val)
            lon_box = LON_val{ibox};
        else
            lon_box = LON_val([1 end]);
        end
        
        lat_box = LAT_val([1 end]);
        
        xvec = [lon_box(1) lon_box(2) lon_box(2) lon_box(1) lon_box(1)];
        yvec = [lat_box(1) lat_box(1) lat_box(2) lat_box(2) lat_box(1)];
        
        if imap==1
            m_plot(xvec,yvec,col_str,'linewidth',box_lwidth);
        else
            plot(xvec,yvec,col_str,'linewidth',box_lwidth);
        end
        
        

                
        
    end
    

    
end



if itext_in_box==1
   
    box_centre_lon = 0.5*(lon_box(1) + lon_box(2) ) + offset_lon;
    box_centre_lat = 0.5*(lat_box(1) + lat_box(2) ) + offset_lat; 
    
    h_text = m_text(box_centre_lon,box_centre_lat,box_text,'HorizontalAlignment','center','VerticalAlignment','middle');
    
end

clear ioverride_box_colour
catch error_box
    clear ioverride_box_colour
    rethrow(error_box);
end

