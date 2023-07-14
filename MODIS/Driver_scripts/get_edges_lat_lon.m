function [lat_out,lon_out] = get_edges_lat_lon(lat,lon)
%[lat_out,lon_out] = get_edges_lat_lon(lat,lon)
%Returns the edges for a 2D lat lon grid - assumes same spacing at the
%edges. lat and lon are 2D grids so that lat(i,j) and lon(i,j) are the lat
%and lon of location (i,j)


%Check for monotonicity by sorting and then checking that nothing has
%changed
[B,I] = sort(lon,2);
% for i=1:size(lon,1)
%     inds = I(i,:);
%     C(i,:) = lon(i,inds);
% end
m=maxALL(abs(lon-B));

%If not monotonic try to phase-shift to make it so - could be that the lons
%start at 0, but then run to 180 and then -180 (as in the UM). In which case converting to 0-360 type would make it montonic.
%Or that they start at 180, then run to close to 360 and then drop to 0 (i.e., Greenwich in the middle with 0 to
%360 notation) - converting to -180 to 180 will sort this type out.
%Wnat to make sure we convert back at the end.
i360=0; %Flags to say what type of lon it was:-0-360 or -180 to 180
i180=0;
if m>0    
    if minALL(lon)>0
        i = find(lon>180);
        lon(i) = lon(i) - 360;
        i360=1;
    else        
        i = find(lon<0);
        lon(i)=lon(i) + 360;
        i180=1;
    end
    
end

%Check again for monotonicity
[B,I] = sort(lon,2);
m=maxALL(abs(lon-B));
if m>0
    error('Cannot make monotonic');
end


lat_out = get_edges_fun(lat);
lon_out = get_edges_fun(lon);

if i360==1
   i = find(lon_out<0);
   lon_out(i) = lon_out(i) + 360;
end
if i180==1
   i0 = find(lon_out>180);
   lon_out(i0) = lon_out(i0) - 360;
end
       
%lat_out=NaN*ones([size(lat,1)+1 size(lat,2)+1]);
%lon_out=NaN*ones([size(lat,1)+1 size(lat,2)+1]);

% %Do latitude first
%         %Find the mid-points of the cell centres - these will be the edges
%         %for the central part of the grid
%         lat2 = 0.5*( lat(1:end-1,:)+lat(2:end,:) );
%         lat2 = 0.5*( lat2(:,1:end-1)+ lat2(:,2:end) );
%         %Put in the central part
%         lat_out(2:end-1,2:end-1)=lat2;
%         %Add the top edge
%         dlat = diff(lat(1:2,1:end-1));
%         lat_out(1,2:end-1) = lat_out(2,2:end-1) - dlat;        
%         %Add the bottom edge
%         dlat = diff(lat(1:2,1:end-1));
%         lat_out(end,2:end-1) = lat_out(end-1,2:end-1) + dlat;  
%         
%         %Now we have defined the full size in the row direction, but not in
%         %the column direction
%         
%         %Add the left edge
%         dlat = diff(lat_out(:,2:3),[],2);
%         lat_out(:,1) = lat_out(:,2) - dlat;
%         
%         %Add the right edge
%         dlat = diff(lat_out(:,end-2:end-1),[],2);
%         lat_out(:,end) = lat_out(:,end-1) + dlat;  
%         
%         
% %Do longitude now
%         %Find the mid-points of the cell centres - these will be the edges
%         %for the central part of the grid
%         lon2 = 0.5*( lon(1:end-1,:)+lon(2:end,:) );
%         lon2 = 0.5*( lon2(:,1:end-1)+ lon2(:,2:end) );
%         %Put in the central part
%         lon_out(2:end-1,2:end-1)=lon2;
%         %Add the top edge
%         dlon = diff(lat(1:2,1:end-1));
%         lat_out(1,2:end-1) = lat_out(2,2:end-1) - dlat;        
%         %Add the bottom edge
%         dlat = diff(lat(1:2,1:end-1));
%         lat_out(end,2:end-1) = lat_out(end-1,2:end-1) + dlat;  
%         
%         %Now we have defined the full size in the row direction, but not in
%         %the column direction
%         
%         %Add the left edge
%         dlat = diff(lat_out(:,2:3),[],2);
%         lat_out(:,1) = lat_out(:,2) - dlat;
%         
%         %Add the right edge
%         dlat = diff(lat_out(:,end-2:end-1),[],2);
%         lat_out(:,end) = lat_out(:,end-1) + dlat; 
        
      
    function lat_out = get_edges_fun(lat)
        lat_out=NaN*ones([size(lat,1)+1 size(lat,2)+1]);
        
        %Do latitude first
        %Find the mid-points of the cell centres - these will be the edges
        %for the central part of the grid
        lat2 = 0.5*( lat(1:end-1,:)+lat(2:end,:) );
        lat2 = 0.5*( lat2(:,1:end-1)+ lat2(:,2:end) );
        %Put in the central part
        lat_out(2:end-1,2:end-1)=lat2;
        %Add the top edge
        dlat = diff(lat(1:2,1:end-1));
        lat_out(1,2:end-1) = lat_out(2,2:end-1) - dlat;        
        %Add the bottom edge
        dlat = diff(lat(1:2,1:end-1));
        lat_out(end,2:end-1) = lat_out(end-1,2:end-1) + dlat;  
        
        %Now we have defined the full size in the row direction, but not in
        %the column direction
        
        %Add the left edge
        dlat = diff(lat_out(:,2:3),[],2);
        lat_out(:,1) = lat_out(:,2) - dlat;
        
        %Add the right edge
        dlat = diff(lat_out(:,end-2:end-1),[],2);
        lat_out(:,end) = lat_out(:,end-1) + dlat; 