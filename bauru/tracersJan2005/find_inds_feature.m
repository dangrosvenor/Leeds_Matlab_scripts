function [yinds,xinds]=find_inds_feature(Grid,X,Y,D)


y=Grid.Y1(2:end-1)';
dy=diff(y(1:2));

ywrapped=[ y(1)-dy*length(y) : dy : y(1)-dy   y   y(end)+dy : dy : y(end)+dy*length(y) ];
y_inds_wrapped=[2:length(y)+1 2:length(y)+1 2:length(y)+1];

[iy,iy2]=findheight(ywrapped,Y-D/2,Y+D/2);
yinds=[y_inds_wrapped(iy:iy2)];


%x indices


y=Grid.X1(2:end-1)';
if length(y)>3
    dy=diff(y(1:2));
    
    ywrapped=[ y(1)-dy*length(y) : dy : y(1)-dy   y   y(end)+dy : dy : y(end)+dy*length(y) ];
    y_inds_wrapped=[2:length(y)+1 2:length(y)+1 2:length(y)+1];
    
    [iy,iy2]=findheight(ywrapped,X-D/2,X+D/2);
    xinds=[y_inds_wrapped(iy:iy2)];
    
end

'end';