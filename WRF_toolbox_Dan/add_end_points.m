function [x,y]=add_end_points(x,y,Lx,Ly)

%remove NaN points
        x(isnan(x))=[];
        y(isnan(y))=[];


%remove all points at the edges of the plot as these go wrong
        x(y==1)=[];
        y(y==1)=[];

        x(y==Ly)=[];
        y(y==Ly)=[];

        y(x==1)=[];
        x(x==1)=[];

        y(x==Lx)=[];
        x(x==Lx)=[];


if length(x)>1 %can only extrapolate if have at least two points
        % find the beginning and end location by linear extrapolation

%find points just off end and beginning of the line
                dx = mean(diff(x));
                x_start = x(1)-dx;
                x_end   = x(end)+dx;

                y_start = interp1(x,y,x_start,'linear','extrap');
                y_end = interp1(x,y,x_end,'linear','extrap');

% interpolate points beyond the screen at the screen edges
                if x_start<1
                        x_start = 1;
                        y_start = interp1(x,y,1,'linear','extrap');
                end

                if y_start>Ly
                        y_start = Ly;
                        x_start = interp1(y,x,Ly,'linear','extrap');
                end

                if x_end>Lx
                        x_end = Lx;
                        y_end = interp1(x,y,Lx,'linear','extrap');
                end

                if y_start<1
                        y_start = 1;
                        x_start = interp1(y,x,1,'linear','extrap');
                end



% are plotting along the lon axis so put in end points accordingly
                x = [x_start x x_end];
                y = [y_start y y_end];

end
