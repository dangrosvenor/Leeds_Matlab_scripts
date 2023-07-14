function [iz_1,iz_2]=findheight_nearest(z,z1,z2,bounds)
%function [iz1,iz2]=findheight_nearest(z,z1,z2) - z=array z1 and z2 are
%lower and upper bounds. The bounds variable is optional and is set to
%'bounds' if we want the returned indices to bound the [z1 z2] range rather
%than be the nearest index to z1 and z2. Default is to return the nearest

if nargin<4
    bounds = '';
end

for i=1:length(z1)
            if z1<z(1)
                iz_1(i) = 1;
            elseif z1>z(end)
                iz_1(i) = length(z);
                continue
            end
            
        %[temp iz_1(i)]=1; %in this case z1 is outside of the z array - i1=1 is the nearest value
        
    iz1=find(z<=z1(i)); %find all indices for z that are lower or equal to z1
    [temp iz]=max(z(iz1)); %find the index for the max value, which will be the one we want, or next to it
    if length(iz)==0

    else
        iz_1(i)=iz1(iz);

        if iz_1(i)~=length(z)
            if abs(z(iz_1(i))-z1(i)) > abs(z(iz_1(i)+1)-z1(i))
                iz_1(i)=iz_1(i)+1;
            end
        end

    end
           
%    iz_1(i)=iz1(iz);
   
    if nargin>=3
        iz2=find(z>=z2(i));
        [temp iz]=min(z(iz2));
        if length(iz2)==0
            [temp iz_2(i)]=length(z);
        else
            iz_2(i)=iz2(iz);
        end
    end
    

end
