function [iz_1,iz_2]=findheight_nearest(z,z1,z2,bounds)
%function [iz1,iz2]=findheight_nearest(z,z1,z2) - z=array z1 and z2 are
%lower and upper bounds. The bounds variable is optional and is set to
%'bounded' if we want the returned indices to bound the [z1 z2] range rather
%than be the nearest index to z1 and z2. Default is to return the nearest.
%Can run without z2 set (just search for the location of one value, z1) -
%for that can either run with just two input arguments, or if want to run
%as 'bounded' then set z2=NaN

if nargin<4
    bounds = '';
end



for i=1:length(z1)
    if nargin~=2 & ~isnan(z2(1)) %if have more than value defining the range
        if z2(i)<z1(i)
            fprintf(1,'***z1 and z2 should be in ascending order***\n');
            return
        end

        if z2(i)<z(1)  %in this case we know that both z1 and z2 are lower than z(1) (since z2>z1)
            iz_1(i)=0;  %set both to zero to denote both below z range
            iz_2(i)=0;
            continue
        end

    end
 
    if z1(i)>z(end) %both larger than range of z
        iz_1(i)=1e20; 
        iz_2(i)=1e20;
        continue
    end



    if z1(i)<z(1)
        switch bounds
            case 'bounded'
                iz_1(i)=0; %setting to zero allows user to know it was out of bounds at lower end
            otherwise
                iz_1(i) = 1;
        end
    elseif z1(i)>z(end)
        switch bounds
            case 'bounded'
                iz_1(i)=1e20; %setting to this allows user to know it was out of bounds at upper end
            otherwise
                iz_1(i) = length(z);
        end


    else %otherwise find the location within

        %[temp iz_1(i)]=1; %in this case z1 is outside of the z array - i1=1 is the nearest value

        iz1=find(z<=z1(i)); %find all indices for z that are lower or equal to z1
        [temp iz]=max(z(iz1)); %find the index for the max value, which will be the one we want, or possibly next to it
        %if we want the nearest value
        iz_1(i)=iz1(iz);

        switch bounds
            case 'bounded'
            otherwise %if want nearest value
                if iz_1(i)~=length(z) %if we are at the end of the array then this will be our closest value
                    %otherwise find which index is closest
                    if abs(z(iz_1(i))-z1(i)) > abs(z(iz_1(i)+1)-z1(i))
                        iz_1(i)=iz_1(i)+1;
                    end
                end
        end


    end



    if nargin~=2 & ~isnan(z2(1)) %if have more than value defining the range
        if z2(i)<z(1)
            switch bounds
                case 'bounded'
                    iz_2(i)=0; %setting to zero allows user to know it was out of bounds at lower end
                otherwise
                    iz_2(i) = 1;
            end
        elseif z2(i)>z(end)
            switch bounds
                case 'bounded'
                    iz_2(i)=1e20; %setting to this allows user to know it was out of bounds at upper end
                otherwise
                    iz_2(i) = length(z);
            end
            
        else


            iz2=find(z>=z2(i));
            [temp iz]=min(z(iz2));
            iz_2(i)=iz2(iz);


            switch bounds
                case 'bounded'
                otherwise %if want nearest value (otherwise we already have our bounded one)
                    if iz_2(i)~=length(z) %if we are at the end of the array then this will be our closest value
                        %otherwise find which index is closest
                        if abs(z(iz_2(i))-z2(i)) > abs(z(iz_2(i)-1)-z2(i))
                            iz_2(i)=iz_2(i)-1;
                        end
                    end
            end




        end


    end

end %for loop
