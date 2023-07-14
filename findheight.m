function [iz_1,iz_2]=findheight(z,z1,z2)
%function [iz1,iz2]=findheight(z,z1,z2) - z=array z1 and z2 are lower and upper bounds
%returns NaNs if z1 to z2 is not in the array z
%returns iz_1=0 if z(1) is between z1 and z2, but z(1)>z1
%returns iz_2=1e20 if z(end) is between z1 and z2, but z(end)<z2

error('Need to fix findheight - use findheight_nearest');


for i=1:length(z1)
    
    if nargin==3

        if z2(i)<z1(i)
            fprintf(1,'***z1 and z2 should be in ascending order***\n');
            return
        end

        if z2(i)<z(1) %in this case we know that both z1 and z2 are lower than z(1) (since z2>z1)
            iz_1(i)=0; %set both to zero to denote both below z range
            iz_2(i)=0;
            continue
        end
        
    end

    if z1(i)>z(end) %both larger than range of z
        iz_1(i)=1e20;
        iz_2(i)=1e20;
        continue
    end

    iz1=find(z>=z1(i));
%    [temp iz]=min(z(iz1));
    [iz]=min(iz1);    
    if iz==1 & z(iz1(iz))>z1
        iz_1(i)=0; %here we know that part of the requested range is in bounds - so user can now make a choice
        %of whether they want to use the indices even though some of them
        %are not in bounds
    else
        iz_1(i)=iz1(iz);
    end
    

    if nargin==3
        iz2=find(z<=z2(i));
        [temp iz]=max(z(iz2));
        if iz==length(z) & z(iz2(iz))<z2
            iz_2(i)=1e20;
        else
            iz_2(i)=iz2(iz);
            
        end

            
        
        
    end

end
