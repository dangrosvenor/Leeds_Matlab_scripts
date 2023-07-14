function dat2 = reduce_high_lats(dat,LAT,operation)


%fractions that want to use (using fractions that give whole number
%divisions into 360 degrees)
fs=[1 1/2 1/3 1/4 1/5 1/6 1/8 1/9];

%mid-points - using these as the boundaries between the different averaging
mfs=0.5*(fs(2:end)+fs(1:end-1));


% lat_N(1)=floor(acos(3/4)*180/pi);
% lat_N(2)=floor(acos(1/2)*180/pi);
% lat_N(3)=floor(acos(1/3)*180/pi);
% lat_N(4)=floor(acos(1/4)*180/pi);
% lat_N(5)=floor(acos(1/5)*180/pi);
% lat_N(6)=floor(acos(1/6)*180/pi);
% lat_N(7)=floor(acos(1/9)*180/pi);

clear lat_N
for i=1:length(mfs)
    lat_N(i)=floor(acos(mfs(i))*180/pi);
end

switch operation
    case 'reduce'

        dat2 = NaN*ones(size(dat));
        

        for ilat=1:180
            i=find(lat_N<abs(LAT(ilat)));
            if length(i)>0
                fsi = fs(i(end)+1);
            else
                fsi=1;
            end
            i0=1;
            for ilon2=1:360*fsi
                        for it=1:size(dat,3)
                            dat2(ilat,i0:i0 + 1/fsi - 1,it) = mean(dat( ilat,i0:i0 + 1/fsi - 1, it ));
                        end
                        i0=i0+1/fsi;
            end
        end


    case 'expand'

        dat2 = NaN*ones([180 360]);

        for ilat=1:180
            i=find(lat_N<abs(LAT(ilat)));
            if length(i)>0
                fsi = fs(i(end)+1);
            else
                fsi=1;
            end
            i0=1;
            for ilon2=1:360*fsi
                for it=1:size(dat,3)
                    dat2(ilat,i0:i0 + 1/fsi - 1 ,it) = dat(ilat,ilon2,it);
                end
                i0=i0+1/fsi;
            end
        end

end


        
        
        
end




