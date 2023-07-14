function sza=sun_pos(time,lat,lon,alt)
%function sza=sun_pos(time,lat,lon)
%***time is in UTC***
%wrapper for the sun_position script (due to its awkward input method)
%This allows it be run using numerical input
%time is in Matlab format time. time can also be in '01-May-2001 12:00'
%format. Or it can be a vector of Matlab date numbers.

if isstr(time)
    time = datenum(time);
end



if length(lat)>2 | ( length(lat)==2 & min(size(lat))==1 )
   siz = size(lat);
   lat=lat(:);
   lon=lon(:);
   time=time(:);
   sza=ones*NaN(size(lat)); %pre allocate to speed up and save memory
   if nargin<4
       alt=zeros([1 length(lat)]);
   else
       alt=alt(:);
   end
else
     if nargin<4
       alt=zeros([1 length(lat)]);
     end
end



if length(lat)==1
    if ~isstr(time)
        nt=length(time);
    else
        nt=1;
    end
    
    lat = repmat(lat,[1 nt]);
    lon = repmat(lon,[1 nt]);  
    alt = repmat(alt,[1 nt]);   
    siz = size(lat);
else
    nt=length(lat);
end

for i=1:nt
    if ~isstr(time)
        if isnan(time(i))
            sza(i)=NaN;
            continue
        else
            t=datestr(time(i),31);
        end
    else
        t=time';
    end


    location.latitude=lat(i);
    location.longitude=lon(i);
    location.altitude=alt(i);

    output=sun_position(t,location);

    sza(i) = output.zenith;

end

sza = reshape(sza,siz);


