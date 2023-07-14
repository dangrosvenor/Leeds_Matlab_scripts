function sza=sun_pos2(Y,M,day_of_month,H,MI,S,lat,lon)
%function sza=sun_pos2(Y,M,day_of_month,H,MI,S,lat,lon)
%***time is in UTC***
%wrapper for the sun_position script (due to its awkward input method)
%This only allows the above format - is quicker than the previous method
%that called datevec.
% 
% if ~isstr(time)
%     nt=length(time);
% else
%     nt=1;
% end
% 
 for i=1:length(Y)
%     if ~isstr(time)
%         t=datestr(time(i),31);
%     else
%         t=time;
%     end


    location.latitude=lat;
    location.longitude=lon;
    location.altitude=0;

    output=sun_position_Dan(Y,M,day_of_month,H,MI,S,location);

    sza(i) = output.zenith;

end

