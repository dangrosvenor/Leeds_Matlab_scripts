function julian = julian_calculation_Dan(year_vec,month_vec,day_vec,hour_vec,min_vec,sec_vec)
% This function compute the julian day and julian century from the local
% time and timezone information. Ephemeris are calculated with a delta_t=0
% seconds. 

utc_vec=0;

% % If time input is a Matlab time string, extract the information from
% % this string and create the structure as defined in the main header of
% % this script.
% if ~isstruct(t_input)
% %    tt = datevec(t_input);
%     utc_vec=0;
% %     year_vec = tt(1);
% %     month_vec = tt(2);
% %     day_vec = tt(3);
% %     hour_vec = tt(4);
% %     min_vec = tt(5);
% %     sec_vec = tt(6);
%  [year_vec,month_vec,day_vec,hour_vec,min_vec,sec_vec] = datevec(t_input);
%     
% else
%     time = t_input;
% end

if(month_vec == 1 | month_vec == 2)
    Y = year_vec - 1;
    M = month_vec + 12;
else
    Y = year_vec;
    M = month_vec; 
end
ut_time = ((hour_vec - utc_vec)/24) + (min_vec/(60*24)) + (sec_vec/(60*60*24)); % time of day in UT time. 
D = day_vec + ut_time; % Day of month in decimal time, ex. 2sd day of month at 12:30:30UT, D=2.521180556


% In 1582, the gregorian calendar was adopted
if(year_vec == 1582)
    if(month_vec == 10)
        if(day_vec <= 4) % The Julian calendar ended on October 4, 1582
            B = 0;    
        elseif(day_vec >= 15) % The Gregorian calendar started on October 15, 1582
            A = floor(Y/100);
            B = 2 - A + floor(A/4);    
        else
            disp('This date never existed!. Date automatically set to October 4, 1582');
            month_vec = 10;
            day_vec = 4; 
            B = 0;
        end
    elseif(month_vec<10) % Julian calendar 
        B = 0;
    else % Gregorian calendar
        A = floor(Y/100);
        B = 2 - A + floor(A/4);
    end
    
elseif(year_vec<1582) % Julian calendar
    B = 0;
else
    A = floor(Y/100); % Gregorian calendar
    B = 2 - A + floor(A/4);
end

julian.day = floor(365.25*(Y+4716)) + floor(30.6001*(M+1)) + D + B - 1524.5;

delta_t = 0; % 33.184;
julian.ephemeris_day = julian.day + (delta_t/86400);

julian.century = (julian.day - 2451545) / 36525; 

julian.ephemeris_century = (julian.ephemeris_day - 2451545) / 36525;

julian.ephemeris_millenium = julian.ephemeris_century / 10; 
