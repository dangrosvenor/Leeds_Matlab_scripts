function [day,day_str] = day_of_year_from_date_func(date_str,date_format)
% [day,day_str] = day_of_year_from_date_func(date_str,date_format)
%date_format is not required if using 'dd-mmm'yyyy' format
%Can input date_str as a number from datenum too
%also see date = date_from_day_of_year_func(day,year)


if nargin==1
    date_format = 'dd-mmm-yyyy';
end

if isstr(date_str)
    date_num = datenum(date_str,date_format);
else
    date_num = date_str;
end

for iday=1:length(date_num)

    [Y,M,D] = datevec(date_num(iday));

    date_str_beg = repmat(['01-Jan-'],[length(Y) 1]);

    date_str2 = cat(2,date_str_beg,num2str(Y));

    day(iday) = floor(date_num(iday) -  datenum(date_str2) + 1);

    if day(iday)<10
        day_str{iday} = ['00' num2str(day(iday))];
    elseif day<100
        day_str{iday} = ['0' num2str(day(iday))];
    else
        day_str{iday} = [num2str(day(iday))];
    end

end

if length(day_str)==1
    day_str = day_str{1};
end