function [date,date_num] = date_from_day_of_year_func(day,year)
%function [date_str,date_num ] = date_from_day_of_year_func(day,year)
% 
%Also see  day = day_of_year_from_date_func(date_str)


for i=1:length(day)
    if isnan(year(i))==1 | isnan(day(i))==1
        date_num(i)=NaN;
        date{i}='NaN';
        continue
    end
date_num(i) = datenum(['01-Jan-' num2str(year(i))]) + day(i) - 1 ;
date{i} = datestr( date_num(i) );
end

if length(day)==1
    date = date{1};
end