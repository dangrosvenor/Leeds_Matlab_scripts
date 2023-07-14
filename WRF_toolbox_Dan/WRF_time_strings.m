function [year,month,day,hour,min,sec,month_text]=WRF_time_strings(tstr)

%tstr=Times(time,:); 

year=tstr(1:4);
month=tstr(6:7);
day=tstr(9:10);
hour=tstr(12:13);
min=tstr(15:16);
sec=tstr(18:19);


month_strings = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
month_text = month_strings{str2num(month)};

