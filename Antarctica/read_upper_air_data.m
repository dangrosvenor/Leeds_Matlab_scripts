%filename='/home/mbexddg5/work/Ant_Jan06/roth';
filename='c:/documents and settings/dan/my documents/Antarctica_general/roth';
fid=fopen(filename,'rt');
upper_text=fgetl(fid);
upper_text2=fgetl(fid);

dat_roth=fscanf(fid,'%f'); %make sure that when requesting the data from the BAS website that
                      %select 'type your own' date format and add a space before the year
                      %otherwise it runs the last data point for each pressure into the year

dat_roth=reshape(dat_roth,[10 length(dat_roth)/10]);

dat_roth(10,:)=dat_roth(10,:)*0.5144444;  %THINK the wind speed is in knots - says so for the manned station data but
 %not specifically for the upper air data : http://www.antarctica.ac.uk/met/metlog/terms_explained.html


%fields are
% 1) YEAR 
% 2) MONTH 
% 3) DAY 
% 4) HOUR 
% 5) PRESSURE 
% 6) HEIGHT 
% 7) TEMPERATURE 
% 8) DEWPOINT 
% 9) WIND_DIRECTION 
%10) WIND_SPEED


'done read of upper air data'                      
                      
