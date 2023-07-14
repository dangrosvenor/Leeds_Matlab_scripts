%This program reads the RASS/SODAR Data


filename='C:\Documents and Settings\dan\My Documents\MATLAB\Jaafars scripts\20100101.dat'; %The path of RASS Data

fid=fopen(filename,'rt');

clear hour


for time=1:24   %This do loop for 24 hours (old RASS data it works upto 24 but new format upto 21)
    time
 
    line_text = fgetl(fid);
    line_text = fgetl(fid);
 
    time_stamp = fscanf(fid,'%f',[1 10]);        % means one column and ten raws which is the header of the file
    hour(time) = time_stamp(5);                 % the location of the hour
 
    for i=1:10                                  %Do loop for the heading of each hour
        line_text = fgetl(fid);                 %means first ten raws in this do loop is text
    end
 
    data=fscanf(fid,'%f',[14 30]); %Do loop for the data (old RASS data 13 column and 30 raws but new format 14 column)
    ibad=find(data==-9999);
    data(ibad)=NaN;
 
    data_structure(time).data=data;
 
    for i=1:2
        line_text = fgetl(fid);
    end
 
end


% 1 ALT
% 2 CT
% 3 SPEED
% 4 DIR
% 5 S DIR
% 6 W
% 7 SW
% 8 INVMI
% 9 STAB
% 10 DTDZ
% 11 ECH T
% 12 T
% 13 ST

fprintf(1,'\nFinished reading RASS data\n');

