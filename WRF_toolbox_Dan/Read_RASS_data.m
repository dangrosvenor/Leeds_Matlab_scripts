filename='C:\Documents and Settings\dan\My Documents\WRF\Yanbu_upper_air_RASS_data\20090313.dat';

fid=fopen(filename,'rt');

clear hour


for time=1:24
    time
    
    line_text = fgetl(fid);
    line_text = fgetl(fid);
    
    time_stamp = fscanf(fid,'%f',[1 10]);
    hour(time) = time_stamp(5);
    
    for i=1:10
        line_text = fgetl(fid);
    end
    
    data=fscanf(fid,'%f',[13 30]);
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
