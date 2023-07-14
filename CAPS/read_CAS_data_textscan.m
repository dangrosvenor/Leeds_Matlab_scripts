%much faster version (compared to repeated use of textread as in read_CAS_data.m)
%stores the data in cell arrays with different names
%datevals_CAS AD_vals_CAS flag_vals_CAS Counts_CAS Bchannel_CAS Psep_CAS datevals_CIP AD_vals_CIP flag_vals_CIP Counts_CIP
%each one has N cells for the n different times read in
%e.g. Counts_CAS{100}{1} will be all the 40 counts for time/particle? 100



clear datevals_CAS AD_vals_CAS flag_vals_CAS Counts_CAS ...
    Bchannel_CAS Psep_CAS datevals_CIP AD_vals_CIP flag_vals_CIP Counts_CIP ...
    

fileName='Y:\BAS_flights\flight99_6thFeb2010\flight99_020610_150001.txt';
%flight99_020610_140001.txt

fid=fopen(fileName,'rt');

%skip 12 lines to get to the beam diameter info
for i=1:12
    textline=fgetl(fid);
end

beam_line = textscan(fid,'%f %f %f',1,'delimiter',',');
beam_diameter = beam_line{1}; %quoted as 0.24 in the datafiles

%Jonny's probe has A=0.18 mm^2. This corresponds to a diameter of 0.48 mm - so twice that of this older probe.
%BAS beam diameter is 0.24 (think is mm). So area = pi*0.24^2/4 = 1.3201e-4 mm^2
%sample_vol = A*air_speed*time
%then divide the particle count in the time period by the sample_vol to get particle concentration
%will need to get airspeed from the other aircraft data?
%time looks to be 1 Hz data (data is counted up and reset every second)

go=1;
last_CAS_missing=0;
last_CIP_missing=0;

%find out where the size range info starts (is two lines below "Laser Power: Low" string)
while go==1
    textline=fgetl(fid);
    if strcmp(textline,'"Laser Power: Low"')==1
        go=0;
    end
end

[range_info_CAS]=textscan(fid,...
        '%d %d %d',1,'delimiter',',');   %data is range, channel count, user table
    
channel_count = range_info_CAS{2};
[bins_CAS]=textscan(fid,...
        '%d %f',channel_count,'delimiter',',');   %data columns are threshold and size 
                %threshold is the voltage threshold and so not really needed
                %size data will then be in bins_CAS{2} and is diameter in microns.
                %Note have 30 bins for forwards scatter (30 for back).
                
[range_info_CAS_back]=textscan(fid,...
        '%d %d %d',1,'delimiter',',');   %data is range, channel count, user table
    
channel_count = range_info_CAS_back{2};
[bins_CAS_back]=textscan(fid,...
        '%d %f',channel_count,'delimiter',',');   %data columns are threshold and size 
                %threshold is the voltage threshold and so not really needed
                %size data will then be in bins_CAS{2} and is diameter in microns.                



%find out where the data starts (is labelled "DATA")
go=1;
while go==1
    textline=fgetl(fid);
    if strcmp(textline,'"DATA"')==1
        go=0;
    end
end
    

idat=0; %the current data record (time index)

go=1;
while go==1 

    idat=idat+1;
    
%%%%%%% Read in data for probe_num=0 --- CAS instrument ---   %%%%%%

    
    [datevals_CAS{idat}]=textscan(fid,...
        '%d %s %d %s %d %d',1,'delimiter',',');   %see above for a list of what these are
    
    [AD_vals_CAS{idat}]=textscan(fid,...   %AD_values are things like instrument voltages, etc. - housekeeping, not needed
        '%f',31,'delimiter',','); 
  
    [flag_vals_CAS{idat}]=...   %flags -see above
        textscan(fid,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',','); 

    %particle counts in channels 0-29. Note have 30 forward and 30 backscatter channels
    [Counts_CAS{idat}]=textscan(fid,...
        '%f',30,'delimiter',','); 
    
    if size(Counts_CAS{idat}{1},1)==0 %when runs out of data to read it returns a zero length variable
        last_CAS_missing=1;
        go=0; %exit the loop
    end

    [Bchannel_CAS{idat}]=textscan(fid,...   %back scatter counts - probably not needed
        '%d',30,'delimiter',',');  

    [Psep_CAS{idat}]=textscan(fid,...       %Jonny thinks this is a histogram of interarrival times - useful to see if have shattering?
        '%f',64,'delimiter',','); 

    
%%%%%%  %now for probe_num=1 --- CIP instrument ---  %%%%

    [datevals_CIP{idat}]=textscan(fid,...
        '%d %s %d %s %d %d',1,'delimiter',','); 

    [AD_values_CIP{idat}]=textscan(fid,...
        '%f',31,'delimiter',','); 

    [flag_vals_CIP{idat}]=...
        textscan(fid,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',','); 

    %particle counts in channels 0-61. Note - this is for the CIP
    [Counts_CIP{idat}]=textscan(fid,...
        '%f',62,'delimiter',','); 
    
    if size(Counts_CIP{idat}{1},1)==0 %when runs out of data to read it returns a zero length variable
        last_CIP_missing=1;
        go=0; %exit the loop
    end
                       
end

fclose(fid);

CAS_bins=bins_CAS{2};

clear CAS_counts CAS_time CAS_psep
L_CAS=length(Counts_CAS);
if last_CAS_missing==1
    L_CAS=L_CAS-1;
end

i=1;
date_str=datevals_CAS{1}{2}{1}(2:11);
CAS_counts(i,:)=Counts_CAS{i}{1};
CAS_psep(i,:)=Psep_CAS{i}{1};
time_str=datevals_CAS{i}{2}{1}(13:20);
hours=str2num(time_str(1:2));
mins=str2num(time_str(4:5));
secs=str2num(time_str(7:8));
CAS_time(i)=hours*3600+mins*60+secs; %calculate the time in seconds
CAS_start_time = CAS_time(i);
CAS_start_time_str = time_str;

for i=2:L_CAS
    CAS_counts(i,:)=Counts_CAS{i}{1};
    CAS_psep(i,:)=Psep_CAS{i}{1};
    time_str=datevals_CAS{i}{2}{1}(13:20);
    hours=str2num(time_str(1:2));
    mins=str2num(time_str(4:5));
    secs=str2num(time_str(7:8));
    CAS_time(i)=hours*3600+mins*60+secs; %calculate the time in seconds
end

%[min,sec,hPa,gpm,T,RH,Dewp]=textscan(fid,...
%    '%d %s %d %s %d %d %f %*[^\n]','headerlines',42);


