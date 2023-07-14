%much faster version (compared to repeated use of textread as in read_CAS_data.m)
%stores the data in cell arrays with different names
%datevals_CAS AD_vals_CAS flag_vals_CAS Counts_CAS Bchannel_CAS Psep_CAS datevals_CIP AD_vals_CIP flag_vals_CIP Counts_CIP
%each one has N cells for the n different times read in
%e.g. Counts_CAS{100}{1} will be all the 40 counts for time 100

clear datevals datevals_CAS AD_vals_CAS flag_vals_CAS Counts_CAS ...
    Bchannel_CAS Psep_CAS datevals_CIP AD_vals_CIP flag_vals_CIP Counts_CIP ...
    press_CAS dynamic_CAS lwc_CAS rh_CAS temp_CAS press_CIP dynamic_CIP lwc_CIP rh_CIP temp_CIP
    
%fileName='Y:\BAS_flights\flight99_6thFeb2010\flight99_020610_150001.txt';
%fileName='Y:\BAS_flights\flight99_6thFeb2010\flight99b_020610_160001.txt';
%fileName='Y:\BAS_flights\flight102\flight102_021110_200000.txt';
%fileName='Y:\BAS_flights\CAPS_cals_02-02-10\2-02-10_CAPS_test_15micron.txt';
%fileName='Y:\BAS_flights\CAPS_cals_02-02-10\2-02-10_CAPS_test_30 micron.txt';
%fileName='Y:\BAS_flights\CAPS_cals_02-02-10\2-02-10_CAPS_test_Lycopodium.txt';
%fileName='Y:\BAS_flights\CAPS_cals_02-02-10\2-02-10_CAPS_test_nebulizer.txt';
%fileName='Y:\BAS_flights\CAPS_cals_02-02-10\2-02-10_CAPS_test_spind.txt'; %blank distribution
%fileName='Y:\BAS_flights\CAPS_cals_02-02-10\2-02-10_CAPS_test_no_beads_2.txt'; %blank distribution
%fileName='Y:\BAS_flights\flight104\flight104_021210_200000.txt';
%fileName='Y:\BAS_flights\flight104\flight104_021210_190001.txt';
%fileName='Y:\BAS_flights\flight104\flight104_021210_210000.txt';
%fileName='Y:\BAS_flights\flight102\flight102_021110_200000.txt';


istr=findstr(fileName,'\flight'); %finds the location of the string '\flight'
if length(istr)==0
    istr=findstr(fileName,'/flight'); %finds the location of the string '\flight'
end
no_flight_data=0;
if length(istr)>0
    flight_no=fileName(istr(2)+7:istr(2)+9); %this should be followed by the flight number
    if strcmp(flight_no(end),'_')==1
        flight_no=flight_no(1:end-1);
    end
else
    no_flight_data=1;
    flight_no='Calibration';
    dat_flt=[];
end

%flight99_020610_140001.txt

fid=fopen(fileName,'rt');

%skip 12 lines to the beam area info
for i=1:12
    textline=fgetl(fid);
end

beam_line = textscan(fid,'%f %f %f',1,'delimiter',',');
beam_area = beam_line{1};


%sample_vol = A*air_speed*time
%then divide the particle count in the time period by the sample_vol to get particle concentration
%will need to get airspeed from the other aircraft data? - can get it from the turbulence probe, but
%CAS instrument should measure it too
%time looks to be 1 Hz data (data is counted up and reset every second)

%The file quotes a beam diameter as 0.24. However, this is wrong and this is the 
%actual beam area, i.e. A=0.24 mm^2 as confirmed by Bill Dawson via Tom LC through email.


last_CAS_missing=0;
last_CIP_missing=0;

%find out where the size range info starts (is two lines below "Laser Power: Low" string)
textline='';
while strcmp(textline,'"Laser Power: Low"')~=1
    textline=fgetl(fid);
end



[range_info_CAS]=textscan(fid,...
        '%d %d %d',1,'delimiter',','); %info on forward scatter sizes:-range,channel count,user table
    
channel_count = range_info_CAS{2};
[bins_CAS]=textscan(fid,...
        '%d %f',channel_count,'delimiter',',');   %data columns are threshold and size 
                %threshold is the voltage threshold and so not really needed
                %size data will then be in bins_CAS{2} and is diameter in microns.
                %Note have 30 bins for forwards scatter (30 for back).
                
[range_info_CAS_back]=textscan(fid,...
        '%d %d %d',1,'delimiter',',');   %info on backscatter sizes:- range, channel count, user table
    
channel_count = range_info_CAS_back{2};
[bins_CAS_back]=textscan(fid,...
        '%d %f',channel_count,'delimiter',',');   %data columns are threshold and size 
                %threshold is the voltage threshold and so not really needed
                %size data will then be in bins_CAS{2} and is diameter in microns.    
                                

%find out where the data starts (is labelled "DATA")
textline='';
while strcmp(textline,'"DATA"')~=1
    textline=fgetl(fid);
end


idat=0; %time index for the list of times read in
idat_CAS=0;
idat_CIP=0; %time indices for the two probes (both should be equal I think)

go=1;
while go==1 

idat=idat+1;

%%%%%%% Read in the time data - also tells us the probe number   %%%%%%
 [datevals{idat}]=textscan(fid,...
        '%d %s %d %s %d %d',1,'delimiter',','); %probe_num,datetime,millisecs,probe_datetime,probe_ms,interval_ms

    if length(datevals{idat}{1})==0
        last_CAS_missing=0;
        break %then have ran out of data, so stop
    end
    
    switch datevals{idat}{1} %datevals{idat}{1} is the instrument number - note sometimes CAS is first, sometims CIP
        case 0 %probe 0 (CAS)
            idat_CAS=idat_CAS+1;
            datevals_CAS{idat_CAS}=datevals{idat};
            
            %read in the CAS data using a function
            %note for future - could have passed e.g. Counts_CAS(idat_CAS) (i.e. non curly brackets) to 
            %the function - then would be able to access data as Counts_CAS{idat_CAS} instead of
            %Counts_CAS{idat_CAS}{1}
            [AD_values_CAS{idat_CAS},flag_vals_CAS{idat_CAS},Counts_CAS{idat_CAS},Bchannel_CAS{idat_CAS}...
                ,Psep_CAS{idat_CAS},last_CAS_missing,go] = read_CAS_probe_0(fid);
            
        
%The ambient data (pressure, temperature, airspeed, etc.) is only used by Tom LC's script from the CIP probe
%(He ignores the CAS AD_values)
            
        case 1 %probe 1 (CIP)
            idat_CIP=idat_CIP+1;
            datevals_CIP{idat_CIP}=datevals{idat};
            
            %read in the CIP data using a function
            [AD_values_CIP{idat_CIP},flag_vals_CIP{idat_CIP},Counts_CIP{idat_CIP}, ...
                last_CIP_missing,go] = read_CIP_probe_1(fid);
            

            if go==1
                press_CIP(idat_CIP) = AD_values_CIP{idat_CIP}{1}(5)*.2525553; %(static?)pressure in mb
                dynamic_CIP(idat_CIP)=AD_values_CIP{idat_CIP}{1}(4)*.0841851; 
                %dynamic pressure - can calculate airspeed from this and the static pressure
                lwc_CIP(idat_CIP)=AD_values_CIP{idat_CIP}{1}(6)*.002442;      %LWC (from hotwire?)
                rh_CIP(idat_CIP)=(AD_values_CIP{idat_CIP}{1}(9)-392)*.0795;   %RH (negative values?)
                temp_CIP(idat_CIP)=(AD_values_CIP{idat_CIP}{1}(14)-2048)*.024414062; %temperature (deg C)
            end
            
                %from Tom LC's script - note this is from IGOR where 0 is the first index for arrays                
                %so a1[n] = AD_values{idat_CAS}{1}(n+1)
                %             press[nele1]=a1[4]*.2525553
                %             dynamic[nele1]=a1[3]*.0841851
                %             lwc[nele1]=a1[5]*.002442
                %             rh[nele1]=(a1[8]-392)*.0795
                %             temp[nele1]=(a1[13]-2048)*.024414062
            
        otherwise
            disp('ERROR - unrecognised instrument')
            break
    end    
                       
end

fclose(fid);

%now pull out some of the data into normal arrays for ease of use

CAS_bins=bins_CAS{2};
CIP_bins = 25*[1:62]'; %25 micron resolution

clear CAS_counts CAS_time CAS_psep AD_CAS stats_CAS CAS_back
if exist('Counts_CAS')
    L_CAS=length(Counts_CAS);
    if last_CAS_missing==1
        L_CAS=L_CAS-1;
    end

    i=1;
    date_str=datevals_CAS{1}{2}{1}(2:11);
    CAS_counts(i,:)=Counts_CAS{i}{1};
    CAS_back(i,:)=Bchannel_CAS{i}{1};
    CAS_psep(i,:)=Psep_CAS{i}{1};
    time_str=datevals_CAS{i}{2}{1}(13:20);
    hours=str2num(time_str(1:2));
    mins=str2num(time_str(4:5));
    secs=str2num(time_str(7:8));
    milli=double(datevals_CAS{1}{3});
    CAS_time(i)=hours*3600+mins*60+secs+milli/1000; %calculate the time in seconds
    CAS_start_time = CAS_time(i);
    CAS_start_time_str = time_str;
    AD_CAS(i,:)=AD_values_CAS{i}{1};
    stats_CAS(i,:)=flag_vals_CAS{i}{1};
    CAS_bins_back = bins_CAS_back{2};
else
    L_CAS=0;
end



for i=2:L_CAS
    CAS_counts(i,:)=Counts_CAS{i}{1};
    CAS_back(i,:)=Bchannel_CAS{i}{1};    
    CAS_psep(i,1:length(Psep_CAS{i}{1}))=Psep_CAS{i}{1};
    time_str=datevals_CAS{i}{2}{1}(13:20);
    hours=str2num(time_str(1:2));
    mins=str2num(time_str(4:5));
    secs=str2num(time_str(7:8));
    milli=double(datevals_CAS{i}{3});
    CAS_time(i)=hours*3600+mins*60+secs+milli/1000; %calculate the time in seconds
    AD_CAS(i,:)=AD_values_CAS{i}{1};
    stats_CAS(i,:)=flag_vals_CAS{i}{1};
end

%%%% now for CIP
clear CIP_counts CIP_time stats_CIP

if ~exist('Counts_CIP')
    L_CIP=0;
else
    L_CIP=length(Counts_CIP);

    if last_CIP_missing==1
        L_CIP=L_CIP-1;
    end

    i=1;
    date_str=datevals_CIP{1}{2}{1}(2:11);
    CIP_counts(i,:)=Counts_CIP{i}{1};
    time_str=datevals_CIP{i}{2}{1}(13:20);
    hours=str2num(time_str(1:2));
    mins=str2num(time_str(4:5));
    secs=str2num(time_str(7:8));
    milli=double(datevals_CIP{1}{3});
    CIP_time(i)=hours*3600+mins*60+secs+milli/1000; %calculate the time in seconds
    CIP_start_time = CIP_time(i);
    CIP_start_time_str = time_str;
    stats_CIP(i,:)=flag_vals_CIP{i}{1};


end





for i=2:L_CIP
    CIP_counts(i,:)=Counts_CIP{i}{1};
    time_str=datevals_CIP{i}{2}{1}(13:20);
    hours=str2num(time_str(1:2));
    mins=str2num(time_str(4:5));
    secs=str2num(time_str(7:8));
    milli=double(datevals_CIP{i}{3});
    CIP_time(i)=hours*3600+mins*60+secs+milli/1000; %calculate the time in seconds
    stats_CIP(i,:)=flag_vals_CIP{i}{1};
end

%automatically sets the right flight data here for the airspeed calc based on the CAS
%fileName
    dat_flt_str = ['dat_flt' flight_no]; %make a string that is the name of the required flight data

if no_flight_data==1
    return  %uncomment this if want to skip the rest of the script (e.g. if don't have flight data)
end
    
%put the data in dat_flt array
    eval(['dat_flt=' dat_flt_str ';']);
    
    TAS_aircraft = interp1(dat_flt(:,1)/1e3,dat_flt(:,4),CIP_time);
    
    

R=6.8557e-2; %cal/gram/K
Cp=0.24; %cal/gram/K
Cv=0.171; %cal/gram/K %units don't matter as are always taking ratios
r=1;
Pstatic = press_CIP;
Pdynamic = dynamic_CIP;  %note the units of Pstatic and Pdynamic don't matter
%as the ratio is taken. But the housekeeping values are mulitplied by different
%numbers - presumably the formula requires the actual pressure values
%rather than the raw housekeeping (AD_values).

%press_CIP(idat_CIP) = AD_values_CIP{idat_CIP}{1}(5)*.2525553; %(static?)pressure in mb
%dynamic_CIP(idat_CIP)=AD_values_CIP{idat_CIP}{1}(4)*.0841851; 
%temp_CIP(idat_CIP)=(AD_values_CIP{idat_CIP}{1}(14)-2048)*.024414062; %temperature (deg C)
                
Mach = sqrt( 2*(Cv/R)*( (1+Pdynamic./Pstatic).^(R/Cp) - 1) );
Ta = (temp_CIP+273.16)./(1+r*Mach.^2*(Cp/Cv-1)/2);
%is the right temperature used here?
TAS = Mach * 20.06 .* sqrt(Ta);

if L_CAS>0
    % need to check that we are using the correct voltage here %
    lwc_cas=AD_CAS(:,28)*0.002442;  %pretty sure that this is the hotwire voltage
    lwcsv_cas=AD_CAS(:,29)*0.002442; % LWC slave voltage - think we don't need this
    lwc_voltage28 = interp1(CAS_time,lwc_cas,CIP_time,[],'extrap');
    lwc_voltage = interp1(CAS_time,lwcsv_cas,CIP_time,[],'extrap');
    LWC_CAS=lwc_calc_from_DMT(temp_CIP,press_CIP,TAS,lwc_voltage28);

else
    CAS_counts=[];
    LWC_CAS=[];
    CAS_time=[];   
    stats_CAS=[];
    AD_CAS28=[];
    CAS_psep=[]; %N.B. - need semicolon for 2-D arrays, don't for 1-D                    
    CAS_back=[]; %N.B. - need semicolon for 2-D arrays, don't for 1-D                                        
end

   

disp('Finished reading CAS data');


%%% CAS housekeeping data description (AD_CAS values), gain, offset
% 1     Pitot Pressure,.0841851,0
% 2     Static Pressure,.2525553,0
% 3     Ambient Temperature,.02442,-2047.502
% 4     Forward Heat Sink Temperature,0,0
% 5     Backward Heat Sink Temperature,0,0
% 6     Forward Block Temperature,0,0
% 7     Back Block Temperature,0,0
% 8     Photodiode 1,.002442,0
% 9     Photodiode 2,.002442,0
% 10	Photodiode 3,.002442,0
% 11	Photodiode 4,.002442,0
% 12	Qualifier TEC Temperature,0,0
% 13	Forward TEC Temperature,0,0
% 14	Backward TEC Temperature,0,0
% 15	Qualifier Heat Sink Temperature,0,0
% 16	Qualifier High Gain Voltage,.002442,0
% 17	Qualifier Mid Gain Voltage,.002442,0
% 18	Qualifier Low Gain Voltage,.002442,0
% 19	Forward High Gain Voltage,.002442,0
% 20	Forward Mid Gain Voltage,.002442,0
% 21	Forward Low Gain Voltage,.002442,0
% 22	Backward High Gain Voltage,.002442,0
% 23	Backward Mid Gain Voltage,.002442,0
% 24	Backward Low Gain Voltage,.002442,0
% 25	Internal Temperature,.4884,-558.9675
% 26	Spare,.002442,0
% 27	Spare,.002442,0
% 28	Hot Wire LWC,.002442,0
% 29	LWC Slave Voltage,.002442,0
% 30	Laser Current Monitor,.19536,0
% 31	Laser Power Monitor,.25,0
