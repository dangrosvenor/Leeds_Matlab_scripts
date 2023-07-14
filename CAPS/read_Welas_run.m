%read Welas instrument


file_case = '29th July BAS tests expt 2';

switch file_case
    case 'test from Chris';
        filedir='C:\Documents and Settings\dan\My Documents\logbook\Antarctica\Flights and instruments_Feb2010\Welas\';
        filename = [filedir '100713_-15v_2_0000.array.pal'];
    case '29th July BAS tests expt 1';
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\';
        filename = [filedir '100729_expt_1_0000.array.pal'];       
    case '29th July BAS tests expt 2';
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\';
        filename = [filedir '100729_expt_2_0000.array.pal'];        
    case '29th July BAS tests expt 3';
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\';
        filename = [filedir '100729_expt_3_0000.array.pal'];                
    case '29th July BAS tests expt 4';
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\';
        filename = [filedir '100729_expt_4_0000.array.pal'];              
    case '29th July BAS tests expt 5';
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\';
        filename = [filedir '100729_expt_5_0000.array.pal'];
end
        
welas_dat= readwelas(filename);

%in case don't have the calibration file - make the sizes up! (just for testing)
i_create_sizes=0;
if i_create_sizes==1
    n=4096;
    minD=2;
    maxD=105;
    dlogD=(log10(maxD)-log10(minD))/(n-1);
    logD=[log10(minD):dlogD:log10(maxD)];
    welas_dat.size = repmat(10.^logD,[100 1])';    
end


start_date_time = datestr(welas_dat.startTime);
%find the datenum for the start of the day so can work from there
start_day_datenum = datenum([start_date_time(1:12) '00:00:00']);

%convert times to run from the start of the day in question (in seconds)
%datenums are in hours from 0 UTC on year zero
welas_dat.time_of_day = (welas_dat.time-start_day_datenum)*24*3600;
welas_dat.conc1(:,length(welas_dat.time_of_day)+1:end)=[];


disp('done read_Welas_run')