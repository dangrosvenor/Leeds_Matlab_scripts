%run the script to read all BAS flight data first
% - is called read_BAS_Feb2010_flights
%gain stage cross-overs are at 1.5 and 10.2 microns for the CAS.

if ~exist('call_from_cas_hotwire_ratio_all_flights') | call_from_cas_hotwire_ratio_all_flights==0 
    %don't set own filedir
    
if ~exist('comp')
    comp='uni';
    comp='laptop';
    comp='UWchallenger';
end

switch comp
    case 'uni'        
        filedir_cas = 'Y:\BAS_flights\';
    case 'UWchallenger'
        filedir_cas = '/home/disk/eos1/d.grosvenor/Ant_flights/';
end

%filedir = [filedir_cas 'flight99_6thFeb2010\'];
%filedir=[filedir_cas 'CAPS_cals_02-02-10\]';
filedir=[filedir_cas 'flight104\'];
%filedir=[filedir_cas 'caps_cals_2010_2\CAPS_cals_08-02-10\'];
%filedir=[filedir_cas 'caps_cals_2010_2\CAPS_cals_10-02-10\'];
%filedir=[filedir_cas 'caps_cals_2010_2\CAPS_cals_14-02-10\'];
%filedir=[filedir_cas 'caps_cals_2010_2\CAPS_cals_02-02-10\'];
%filedir=[filedir_cas 'caps_cals_2010_2\CAPS_tests_31-01-10\'];
%filedir=[filedir_cas '29thJune2010_CAS_comparisons\bas_cas_29-06-10\'];
%N.B. files 6,7,8 and 9 correspond to the calibrations of 15,20,30 and 40 microns

%filedir=[filedir_cas '29thJuly2010_BAS_chamber_comparisons\BAS CAS\'];

%think that the one below contains all the requried ones
%filedir=[filedir_cas '29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\cas_probes\bas_cas\'];

switch comp
    case 'UWchallenger'
        ifind=strfind(filedir,'\');
        filedir(ifind)='/';
end

else %otherwise set own filedir and clear the flag for any future runs
    clear call_from_cas_hotwire_ratio_all_flights
end

ijust_plot=0; %flag to just plot the data (already read in)
ijust_read=1; %flag to just read and not plot
iselect_cas_file=0; %flag to say that want to select a cas file - will ask to select a file. Otherwise reads all

clear files dir

if ijust_plot==1
    files.name='one_file.txt';
else
    files=dir(filedir);
    CAS_time_all=[];
    CAS_counts_all=[];
    CIP_time_all=[];
    CIP_counts_all=[];
    TAS_all=[];
    LWC_CAS_all=[];
    stats_CAS_all=[];
    stats_CIP_all=[]; 
    hotwire_raw_all=[];
    CAS_psep_all=[];
    CAS_back_all=[];    
    clear CAS_total_number

end


cas_file=0; %default - read all CAS files
icount=0;
if iselect_cas_file==1
    for ifile=1:length(files)        
        %    if findstr(files(ifile).name,'020610') & findstr(files(ifile).name,'.txt')
        if length(files(ifile).name)>3
            if strfind(files(ifile).name(end-3:end),'.txt') & length(strfind(files(ifile).name,'_s1.txt'))==0 ...
                    & length(strfind(files(ifile).name,'ncd'))==0
                icount=icount+1;
                fileName=[filedir files(ifile).name];
                disp([num2str(icount) ')' fileName]);
            end
        end

    end
    
cas_file = input('Select file : ') ;

end




icas_count=0;
for ifile=1:length(files)
%    if findstr(files(ifile).name,'020610') & findstr(files(ifile).name,'.txt')

if length(files(ifile).name)>3
    text_file = strfind(files(ifile).name(end-3:end),'.txt');
else
    text_file = 0;
end

    if text_file & length(strfind(files(ifile).name,'_s1.txt'))==0 ...
            & length(strfind(files(ifile).name,'cpc'))==0 & length(strfind(files(ifile).name,'ncd'))==0
        icas_count=icas_count+1;
        
        if icas_count==cas_file | cas_file==0 | ijust_plot==1
            
            if ijust_plot==0
                fileName=[filedir files(ifile).name];
                disp(fileName);
                read_CAS_data_textscan2;

                %combine the data from each file into single arrays
                CAS_time_all = [CAS_time_all CAS_time];
                CAS_counts_all = [CAS_counts_all; CAS_counts];
                
                if exist('CIP_time') %CIP specific data
                    CIP_time_all = [CIP_time_all CIP_time];
                    CIP_counts_all = [CIP_counts_all; CIP_counts];
                    stats_CIP_all = [stats_CIP_all; stats_CIP]; %N.B. - need semicolon for 2-D arrays, don't for 1-D
                end
                
                stats_CAS_all = [stats_CAS_all; stats_CAS]; %N.B. - need semicolon for 2-D arrays, don't for 1-D

                
                CAS_psep_all = [CAS_psep_all; CAS_psep]; %N.B. - need semicolon for 2-D arrays, don't for 1-D                    
                CAS_back_all = [CAS_back_all; CAS_back]; %N.B. - need semicolon for 2-D arrays, don't for 1-D                                    
                
                if no_flight_data==0
                    TAS_all = [TAS_all TAS];
                    LWC_CAS_all = [LWC_CAS_all LWC_CAS];
                    hotwire_raw_all = [hotwire_raw_all; lwc_cas]; %N.B. - need semicolon for 2-D arrays, don't for 1-D                                                    
                end
                
            end
            
            if ijust_read==0 | ijust_plot==1

                itype='timh';
                multisaveplot;

                itype='prof'; %(or timeseries)
                itimser = 'CAS plots';

                time_graph = {'Total number CAS'};
                multisaveplot;

                time_graph = {'Total number CIP'};
                multisaveplot;

                time_graph = {'LWC_dist_CAS','LWC_hotwire','LWC_dist_CAS_cutoff'};
                multisaveplot;

%                time_graph = {'Last bin number CAS'};
                time_graph = {'Number in latter bins for CAS'};
                multisaveplot;

%                time_graph = {'Second bin number CIP'};
                time_graph = {'Second bin dN/dlogD CIP'};
                multisaveplot;

                itype='scatter'; %(or timeseries)
                
                plot_type='CAS LWC vs hotwire LWC';
                multisaveplot;
                
                plot_type='CAS LWC vs conc.';
                multisaveplot;

                plot_type='Hotwire LWC vs conc.';
                multisaveplot;

                plot_type='Hotwire LWC vs mean';
                multisaveplot;

                plot_type='CAS LWC vs mean'; 
                multisaveplot;
                
                plot_type='CAS LWC vs height';
                multisaveplot;
                
                plot_type='Hotwire LWC vs height';
                multisaveplot;
                
                plot_type='Reject diags';
                multisaveplot;
                
            end

        end
    end
    
end

try %in case the files aren't named so that flight_no is retrieved
    eval(['dat_flt = dat_flt' flight_no ';']);  %put the data for the required flight here
    eval(['time_flt = time_flt' flight_no ';']);  %put the data for the required flight here
    eval(['qv_flt_fp = qv_flt' flight_no '_fp;']);  %put the data for the required flight here
    eval(['qv_flt_humi = qv_flt' flight_no '_humi;']);  %put the data for the required flight here
    
catch
    disp('*** WARNING *** flight data has not been successfully put into dat_flt');    
end
