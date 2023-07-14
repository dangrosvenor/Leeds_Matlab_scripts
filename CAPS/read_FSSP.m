%read in the FSSP output (.min files) from the Manchester probe

iselect_cas_file=0; %flag to say that want to select a file - will ask to select a file
ijust_plot=0;
    
%

% file_choose='29th July BAS chamber comparisons, 50um bead calibraion';
% file_choose='29th July BAS chamber comparisons, 60um bead calibraion';
 file_choose='29th July BAS chamber comparisons, experiment 1';
% file_choose='29th July BAS chamber comparisons, experiment 2';
% file_choose='29th July BAS chamber comparisons, experiment 3';
% file_choose='29th July BAS chamber comparisons, experiment 4';
% file_choose='29th July BAS chamber comparisons, experiment 5';

switch file_choose
    case '29th July BAS chamber comparisons, experiment 1'
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_expt_1\';
        fileName=[filedir '0729100.min'];       
%        fileName=[filedir '07291252.min'];   
    case '29th July BAS chamber comparisons, experiment 2'
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_expt_2\';        
        %Paul's .mat file for expt 2 for comparison
        %comparison with this comfirms that this script reads
        %in ok and that the numbers are concentrations and bins
        %are in microns
        fileName_Paul = [filedir 'FSSP-exp2.mat'];
        fileName=[filedir '0729123.min'];    
    case '29th July BAS chamber comparisons, 50um bead calibraion'
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_50um_beads\';                
%%%  NO DATA     fileName=[filedir '0729345.min'];  

%         fileName=[filedir '0729346.min'];  %one second of data
%  but has a defined peak

%%%  NO DATA     fileName=[filedir '0729347.min'];  


%        fileName=[filedir '0729350.min'];   %most data here       
        
    case '29th July BAS chamber comparisons, 60um bead calibraion' 
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_60um_beads\';                
        fileName=[filedir '0729357.min'];  
        
    case '29th July BAS chamber comparisons, experiment 3' 
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_expt_3\';
        
    case '29th July BAS chamber comparisons, experiment 4' 
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_expt_4\';        
    
    case '29th July BAS chamber comparisons, experiment 5' 
        filedir='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\fssp\100729_expt_5\';        
        

end



clear files dir

if ijust_plot==1
    files.name='one_file.txt';
else
    files=dir(filedir);
    data_FSSP=[];
    bins_FSSP=[];   
    time_FSSP=[];
end


cas_file=0; %default - read all CAS files
icount=0;
if iselect_cas_file==1
    for ifile=1:length(files)        
        %    if findstr(files(ifile).name,'020610') & findstr(files(ifile).name,'.txt')
        if length(files(ifile).name)>3
            if strfind(files(ifile).name(end-3:end),'.min')
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
    text_file = strfind(files(ifile).name(end-3:end),'.min');
else
    text_file = 0;
end

    if text_file 
        icas_count=icas_count+1;
        
        if icas_count==cas_file | cas_file==0
            

                fileName=[filedir files(ifile).name];
                disp(fileName);
                [bins_FSSP,data_FSSP_file,time_FSSP_file]=read_FSSP_file(fileName);

                %combine the data from each file into single arrays
                data_FSSP = [data_FSSP data_FSSP_file];
                time_FSSP = [time_FSSP time_FSSP_file];                
                

        end
    end
    
end

%sort the data in time order (as files may not have been read in time order)
[time_FSSP, I] = sort(time_FSSP);
data_FSSP = data_FSSP(:,I);





