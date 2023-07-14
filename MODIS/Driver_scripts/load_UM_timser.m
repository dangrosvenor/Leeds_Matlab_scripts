function [dat_UM_timser,flux_UM_timser] = load_UM_timser(file_tim,append_str)

%Default outputs
dat_UM_timser.timeseries_UM=[];
flux_UM_timser.timeseries_UM=[];

switch append_str
    case '.txt'  %If we have a .txt created by Python on Postproc, etc. then read in here
        file_tim = [file_tim append_str];
        fid = fopen(file_tim,'rt');
        if fid==-1
            error(['DPG - file does not exist (' file_tim ')']);
        end
        header = textscan(fid,'%s %s %s %s %s %s %s %f',1);
        sf = header{8};
        
%        first_row = textscan(fgetl(fid),'%f');
%        nvars = length(first_row{1});

% Read the file in line-by-line for unknown number of entries on each line
        dat_all = textscan(fid,'%[^\n]');  %  %[^...] - reads characters not matching characters between the
%                  brackets until first matching character.
% So, it reads everything up to a new line and then starts a new cell array
% entry
        ifluxes=0;
        %Loop over all lines
        for i=1:length(dat_all{1})
            dat_line = textscan(dat_all{1}{i},'%f');
            dat(i,:) = dat_line{1}(1:2);
            if length(dat_line{1})>2
                ifluxes=1;
                dat_flux(i,:) = dat_line{1}(3:end);
            end
        end
            
            
            
%        dat = fscanf(fid,'%f %f',[2 Inf]);       
        %could also do this:-
%        dat=textscan(fid,'%s%s'); 

        fclose(fid);
        
        dat_UM_timser.timeseries_UM = dat(:,2)/sf;
        dat_UM_timser.time_UM = datenum('01-Jan-1970') + dat(:,1)/24;
        
        if ifluxes==1
            flux_UM_timser.timeseries_UM = dat_flux(:,:)/sf;
            flux_UM_timser.time_UM = datenum('01-Jan-1970') + dat(:,1)/24;
        end

        
    otherwise
        dat_UM_timser = load(file_tim);
        
end