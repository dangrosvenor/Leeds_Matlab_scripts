clear datevals_CAS AD_vals_CAS flag_vals_CAS Counts_CAS Bchannel_CAS Psep_CAS datevals_CIP AD_vals_CIP flag_vals_CIP Counts_CIP

fileName='Y:\BAS_flights\flight99_6thFeb2010\flight99_020610_150001.txt';
fid=fopen(fileName,'rt');

%find out where the data starts (is labelled "DATA")
go=1;
i=0;
while go==1
    i=i+1;
    textline=fgetl(fid);
    if strcmp(textline,'"DATA"')==1
        go=0;
    end
end
    


%[probe_num{idat},datetime{idat},millisecs{idat},probe_datetime{idat},probe_ms{idat},interval_ms{idat}]=fscanf(fid,...
%        '%d %s %d %s %d %d',1);
    
%[dat]=fscanf(fid,'%d,%s,%d,%s,%d,%d',1); %need to explicitly tell it what to read in






%i now determines how many lines to skip (comes after 'headerlines', 
%which skips i lines in file)
idat=0;

go=1;
while go==1 

    idat=idat+1
%%%%%%% Read in data for probe_num=0 --- CAS instrument ---   %%%%%%


%    [probe_num{idat},datetime{idat},millisecs{idat},probe_datetime{idat},probe_ms{idat},interval_ms{idat}]=textscan(fileName,...
%        '%d %s %d %s %d %d',1,'delimiter',','); i=i+1;
    
    [datevals_CAS{idat}]=textscan(fid,...
        '%d %s %d %s %d %d',1,'delimiter',',');   %see above for a list of what these are

%    [AD_values{idat}]=textscan(fileName,...
%        '%f',31,'delimiter',','); i=i+1;
    
    [AD_vals_CAS{idat}]=textscan(fid,...   %AD_values are things like instrument voltages, etc. - housekeeping, not needed
        '%f',31,'delimiter',','); 

%    [RejectedDofCount,RejectedATCount,AverageTransit,FIFOfull,ResetFlag,...
%        ADCoverflow,BackOverflow,OversizeRejects,EndRejects,SumOfParticles,SumOfTransit]=...
%        textscan(fileName,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',','); i=i+1;
    
    [flag_vals_CAS{idat}]=...   %flags -see above
        textscan(fid,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',','); 

    %particle counts in channels 0-39. Note - is in 4 lines so need to increment i accordingly
    [Counts_CAS{idat}]=textscan(fid,...
        '%d',40,'delimiter',','); 

    [Bchannel_CAS{idat}]=textscan(fid,...   %back scatter counts - probably not needed
        '%d',20,'delimiter',',');  %in 2 lines

    [Psep_CAS{idat}]=textscan(fid,...       %Jonny thinks this is a histogram of interarrival times - useful to see if have shattering
        '%d',64,'delimiter',','); 

%%%%%%  %now for probe_num=1 --- CIP instrument ---  %%%%

    [datevals_CIP{idat}]=textscan(fid,...
        '%d %s %d %s %d %d',1,'delimiter',','); 

    [AD_values_CIP{idat}]=textscan(fid,...
        '%f',31,'delimiter',','); 

    [flag_vals_CIP{idat}]=...
        textscan(fid,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',','); 

    %particle counts in channels 0-61. Note - this is for the CIP
    [Counts_CIP{idat}]=textscan(fid,...
        '%d',62,'delimiter',','); 
    
    if size(Counts_CIP{idat}{1},1)==0 %when runs out of data to read it returns a zero length variable
        go=0; %exit the loop
    end
                       
end

%[min,sec,hPa,gpm,T,RH,Dewp]=textscan(fid,...
%    '%d %s %d %s %d %d %f %*[^\n]','headerlines',42);


