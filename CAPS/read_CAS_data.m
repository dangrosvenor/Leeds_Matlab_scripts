clear probe_num datetime millisecs probe_datetime probe_ms interval_ms

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
    
%i now determines how many lines to skip (comes after 'headerlines', 
%which skips i lines in file)
idat=0;

go=1;
while go==1

    idat=idat+1
%%%%%%% Now read in data for probe_num=0 %%%%%%


    [probe_num{idat},datetime{idat},millisecs{idat},probe_datetime{idat},probe_ms{idat},interval_ms{idat}]=textread(fileName,...
        '%d %s %d %s %d %d',1,'delimiter',',','headerlines',i); i=i+1;

    [AD_values{idat}]=textread(fileName,...
        '%f',31,'delimiter',',','headerlines',i); i=i+1;

    [RejectedDofCount,RejectedATCount,AverageTransit,FIFOfull,ResetFlag,...
        ADCoverflow,BackOverflow,OversizeRejects,EndRejects,SumOfParticles,SumOfTransit]=...
        textread(fileName,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',',','headerlines',i); i=i+1;

    %particle counts in channels 0-39. Note - is in 4 lines so need to increment i accordingly
    [Counts{idat}]=textread(fileName,...
        '%d',40,'delimiter',',','headerlines',i); i=i+4;

    [Bchannel{idat}]=textread(fileName,...
        '%d',20,'delimiter',',','headerlines',i); i=i+2; %in 2 lines

    [Psep{idat}]=textread(fileName,...
        '%d',64,'delimiter',',','headerlines',i); i=i+1;

%%%%%%  %now for probe_num=1 %%%%

    [probe_num_1{idat},datetime_1{idat},millisecs_1{idat},probe_datetime_1{idat},probe_ms_1{idat},interval_ms_1{idat}]=textread(fileName,...
        '%d %s %d %s %d %d',1,'delimiter',',','headerlines',i); i=i+1;

    [AD_values_1{idat}]=textread(fileName,...
        '%f',31,'delimiter',',','headerlines',i); i=i+1;

    [RejectedDofCount,RejectedATCount,AverageTransit,FIFOfull,ResetFlag,...
        ADCoverflow,BackOverflow,OversizeRejects,EndRejects,SumOfParticles,SumOfTransit]=...
        textread(fileName,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',',','headerlines',i); i=i+1;   

    %particle counts in channels 0-39. Note - is in 4 lines so need to increment i accordingly
    [Counts_1{idat}]=textread(fileName,...
        '%d',40,'delimiter',',','headerlines',i); i=i+4;
    
    %22 integers on one line - Bchannel? Psep?
    [Unknown{idat}]=textread(fileName,...
        '%d',22,'delimiter',',','headerlines',i); i=i+1; %                       


end

%[min,sec,hPa,gpm,T,RH,Dewp]=textread(fileName,...
%    '%d %s %d %s %d %d %f %*[^\n]','headerlines',42);


