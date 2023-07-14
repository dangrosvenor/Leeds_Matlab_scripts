function [AD_vals_CAS,flag_vals_CAS,Counts_CAS,Bchannel_CAS,Psep_CAS,last_CAS_missing,go] ...
    = read_CAS_probe_0(fid)

      
    
    [AD_vals_CAS]=textscan(fid,...   %AD_values are things like instrument voltages, etc. - housekeeping, not needed
        '%f',31,'delimiter',','); 
    
%    [flag_vals_CAS]=... %[RejectedDofCount,RejectedATCount,AverageTransit,FIFOfull,ResetFlag,ADCoverflow,BackOverflow,OversizeRejects,EndRejects,SumOfParticles,SumOfTransit]
%        textscan(fid,'%d %d %d %d %d %d %d %d %d %d %d',1,'delimiter',','); 
    
    [flag_vals_CAS]=... %[RejectedDofCount,RejectedATCount,AverageTransit,FIFOfull,ResetFlag,ADCoverflow,BackOverflow,OversizeRejects,EndRejects,SumOfParticles,SumOfTransit]
        textscan(fid,'%f',11,'delimiter',','); 

    %particle counts in channels 0-29. Note have 30 forward and 30 backscatter channels
    [Counts_CAS]=textscan(fid,...
        '%f',30,'delimiter',','); 
    
    if size(Counts_CAS{1},1)==0 %when runs out of data to read it returns a zero length variable
        last_CAS_missing=1
        go=0; %exit the loop
    else
        last_CAS_missing=0;
        go=1;
    end

    [Bchannel_CAS]=textscan(fid,...   %back scatter counts
        '%f',30,'delimiter',',');  %in 2 lines

    [Psep_CAS]=textscan(fid,...       %Histogram of interarrival times? Could be useful to see if have shattering if so
        '%f',64,'delimiter',','); 