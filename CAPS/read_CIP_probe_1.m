function [AD_values_CIP,flag_vals_CIP,Counts_CIP,last_CIP_missing,go] = read_CIP_probe_1(fid)

%housekeeping instrument information
[AD_values_CIP]=textscan(fid,...
        '%f',31,'delimiter',','); 

[flag_vals_CIP]=... %[RejectedDofCount,RejectedATCount,AverageTransit,FIFOfull,ResetFlag,ADCoverflow,BackOverflow,OversizeRejects,EndRejects,SumOfParticles,SumOfTransit]
        textscan(fid,'%f',11,'delimiter',','); 

    %particle counts in channels 0-61. Note - this is for the CIP
[Counts_CIP]=textscan(fid,...
        '%f',62,'delimiter',','); 
    
    if size(Counts_CIP{1},1)==0 %when runs out of data to read it returns a zero length variable
        last_CIP_missing=1;
        go=0; %exit the loop
    else
        last_CIP_missing=0;
        go=1;
    end