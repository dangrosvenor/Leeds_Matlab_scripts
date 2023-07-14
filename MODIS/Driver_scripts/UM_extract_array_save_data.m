function [arr]=UM_extract_array_save_data(struc_in,var_str,irun_in,itime_in)

if strcmp(itime_in,'all_times')==1
    itime_in = 1:size(struc_in,2);
end

arr = NaN*ones([length(irun_in) length(itime_in)]);

for irun=1:length(irun_in)
for itime=1:length(itime_in)
    irun2=irun_in(irun);
    itime2=itime_in(itime);
%    struc_in{irun2,itime2}.P_mean{1}
    eval_str = ['val = struc_in{irun2,itime2}.' var_str '{1};'];
    eval(eval_str); 
    if length(val)==1
        arr(irun,itime)=val;
    end
end
end