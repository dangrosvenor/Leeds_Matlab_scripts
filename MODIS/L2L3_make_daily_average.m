function L2L3_make_daily_average(load_filename,save_filename,var_str,ihtot)
%load the data
eval(['load(load_filename,''' var_str ''');']);

var_str2 = [var_str '.timeseries3'];
%screen the data
var_dat = eval(var_str2);
var_dat(ihtot)=NaN;

if length(strfind(var_str,'Maximum'))>0
    eval([var_str2  ' = max(var_dat,[],3);']);
elseif length(strfind(var_str,'Minimum'))>0
    eval([var_str2 ' = min(var_dat,[],3);']);
else
    eval([var_str2 ' = meanNoNan(var_dat,3);']);
end

%save to the specified file
if exist(save_filename)==2
    iappend=1;
else
    iappend=0;
end


if iappend==1
    eval(['save(save_filename,''' var_str ''',''-APPEND'',''-V7.3'');']);
else
    eval(['save(save_filename,''' var_str ''',''-V7.3'');']);
end


