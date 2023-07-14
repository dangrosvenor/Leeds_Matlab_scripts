function periods = consequetive_periods_find_func(i,if_str,vector01)
%look for consecutive periods
if length(i)==0
    periods=[];
    return
end
if nargin==1
    if_str = 'if go==1';
    go=1;
end
idiff=diff(i);
icon = find(idiff~=1);
%account the last part that runs from icon(end) to the end of the array
if length(icon)>0
    Lcon=length(icon)+1;
else
    %if there were no icon values (i.e. the whole period was continuous, or was just one
    %value)
    Lcon=1;
end
icon(Lcon) = length(idiff)+1;


istart=1;
iperiods=1;
for iperiods2=1:length(icon)
    ivals = i(istart:icon(iperiods2));  
    eval_str=[if_str ';'...
                 'periods(iperiods).p = ivals'';'...
                 'iperiods=iperiods+1;'...
              'end'];
    eval(eval_str);
    istart=icon(iperiods2)+1;
end
% if length(iperiods2)==0 %if there were no icon values (i.e. the whole period was continuous)
%     %or if there was only one value in i (then would be no icon diffs)
%     periods(1).p=i(istart):i(end);
% end
% if ~exist('periods')
%     Lperiods = 0;
% else
%     Lperiods = length(periods);
% end
% if istart<=i(end)
%     periods(Lperiods+1).p = i(istart):i(end);
% end