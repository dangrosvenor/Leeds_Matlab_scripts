switch weight
    case 'own'
        iweight=i;
    case 'choose'
        %use iweight
end
w = SW_up_TOA(iweight).timeseries_UM(:)';
switch weight
    case 'none'
        w(:)=1;
    case 'SW_down_TOA'
        w = SW_down_TOA.timeseries_UM(1:end-1);
    case 'Time range'
        w(:)=0;
        for it=1:length(time_range_weight)
            itfind_temp = find(SW_up_TOA(1).time_UM >=time_range_weight{it}(1) & SW_up_TOA(1).time_UM <= time_range_weight{it}(2));
            if it>1
                itfind = cat(2,itfind,itfind_temp);
            else
                itfind = itfind_temp;
            end
        end
        w(itfind)=1;
end