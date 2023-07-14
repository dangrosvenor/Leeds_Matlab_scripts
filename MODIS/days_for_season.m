switch time_mean_str
    case 'DJF'
        days_required_for_mean = [336:366 1:60];
    case 'MAM'
        days_required_for_mean = [61:152];
    case 'JJA'
        days_required_for_mean = [153:244];
    case 'SON'
        days_required_for_mean = [245:335];
end

%don't need to do for 'ALL' as this is done by time_inds_modisL3_timeseries3