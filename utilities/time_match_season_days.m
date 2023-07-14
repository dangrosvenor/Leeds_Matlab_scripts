function iseason = time_match_season_days(days,season)
%function iseason = time_match_season_days(days,season)

switch season
    case 'DJF'
        idays = [days_of_month(12) days_of_month(1) days_of_month(2)];
    case 'MAM'
        idays = [days_of_month(3) days_of_month(4) days_of_month(5)];
    case 'JJA'
        idays = [days_of_month(6) days_of_month(7) days_of_month(8)];
    case 'SON'
        idays = [days_of_month(9) days_of_month(10) days_of_month(11)];
end


icount=0;
for i=1:length(idays)
    id = find(days==idays(i));
    if length(id)==1
        icount=icount+1;
        iseason(icount) = id;
    elseif length(id)>1
        error('More than one day found!')
    end
end
