init=100; %initial amount (k)
N=100; %max no. years
p=6; %annual interest 
R=700; %monthly total paid

[years,cost] = mortgage_calc_FUNC(init,N,p,R);


%Work out how much would have made if had invested the money in an ISA with
%the same rate as the mortgage
min_pay = 550;
d_pay = R-min_pay;
[years_min,cost_min] = mortgage_calc_FUNC(init,N,p,min_pay);
[savings_min] = saver_FUNC(years_min,d_pay,p);

%Now work out how much would make if invested the total for the remaining
%years in the quick pay off case
[savings] = saver_FUNC(years_min-years,R,p);

%tot_cost_min = R*12*years_min; %In both cases have spent 700 a month for years_min
ratio = savings ./ savings_min





