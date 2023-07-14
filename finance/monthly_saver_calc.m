cash=25;
save_per_month=300; %per month amount
years=1;
in=cash;

rate=7; %

for i=1:12*years
    cash = cash + cash*rate/100/12;
    cash = cash + save_per_month;
    in = in + save_per_month; %amount actually put in
end

gain = (cash/in-1) *100 / years


