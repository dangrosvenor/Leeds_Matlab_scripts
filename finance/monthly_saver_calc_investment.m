cash=1700;
save_per_month=200; %per month amount

in=cash;

rate=7.75; %

for i=1:12
    cash = cash + cash*rate/100/12;
    cash = cash + save_per_month;
    in = in + save_per_month; %amount actually put in
end

cash_1 = cash; %money after 1 year of superISA
in_1 =in;
gain = (cash/in-1) *100;


%%%%% 5 year best case %%%%%%%%%%%%
rate=6; %
for i=1:12*4
    cash = cash + cash*rate/100/12;
    cash = cash + save_per_month;
    in = in + save_per_month; %amount actually put in
end

inv = 5*12*200 + 1500;
in5 = in + inv; %amount added in investment
b5 = inv*1.5;
gain_b5 = ( ( (b5 + cash) / in5 - 1 )*100 ) /5

%%%%% 5 year worst case %%%%%%%%%%%%
w5 = inv*1.17;
gain_w5 = ( ( (w5 + cash) / in5 - 1 )*100 ) /5


%%%%%%%%%%%3.5 years %%%%%%%%%%%%%%%
%%%%% 3.55 year best case %%%%%%%%%%%%
cash=cash_1;
in=in_1;
rate=6; %
for i=1:12*2.5
    cash = cash + cash*rate/100/12;
    cash = cash + save_per_month;
    in = in + save_per_month; %amount actually put in
end

inv = 3.5*12*200 + 1500;
in_3_5 = in + inv; %amount added in investment
b_3_5 = inv*1.35;
gain_b3_5 = ( ( (b_3_5 + cash) / in_3_5 - 1 )*100 ) /3.5

%%%%% 5 year worst case %%%%%%%%%%%%
w_3_5 = inv*1.05;
gain_w_3_5 = ( ( (w_3_5 + cash) / in_3_5 - 1 )*100 ) /3.5