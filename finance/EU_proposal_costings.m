pay=[27466 28290 29138 30013 30912 31840]; %my salary for years following May 2008
pay=[33309 34324 35367 36445 37552 38694]; %my salary (total cost) for years following May 2008
year=[2008:2012];
xr=1.33929; %exchange rate for euros
%xr=1.44; if xr goes up we get fewer pounds per euro so is bad...
travel = 4000;
consumables = 4000; %consumables
start_year=2010;
months=18;

pay2=70000;
pay2=57000;

months_tom=4;
pay_tom=months_tom/12 * pay2;


years=months/12;

iyear=findheight(year,start_year);

whole_years = floor(years);
rem_years = rem(years,whole_years);
totpay=sum( pay(iyear:iyear+whole_years-1 ) ) + rem_years*pay(iyear+whole_years);

%sum(pay(2:3)*xr)

totA=( (totpay+pay_tom)*xr + travel + consumables )*1.6/1000  / 0.75; %1.6 is for the 60 % charge from the university.
                                        % divide by 0.75 is so that when take 75 % get the total needed (EU only gives 75%)

totB=( (totpay+pay_tom)*xr + travel + consumables )*1.6/1000 + 6*xr

%tot=totA + 6*xr  % add on £6k management/audit fees
 
tot=totB/0.75;


travel= 3000;
consumables =3000;
tot2 = 180000; %tot in euros
tot2 = 80000;
p= ( tot2*0.75/1.6 - travel - consumables )/xr/1000;   %inverse calculation to get the pay out of a given total
                                        
                                        