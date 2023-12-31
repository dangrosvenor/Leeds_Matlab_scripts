function [profit preturn break_even]= share_profit(buy,sell,n)
%function profit = share_profit(buy,sell,n) - n=number of shares, amounts in �s
%profit = n*sell - n*buy - 0.5/100*n*buy - 30;  %30 is the comission (15*2)

cost = n*buy + 0.5/100*n*buy + 30;
profit = n*sell - cost;
preturn = profit/cost *100;
break_even = cost/n; %since for break-even cost = n*sell (profit=0)


%LLOY - 11,278 @ 0.6022 average
%BARC - 1737   @ 1.82 average


exe=0; %don't set to one as causes infinite loop! 
%highlight and execute the section
if exe==1
%sell_LLOY=0.766;  %summer 2010
%sell_BARC=3.4; %summer 2010

sell_LLOY=0.7616;  %20th Sep 2010
sell_BARC=3.06; %20th Sep 2010

sell_LLOY=0.7275; %12th Oct 2010
sell_BARC=2.935;  %12th Oct 2010

sell_LLOY=0.7090; %
sell_BARC=2.8195;  %

sell_BARC=3.368;  %sold price 17th Feb, 2011


[profit_LLOY preturn_LLOY break_even_LLOY]= share_profit(0.6022,sell_LLOY,11278)
[profit_BARC preturn_BARC break_even_BARC]= share_profit(1.82,sell_BARC,1737)
days=datenum(date)-datenum('01-Oct-2008');
percent_LLOY_per_year = preturn_LLOY/(days/365)
percent_BARC_per_year = preturn_BARC/(days/365)
total_profit = profit_BARC+profit_LLOY
overall=100*(profit_BARC+profit_LLOY)/(1737*1.82 + 11278*0.6022)  %overall percentage profit
overall_per_year = overall/(days/365)
end
%2900 now or 3740 later if get to 4.00? 840 difference or almost 30%. Dividend consideration?
%Price of 4.15 brings 1112, which is 38.5% more.