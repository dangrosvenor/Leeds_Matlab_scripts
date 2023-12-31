%working in billions

n_raised=1.407+0.169; %1407 million shares raised at 282p and 169 million at 296p, 25th June, 2008.

%prices on 28th April, 2009
market_cap = 18.90391;
share_price = 2.2825;
n_current = market_cap/share_price; %current number of shares = 8.3 billion

n_previous = n_current-n_raised; %before the 25th June input of shares = 6.7 billion

market_cap_previous = n_previous*2.82; %=18.9 billion around 25th June -- exactly the same as that on 28th April, 2009!

%so there is scope for growth provided market conditions can return to those before the crash
%since share price was approaching 500p at the start of 2008
%assuming the same market capacity (BUT has there been any more share diution since then - what was market cap then?)
%number of shares look to have been approx similar 
%eps=[59.3 68.9 71.9 54.4 51.2 42.3]/100; = earnings (net profit) per share for 2008, 2007, 2006, 2005 & 2004 - from Yahoo finance
%net_prof=[4.382 4.417 4.571 3.447 3.268 2.744]; = net profit for those years
%so no. shares = net_prof./eps = [7.3895   6.4107    6.3574    6.3364    6.3828    6.4870] billion shares - approx constant except for 2008 when had
%shares issued - 2008 figure is probably some kind of time average as there were approx 8.3 billion after June 2008 - in fact 0.5*8.3 + 0.5*6.4 = 7.35
%(since June is approx halfway through the year)

aim_price = n_previous*4.75 / n_current   %=�3.84


%2007 - sp max of 8.90 low of 4.90
%2008 - sp max of 500 low of 1.50

aim_2008 = 5 * 6.4 / n_current ; %using start of 2008. gives �3.86
aim_2007 = 8.9 * 6.36 / n_current ; %using start of 2008. gives �6.83