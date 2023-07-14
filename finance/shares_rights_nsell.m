function [nsell,nnew,ratio,av_price,av_price_new]=shares_rights_nsell(N1,current_price,extra_cash_invest)

ratio_available=1.34; %number of shares available per share held


N2= (extra_cash_invest + current_price*N1 ) / (ratio_available*0.37 + current_price );

nsell = N1-N2;

nnew = N2 + N2*ratio_available;

ratio = nnew / (N1+N1*ratio_available); %ratio of new shares to total amount possible

av_price = (0.91 + ratio_available*0.37 ) / (1+ratio_available); %average price is independant of how many are bought

av_price_new = ( N1*0.91 + extra_cash_invest ) / nnew ; %average price taking into account any profit made from the inital sale
                                                        %total money spent divided by final no. shares

