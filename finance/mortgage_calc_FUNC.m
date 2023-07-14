function [years,cost_tot] = mortgage_calc_FUNC(init,N,p,R)

%x(1)=100; %initial amount (k)
%N=30; %no. years
%p=6; %annual interest 
%R=700; %monthly repayment

x(1)=init; %initial amount (k)
cost(1)=0;
f=100;

ann_pay=R*12/1000;
for i=2:12*f*(N+1)
    x(i)=( x(i-1) - ann_pay/12/f )*(1+p/100/12/f);
    cost(i) = cost(i-1) + ann_pay/12;
end

b=find(x<0);
years=(b(1)-1)/12/f;
cost_tot = cost(b(1))*1000;
%20 years to pay off �80k at 550 per month



