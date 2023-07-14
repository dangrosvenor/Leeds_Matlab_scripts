x(1)=100; %initial amount (k)
N=30; %no. years
p=6; %annual interest 
R=700; %monthly repayment


ann_pay=R*12/1000;
for i=2:N+1
    x(i)=( x(i-1) - ann_pay )*(1+p/100);
end

b=find(x<0);
years=b(1)-1
%20 years to pay off �80k at 550 per month

%Work out how much would have made if had invested the money in an ISA with
%the same rate as the mortgage

