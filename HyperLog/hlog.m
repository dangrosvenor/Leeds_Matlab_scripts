function h=hlog(x)

guess=0;
b=350;
d=1;
r=1;

for i=1:length(x)
    h(i)=fzero(@invHlogROOT,guess,[],x(i),b,d,r); %2nd arg=first guess, 4th=b, 5th=b, 6th=r
end