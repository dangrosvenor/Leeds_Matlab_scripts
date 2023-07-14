function EH=invHlog(x)

b=333;
d=3;
r=1;

if x>=0
	EH=10.^(d.*x./r) + b.*d.*x./r - 1;
else
    EH=-10.^-(d.*x./r) + b.*d.*x./r + 1;
end

