function EH=invHlogROOT(y,x,b,d,r)


	if y>=0
		EH=10.^(d.*y./r) + b.*d.*y./r - 1 - x; %minus x so that returns the value that should be zero for root finding
	else
        EH=-10.^-(d.*y./r) + b.*d.*y./r + 1 - x;
	end


