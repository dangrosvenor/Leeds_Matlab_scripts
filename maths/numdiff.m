function diff=numdiff(func,x,y,dx,xy)

switch xy
case 'x'
	diff=( feval(func,x+dx/2,y) - feval(func,x-dx/2,y) ) ./ dx;
case 'y'
	diff=( feval(func,x,y+dx/2) - feval(func,x,y-dx/2) ) ./ dx;   
otherwise
    fprintf(1,'\n*** need to enter ''x'' or ''y'' for the choice of variable (1 or 2 as entered) to differentiate over ***');
end