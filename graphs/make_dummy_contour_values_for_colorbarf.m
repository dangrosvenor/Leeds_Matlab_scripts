%make dummy cbfA(i).c and cbfB(i).c

%cA contains the useful data cA(1).c(1,1) is the contour value (291 K here). N=cA(1).c(2,1) is then the number of x,y points that follow
%cA(1).c(1,2:N) will then be all the x values and (2,2:N) all the y ones


cb_dat = cbfA(i).c;

range_dat = maxcovOvr-mincovOvr;
%ncont is nnumber of contours required
dx_dat = sigfig(range_dat/np,1);
conts_dat = [mincovOvr:dx_dat:maxcovOvr];

for i_dat=1:length(conts_dat)
    cbfA_dummy(i).c(i_dat*2-1,1) = conts_dat(i_dat);
    cbfA_dummy(i).c(i_dat*2,1) = 1;
    cbfA_dummy(i).c(i_dat*2-1,2) = 0;
    cbfA_dummy(i).c(i_dat*2,2) = 0;
end

