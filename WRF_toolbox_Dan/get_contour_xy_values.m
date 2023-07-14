[cA.c cB.c]=contour(timesTH(i).t,height,pdat(i).p,[291 291],'r'); %x,y data for the 291 K potemp contour (run plot.. first
    %for vertical cross section
    
%cA contains the useful data cA(1).c(1,1) is the contour value (291 K here). N=cB(1).c(2,1) is then the number of x,y points that follow
%cA(1).c(1,2:N) will then be all the x values and (2,2:N) all the y ones

Y_291=cA(1).c(1,2:end);