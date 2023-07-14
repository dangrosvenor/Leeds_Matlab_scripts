function [conts,fcont,x1s,x2s]=follow_contour(cbfA,x1,x2,data,outcont)

if ~exist('outcont');
    outcont=0;
end

i=1;
jc=1; %counter for column position in cbfA.c matrix
icont=1; %counter for number of contours
while jc<size(cbfA(i).c,2)
    np=cbfA(i).c(2,jc);
    icont=icont+1;
    for ipoint=1:np %loop through all points in the contour
        conts(icont)=cbfA(i).c(1,jc);
        
        xval=cbfA(i).c(1,jc+ipoint);
        x2val=cbfA(i).c(2,jc+ipoint);
        
        if isnan(xval) | isnan(x2val)
            vcont(icont,ipoint)=NaN;
        else
            xfind=xval; %value on x1 axis of point on contour
            x2find=x2val; %value on x2 axis
            i1=findheight(x1,xfind); %index in x1 vector
            i2=findheight(x2,x2find); %x2 index
            
            fcont(icont,ipoint)=data(i1,i2); %store value of data at point on contour in vcont matrix
            
            if icont==outcont
                x1s(ipoint)=xfind;
                x2s(ipoint)=x2find;
            end
        end
        
    end

  
    
    jc=jc+cbfA(i).c(2,jc)+1;
end






