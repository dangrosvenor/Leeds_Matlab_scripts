function c=get_cbar_colours(im,map,i,ix,di,method)

%method='constant_spacing';
%method='look_for_new';

cstart=im(i,ix,:); %new colour


switch method
case 'constant_spacing'
    c(1,1,:)=cstart;
    icount=2;
    for iy=i:-di:1
        c(1,icount,:)=im(iy,ix,:); %new colour
        icount=icount+1;
    end
    
    
case 'look_for_new'==1

    %cnew(1,:,:)=[999 999 999];
    %cold(1,:,:)=[999 999 999];
    cold=cstart;
    cnew=cstart;
    
    blackf=0;
    start=1;
    
    j=1;
    b=1;
    c(1,1,:)=cstart; %store the first colour here
    %~all(a) means that cnew is not the same as cold
    %all(x) returns a one if all are non-zero so ~all returns one if there are any zero elements
    while (~all(a) & i>=2)|start==1 %until reach original colour (on assumption that 
        %after colorbar goes back to the orginal colour
        while (all(b) & i>=2)  %until colour changes
            cold=cnew;
            i=i-1;   %move up image
            cnew=im(i,ix,:); %new colour
            b=(cnew==cold);
        end
        j=j+1;
        c(1,j,:)=cnew;  %store new colours
        cold=cnew;
        start=0;
        a=(cnew==cstart);
        b=1; %reset b so we re-enter the inner while loop
    end %finds colour numbers of radar scale
    
end