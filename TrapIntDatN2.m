%finds, displays and stores the average, total integral and max of dat for column con

    if j==1;
        clear avall maxall totall;
    end;
    size(dat);
    npoints=ans(1)-1;
    totsum=0;
    negflag=0;
    posflag=0;
    for i=1:npoints;
        if dat(i,con)<0 
            negflag=1;
            %dat(i,con)=0; if don't want to count neg values
        end
        if dat(i,con)>0 
            posflag=1;
        end
        spacing=dat(i+1,1)-dat(i,1);
        totsum=totsum+0.5*spacing*(dat(i+1,con)+dat(i,con));
    end;
    
    %fprintf(1,'data%d tot=%E av.=%E\n',j,totsum,(totsum/dat(npoints,1)) );
    totall(j)=totsum;
    avall(j)=totsum/dat(npoints,1);
    if negflag==1 & posflag==0
        maxall(j)=min(dat(:,con));
    else
        maxall(j)=max(dat(:,con));
    end
   

