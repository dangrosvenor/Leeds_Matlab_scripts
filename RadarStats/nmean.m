function nme=nmean(np,vals)
    
        midvals=(vals(1:end-1)+vals(2:end))/2;
        totn=sum(np);
        if totn==0
            nme=0;
        else
            nme=sum(midvals.*np)/totn;
        end
    