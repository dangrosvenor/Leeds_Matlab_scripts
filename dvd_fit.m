a=[1.892 1.988];
targ=4.488;

n=length(a);
d=0;
best=1e30;

    for i=1:n-2
        for j=i+1:n
            for k=j+1:n
                tot=sum(a([i j k]));
                if (tot-targ) < 0 & abs(tot-targ)<abs(best-targ)
                    best=tot;
                    besti=[i j k];
                end
            end
        end
    end
        best2=1e30;
        for j=1:n-1
            for k=j+1:n
                tot=sum(a([j k]));
                if (tot-targ) < 0 & abs(tot-targ)<abs(best2-targ)
                    best2=tot;
                    besti2=[j k];
                end
            end
        end
   
    best4=1e30;
   for i4=1:n-3 
    for i=i4+1:n
        for j=i+1:n
            for k=j+1:n
                tot=sum(a([i4 i j k]));
                if (tot-targ) < 0 & abs(tot-targ)<abs(best4-targ)
                    best4=tot;
                    besti4=[i4 i j k];
                end
            end
        end
    end
   end 