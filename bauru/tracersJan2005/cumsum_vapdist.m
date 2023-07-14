bins=[0:0.1:30];

ih1=findheight(GridDan(idir).Z/1000+0.62,16);
ih2=findheight(GridDan(idir).Z/1000+0.62,19);

zs=GridDan(1).Z(ih1:ih2)/1000 + 0.62;
            
cmean=cumsum(bins(2:end).*squeeze(vapdist(1).v(1:end,13,1))') ...
    ./cumsum(squeeze(vapdist(1).v(1:end,13,1))');


clear cmean2

it1=1;
it2=51;


bins=(bins(2:end)+bins(1:end-1))/2;

bins2=repmat(bins(it1:it2),[size(vapdist(1).v,3) 1]);
for it=1:size(vapdist(1).v,2)
        cmean2(:,it)=sum(bins2.*squeeze(vapdist(1).v(it1:it2,it,:))' , 2) ...
    ./sum(squeeze(vapdist(1).v(it1:it2,it,:))' ,2);

    for k=1:size(vapdist(1).v,3)

        [a b]=max(vapdist(1).v(:,it,k),[],1);
        sum_bel(k,it)=-sum((bins(1:b-1)-bins(b)).*squeeze(vapdist(1).v(1:b-1,it,k))' );
        sum_abv(k,it)=sum((bins(b+1:end)-bins(b)).*squeeze(vapdist(1).v(b+1:end,it,k))' );
        
        meanvapd(k,it)=sum( bins(1:end).*squeeze(vapdist(1).v(1:end,it,k) )' )...
                            ./sum(vapdist(1).v(1:end,it,k) );
    end
        
end

a=find(vapdist(1).v>1e-4);
[b,c,d]=ind2sub(size(vapdist(1).v),a);

figure
plot(sum_abv(:,13),zs,'r');
hold on
plot(sum_bel(:,13),zs,'b');


    

