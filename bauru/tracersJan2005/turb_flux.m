%eddy flux from TwoD

th=TwoD.TH1 + repmat( GridDan(1).THREF , [1 size(TwoD.W,2)] );
thmean=mean(th,2);

thdash=th - repmat( thmean , [1 size(TwoD.W,2)] );

wdash=TwoD.W - repmat( mean(TwoD.W,2) , [1 size(TwoD.W,2)] );

for k=1:length(GridDan(1).Z)
    turb(k)=mean(TwoD.W(k,:).*TwoD.TH1(k,:)); %mean of w'*th'
end

turb2=mean( wdash .* thdash , 2);

wth=mean( TwoD.W .* th , 2); 



vdash= TwoD.V - repmat( mean(TwoD.V,2) , [1 size(TwoD.W,2)] );
vturb = mean( wdash .* vdash , 2);