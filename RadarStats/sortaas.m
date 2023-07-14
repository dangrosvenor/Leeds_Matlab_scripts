% sorts radar 5 in one file
iload=0;

if iload==1
	'loading rad stats...'
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\slices\sliceDat';    


end
clear aa3;

hdiff=hrs(2:end)-hrs(1:end-1);
fa=find(hdiff>2); %find where step up is large signifying the next slice - usually steps from 00hrs to 14hrs
fa(end+1)=length(aa);
fa(2:end+1)=fa(1:end);
fa(1)=0;

[x1 x2]=size(aa(1).i1);
for i=1:length(fa)-1
    for j=fa(i)+1:fa(i+1)
        if size(aa(j).i2,1)>0
            iz=find(aa(j).i2>0);
            aa3(i).i2max(j-fa(i))=min(min(aa(j).i2(iz))); %only intersested in 10dBz echo and y direction (i2)
        end    
    end                                
end




